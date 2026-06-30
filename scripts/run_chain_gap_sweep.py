#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONDA = Path("/home/senescence/miniconda3/bin/conda")
ENV = "vallescope_dev"
CONDA_BIN = Path("/home/senescence/miniconda3/envs/vallescope_dev/bin")
VSG_BENCH = CONDA_BIN / "vsg-benchmark"
SIM = Path(
    "/home/senescence/Documents/storage2/hiura/docker/bulk_env/data/00_analysis/"
    "Repthoscope.v1.0.0/vsg-benchmark/simulations/chr1_alphaSat_1M"
)
REFERENCE_TRUTH = SIM / "truth/reference/visor_reference.fa"
TRUTH_VCFS = [
    SIM / "truth/vcfs/sample01.truth.vcf",
    SIM / "truth/vcfs/sample02.truth.vcf",
    SIM / "truth/vcfs/sample03.truth.vcf",
]
INPUT_FASTAS = [
    ROOT / "results/assignment_sweep_20260630/input_fastas/reference.fa",
    ROOT / "results/assignment_sweep_20260630/input_fastas/sample01.fa",
    ROOT / "results/assignment_sweep_20260630/input_fastas/sample02.fa",
    ROOT / "results/assignment_sweep_20260630/input_fastas/sample03.fa",
]
GAPS = [500, 2000, 5000, 10000]


def run_logged(command: list[str], log_path: Path) -> tuple[float, int]:
    start = time.monotonic()
    with log_path.open("w") as log:
        log.write("$ " + " ".join(command) + "\n\n")
        log.flush()
        proc = subprocess.run(
            command,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return time.monotonic() - start, proc.returncode


def benchmark_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def micro(rows: list[dict[str, str]]) -> tuple[int, int, int, float, float, float]:
    tp = sum(int(row["tp_call"]) for row in rows)
    fp = sum(int(row["fp"]) for row in rows)
    fn = sum(int(row["fn"]) for row in rows)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def main() -> int:
    outdir = ROOT / "results/chain_gap_sweep_20260630"
    outdir.mkdir(parents=True, exist_ok=True)
    path_value = f"{ROOT / 'build'}:{CONDA_BIN}:{os.environ.get('PATH', '')}"

    rows_out: list[dict[str, object]] = []
    for gap in GAPS:
        run_name = f"vallescope2_chain_gap{gap}_20260630"
        graph = ROOT / "results/graphs" / f"{run_name}.vallescope2.gfa"
        bench_out = ROOT / "results/simulation" / f"{run_name}_graph_vcf"
        bench_summary = bench_out / "benchmark_summary.tsv"
        construct_log = outdir / f"gap{gap}.construct.log"
        benchmark_log = outdir / f"gap{gap}.benchmark.log"

        if not graph.exists():
            construct = [
                str(CONDA),
                "run",
                "-n",
                ENV,
                "env",
                f"PATH={path_value}",
                "PYTHONPATH=src",
                "vsg-graph",
                "construct",
                "--fasta",
                str(INPUT_FASTAS[0]),
                "--fasta",
                str(INPUT_FASTAS[1]),
                "--fasta",
                str(INPUT_FASTAS[2]),
                "--fasta",
                str(INPUT_FASTAS[3]),
                "--name",
                run_name,
                "--output",
                "results/graphs",
                "--run-directory",
                "runs",
                "--aligner",
                "vallescope2",
                "--threads",
                "8",
                "--aligner-extra",
                f"--max-chain-gap {gap}",
            ]
            sec, code = run_logged(construct, construct_log)
            if code != 0:
                print(f"gap={gap}: construct failed; see {construct_log}")
                continue
            print(f"gap={gap}: construct {sec:.1f}s")
        else:
            print(f"gap={gap}: reuse existing graph")

        if not bench_summary.exists():
            benchmark = [
                str(VSG_BENCH),
                "simulate",
                "benchmark-graph-vcf",
                "--gfa",
                str(graph),
                "--reference-path",
                "chr1_124500000-125500000",
                "--reference-fasta",
                str(REFERENCE_TRUTH),
                "--truth-vcf",
                str(TRUTH_VCFS[0]),
                "--truth-vcf",
                str(TRUTH_VCFS[1]),
                "--truth-vcf",
                str(TRUTH_VCFS[2]),
                "--output",
                str(bench_out),
                "--all-snarls",
                "--threads",
                "8",
                "--vg",
                str(CONDA_BIN / "vg"),
                "--bcftools",
                str(CONDA_BIN / "bcftools"),
                "--bgzip",
                str(CONDA_BIN / "bgzip"),
                "--tabix",
                str(CONDA_BIN / "tabix"),
                "--truvari",
                str(CONDA_BIN / "truvari"),
            ]
            sec, code = run_logged(benchmark, benchmark_log)
            if code != 0:
                print(f"gap={gap}: benchmark failed; see {benchmark_log}")
                continue
            print(f"gap={gap}: benchmark {sec:.1f}s")
        else:
            print(f"gap={gap}: reuse existing benchmark")

        rows = benchmark_rows(bench_summary)
        tp, fp, fn, precision, recall, f1 = micro(rows)
        paf = ROOT / "runs" / run_name / "all_to_all_paf/vallescope2/all.paf"
        paf_records = 0
        if paf.exists():
            with paf.open() as handle:
                paf_records = sum(1 for _ in handle)
        rows_out.append(
            {
                "max_chain_gap": gap,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": f"{precision:.6f}",
                "recall": f"{recall:.6f}",
                "f1": f"{f1:.6f}",
                "paf_records": paf_records,
            }
        )
        print(
            f"gap={gap}: F1={f1:.4f}, P={precision:.4f}, R={recall:.4f}, "
            f"TP/FP/FN={tp}/{fp}/{fn}, PAF={paf_records}"
        )

    out = outdir / "summary.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows_out[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
