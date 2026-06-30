#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONDA = Path("/home/senescence/miniconda3/bin/conda")
ENV = "vallescope_dev"
VSG_BENCH = Path("/home/senescence/miniconda3/envs/vallescope_dev/bin/vsg-benchmark")
SIM = Path(
    "/home/senescence/Documents/storage2/hiura/docker/bulk_env/data/00_analysis/"
    "Repthoscope.v1.0.0/vsg-benchmark/simulations/chr1_alphaSat_1M"
)

REFERENCE = SIM / "truth/reference/visor_reference.fa"
SAMPLES = [
    SIM / "haplotypes_visor_run2/sample01/h1.fa",
    SIM / "haplotypes_visor_run2/sample02/h1.fa",
    SIM / "haplotypes_visor_run2/sample03/h1.fa",
]
PREPARED_FASTAS = [
    ("reference", REFERENCE, "chr1_124500000-125500000"),
    ("sample01", SAMPLES[0], "sample01"),
    ("sample02", SAMPLES[1], "sample02"),
    ("sample03", SAMPLES[2], "sample03"),
]
TRUTH_VCFS = [
    SIM / "truth/vcfs/sample01.truth.vcf",
    SIM / "truth/vcfs/sample02.truth.vcf",
    SIM / "truth/vcfs/sample03.truth.vcf",
]


CONDITIONS = [
    ("baseline", ""),
    ("shared1", "--min-shared-tmus 1"),
    ("shared5", "--min-shared-tmus 5"),
    ("score_m30", "-mcs -30"),
    ("pm1", "-pm 1"),
]


def run_logged(command: list[str], log_path: Path, env: dict[str, str]) -> tuple[float, int]:
    start = time.monotonic()
    with log_path.open("w") as log:
        log.write("$ " + " ".join(command) + "\n\n")
        log.flush()
        proc = subprocess.run(
            command,
            cwd=ROOT,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    elapsed = time.monotonic() - start
    return elapsed, proc.returncode


def parse_benchmark_summary(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open() as handle:
        return json.load(handle)


def prepare_named_fastas(outdir: Path) -> list[Path]:
    input_dir = outdir / "input_fastas"
    input_dir.mkdir(parents=True, exist_ok=True)
    prepared: list[Path] = []
    for stem, source, record_name in PREPARED_FASTAS:
        destination = input_dir / f"{stem}.fa"
        with source.open() as src, destination.open("w") as dst:
            first_header = True
            for line in src:
                if first_header and line.startswith(">"):
                    dst.write(f">{record_name}\n")
                    first_header = False
                else:
                    dst.write(line)
        source_index = Path(f"{source}.fai")
        if source_index.exists():
            destination_index = Path(f"{destination}.fai")
            if destination_index.exists():
                destination_index.unlink()
        prepared.append(destination)
    return prepared


def micro_scores(rows: list[dict[str, str]]) -> tuple[int, int, int, float, float, float]:
    tp = sum(int(row["tp_call"]) for row in rows)
    fp = sum(int(row["fp"]) for row in rows)
    fn = sum(int(row["fn"]) for row in rows)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def main() -> int:
    outdir = ROOT / "results/assignment_sweep_20260630"
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "summary.tsv"
    input_fastas = prepare_named_fastas(outdir)

    env = os.environ.copy()
    conda_bin = Path("/home/senescence/miniconda3/envs/vallescope_dev/bin")
    env["PATH"] = f"{ROOT / 'build'}:{conda_bin}:{env.get('PATH', '')}"
    env["PYTHONPATH"] = "src"

    rows_out: list[dict[str, object]] = []
    for label, extra in CONDITIONS:
        run_name = f"vallescope2_assignment_{label}_20260630"
        graph = ROOT / "results/graphs" / f"{run_name}.vallescope2.gfa"
        bench_out = ROOT / "results/simulation" / f"{run_name}_graph_vcf"
        construct_log = outdir / f"{label}.construct.log"
        benchmark_log = outdir / f"{label}.benchmark.log"

        aligner_extra = extra.strip()
        construct = [
            str(CONDA),
            "run",
            "-n",
            ENV,
            "env",
            f"PATH={env['PATH']}",
            "PYTHONPATH=src",
            "vsg-graph",
            "construct",
            "--fasta",
            str(input_fastas[0]),
            "--fasta",
            str(input_fastas[1]),
            "--fasta",
            str(input_fastas[2]),
            "--fasta",
            str(input_fastas[3]),
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
        ]
        if aligner_extra:
            construct += ["--aligner-extra", aligner_extra]

        construct_sec, construct_code = run_logged(construct, construct_log, env)
        if construct_code != 0:
            print(f"{label}: construct failed; see {construct_log}", file=sys.stderr)
            continue

        benchmark = [
            str(VSG_BENCH),
            "simulate",
            "benchmark-graph-vcf",
            "--gfa",
            str(graph),
            "--reference-path",
            "chr1_124500000-125500000",
            "--reference-fasta",
            str(REFERENCE),
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
            "/home/senescence/miniconda3/envs/vallescope_dev/bin/vg",
            "--bcftools",
            "/home/senescence/miniconda3/envs/vallescope_dev/bin/bcftools",
            "--bgzip",
            "/home/senescence/miniconda3/envs/vallescope_dev/bin/bgzip",
            "--tabix",
            "/home/senescence/miniconda3/envs/vallescope_dev/bin/tabix",
            "--truvari",
            "/home/senescence/miniconda3/envs/vallescope_dev/bin/truvari",
        ]
        benchmark_sec, benchmark_code = run_logged(benchmark, benchmark_log, env)
        if benchmark_code != 0:
            print(f"{label}: benchmark failed; see {benchmark_log}", file=sys.stderr)
            continue

        bench_rows = parse_benchmark_summary(bench_out / "benchmark_summary.tsv")
        tp, fp, fn, micro_p, micro_r, micro_f1 = micro_scores(bench_rows)
        paf = ROOT / "runs" / run_name / "all_to_all_paf/vallescope2/all.paf"
        paf_records = 0
        if paf.exists():
            with paf.open() as handle:
                paf_records = sum(1 for _ in handle)
        base_meta = read_json(ROOT / "runs" / run_name / "all_to_all_paf/vallescope2/base_alignment.meta.json")
        assignment_meta = read_json(ROOT / "runs" / run_name / "all_to_all_paf/vallescope2/assignments.meta.json")

        for row in bench_rows:
            rows_out.append(
                {
                    "condition": label,
                    "assignment_params": aligner_extra or "default",
                    "construct_sec": f"{construct_sec:.2f}",
                    "benchmark_sec": f"{benchmark_sec:.2f}",
                    "paf_records": paf_records,
                    "anchor_guided": base_meta.get("anchor_guided_alignment_count", "."),
                    "whole_bundle": base_meta.get("whole_bundle_alignment_count", "."),
                    "assignment_primary": assignment_meta.get("primary_count", "."),
                    "assignment_ambiguous_query": assignment_meta.get("ambiguous_query_count", "."),
                    "assignment_unmatched": assignment_meta.get("unmatched_count", "."),
                    "micro_tp": tp,
                    "micro_fp": fp,
                    "micro_fn": fn,
                    "micro_precision": f"{micro_p:.6f}",
                    "micro_recall": f"{micro_r:.6f}",
                    "micro_f1": f"{micro_f1:.6f}",
                    **row,
                }
            )
        print(
            f"{label}: micro F1={micro_f1:.4f}, P={micro_p:.4f}, R={micro_r:.4f}, "
            f"TP/FP/FN={tp}/{fp}/{fn}, construct={construct_sec:.1f}s, benchmark={benchmark_sec:.1f}s"
        )

    fieldnames = [
        "condition",
        "assignment_params",
        "construct_sec",
        "benchmark_sec",
        "paf_records",
        "anchor_guided",
        "whole_bundle",
        "assignment_primary",
        "assignment_ambiguous_query",
        "assignment_unmatched",
        "micro_tp",
        "micro_fp",
        "micro_fn",
        "micro_precision",
        "micro_recall",
        "micro_f1",
    ]
    for row in rows_out:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
