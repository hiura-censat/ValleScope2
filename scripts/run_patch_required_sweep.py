#!/usr/bin/env python3
"""Temporary ValleScope2 parameter sweep with gap patching enabled."""

import csv
import os
import subprocess
import time
from pathlib import Path


ROOT = Path.cwd()
CONDA = "/home/senescence/miniconda3/bin/conda"
ENV_BIN = Path("/home/senescence/miniconda3/envs/vallescope_dev/bin")
TRUTH_ROOT = Path(
    "/home/senescence/Documents/storage2/hiura/docker/bulk_env/data/00_analysis/"
    "Repthoscope.v1.0.0/vsg-benchmark/simulations/chr1_alphaSat_1M/truth"
)
FASTAS = [
    "results/assignment_sweep_20260630/input_fastas/reference.fa",
    "results/assignment_sweep_20260630/input_fastas/sample01.fa",
    "results/assignment_sweep_20260630/input_fastas/sample02.fa",
    "results/assignment_sweep_20260630/input_fastas/sample03.fa",
]
TRUTH_VCFS = [
    TRUTH_ROOT / "vcfs/sample01.truth.vcf",
    TRUTH_ROOT / "vcfs/sample02.truth.vcf",
    TRUTH_ROOT / "vcfs/sample03.truth.vcf",
]


CONDITIONS = [
    # Current good chain settings, patching forced on.
    ("g500_r105_a20_b1_pg5000_w300_id85_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 5000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg10000_w300_id85_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 10000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg20000_w300_id85_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg70000_w300_id85_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg10000_w300_id90_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 10000 "
     "--patch-window-bp 300 --min-patch-identity 0.90 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg20000_w300_id90_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.90 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg70000_w300_id90_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.90 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg20000_w300_id95_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.95 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg70000_w300_id95_noce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.95 "
     "--max-bundle-align-bp 50000 --no-chain-extension --min-copy-support-anchors 1"),

    # Same, but chain extension on.
    ("g500_r105_a20_b1_pg10000_w300_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 10000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg20000_w300_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg70000_w300_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg20000_w300_id90_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.90 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg70000_w300_id90_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.90 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),

    # Larger patch windows.
    ("g500_r105_a20_b1_pg20000_w1000_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 1000 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg70000_w1000_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 1000 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg20000_w1000_id90_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 1000 --min-patch-identity 0.90 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b1_pg70000_w1000_id90_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 1000 --min-patch-identity 0.90 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),

    # Try slightly stricter chain filters while keeping patching on.
    ("g500_r105_a25_b1_pg20000_w300_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 25 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a25_b1_pg70000_w300_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 25 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b2_pg20000_w300_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 2 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g500_r105_a20_b2_pg70000_w300_id85_ce",
     "--max-chain-gap 500 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 2 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),

    # Slightly looser chain gap with patching.
    ("g750_r105_a20_b1_pg20000_w300_id85_ce",
     "--max-chain-gap 750 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
    ("g750_r105_a20_b1_pg70000_w300_id85_ce",
     "--max-chain-gap 750 --chain-max-gap-ratio 1.05 --min-chain-anchors 20 "
     "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 "
     "--patch-window-bp 300 --min-patch-identity 0.85 "
     "--max-bundle-align-bp 50000 --chain-extension --min-copy-support-anchors 1"),
]


def run(cmd, log):
    start = time.time()
    with open(log, "w") as handle:
        proc = subprocess.run(cmd, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    return proc.returncode, time.time() - start


def aggregate(summary_path):
    rows = list(csv.DictReader(open(summary_path), delimiter="\t"))
    tp = sum(int(r["tp_call"]) for r in rows)
    fp = sum(int(r["fp"]) for r in rows)
    fn = sum(int(r["fn"]) for r in rows)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def main():
    out_root = ROOT / "results/patch_required_sweep_20260704"
    log_root = out_root / "logs"
    log_root.mkdir(parents=True, exist_ok=True)
    summary = out_root / "summary.tsv"
    fields = [
        "idx", "name", "extra", "construct_seconds", "benchmark_seconds",
        "tp", "fp", "fn", "precision", "recall", "f1", "status",
    ]
    done = set()
    if summary.exists():
        for row in csv.DictReader(open(summary), delimiter="\t"):
            done.add(row["name"])
    write_header = not summary.exists()
    with open(summary, "a", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        if write_header:
            writer.writeheader()
        for idx, (name, extra) in enumerate(CONDITIONS):
            if name in done:
                continue
            print("RUN", idx, name, extra, flush=True)
            run_name = f"vs2_patchreq_{name}"
            path_value = f"{ROOT / 'build'}:{ENV_BIN}:{os.environ.get('PATH', '')}"
            construct_log = log_root / f"{name}.construct.log"
            bench_log = log_root / f"{name}.bench.log"
            construct_cmd = [
                CONDA, "run", "-n", "vallescope_dev", "env",
                f"PATH={path_value}",
                "PYTHONPATH=src",
                "vsg-graph", "construct",
            ]
            for fasta in FASTAS:
                construct_cmd.extend(["--fasta", fasta])
            construct_cmd.extend([
                "--name", run_name,
                "--output", "results/graphs",
                "--run-directory", "runs",
                "--aligner", "vallescope2",
                "--threads", "8",
                "--aligner-extra", extra,
            ])
            rc, construct_seconds = run(construct_cmd, construct_log)
            row = {
                "idx": idx,
                "name": name,
                "extra": extra,
                "construct_seconds": f"{construct_seconds:.2f}",
                "benchmark_seconds": ".",
                "tp": ".", "fp": ".", "fn": ".",
                "precision": ".", "recall": ".", "f1": ".",
                "status": "construct_failed" if rc else "ok",
            }
            if rc == 0:
                gfa = ROOT / "results/graphs" / f"{run_name}.vallescope2.gfa"
                bench_dir = ROOT / "results/simulation" / f"{run_name}_type_harmonized"
                bench_cmd = [
                    str(ENV_BIN / "vsg-benchmark"),
                    "simulate", "benchmark-graph-vcf",
                    "--gfa", str(gfa),
                    "--reference-path", "chr1_124500000-125500000",
                    "--reference-fasta", str(TRUTH_ROOT / "reference/visor_reference.fa"),
                ]
                for vcf in TRUTH_VCFS:
                    bench_cmd.extend(["--truth-vcf", str(vcf)])
                bench_cmd.extend([
                    "--output", str(bench_dir),
                    "--all-snarls",
                    "--threads", "8",
                    "--vg", str(ENV_BIN / "vg"),
                    "--bcftools", str(ENV_BIN / "bcftools"),
                    "--bgzip", str(ENV_BIN / "bgzip"),
                    "--tabix", str(ENV_BIN / "tabix"),
                    "--truvari", str(ENV_BIN / "truvari"),
                    "--type-ignore",
                    "--dup-to-ins",
                ])
                rc, benchmark_seconds = run(bench_cmd, bench_log)
                row["benchmark_seconds"] = f"{benchmark_seconds:.2f}"
                if rc == 0:
                    tp, fp, fn, precision, recall, f1 = aggregate(
                        bench_dir / "benchmark_summary.tsv")
                    row.update({
                        "tp": tp,
                        "fp": fp,
                        "fn": fn,
                        "precision": f"{precision:.6f}",
                        "recall": f"{recall:.6f}",
                        "f1": f"{f1:.6f}",
                        "status": "ok",
                    })
                    print(
                        "RESULT", idx, "P", f"{precision:.4f}",
                        "R", f"{recall:.4f}", "F1", f"{f1:.4f}",
                        "TP", tp, "FP", fp, "FN", fn,
                        flush=True)
                else:
                    row["status"] = "benchmark_failed"
            writer.writerow(row)
            handle.flush()


if __name__ == "__main__":
    main()
