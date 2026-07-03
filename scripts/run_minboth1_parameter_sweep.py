#!/usr/bin/env python3
"""Temporary parameter sweep for ValleScope2 min-chain-both-anchors experiments."""

import csv
import json
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
    ("baseline_minboth1", "--min-chain-both-anchors 1"),
    ("mcs10_pm10", "--min-chain-both-anchors 1 --min-candidate-score 10 --primary-margin 10"),
    ("mcs10_pm12", "--min-chain-both-anchors 1 --min-candidate-score 10 --primary-margin 12"),
    ("mcs15_pm10", "--min-chain-both-anchors 1 --min-candidate-score 15 --primary-margin 10"),
    ("mcs20_pm10", "--min-chain-both-anchors 1 --min-candidate-score 20 --primary-margin 10"),
    ("mca30_mcs10", "--min-chain-both-anchors 1 --min-chain-anchors 30 --min-candidate-score 10"),
    ("mca40_mcs10", "--min-chain-both-anchors 1 --min-chain-anchors 40 --min-candidate-score 10"),
    ("gap500_mcs10", "--min-chain-both-anchors 1 --max-chain-gap 500 --min-candidate-score 10"),
    ("gap500_mca30_mcs10", "--min-chain-both-anchors 1 --max-chain-gap 500 --min-chain-anchors 30 --min-candidate-score 10"),
    ("ratio105_mcs10", "--min-chain-both-anchors 1 --chain-max-gap-ratio 1.05 --min-candidate-score 10"),
    ("bundle25k_mcs10", "--min-chain-both-anchors 1 --max-bundle-align-bp 25000 --min-candidate-score 10"),
    ("patch20k_id90", "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 --min-patch-identity 0.90"),
    ("patch70k_id90", "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 --min-patch-identity 0.90"),
    ("patch20k_id95_mcs10", "--min-chain-both-anchors 1 --max-patch-gap-bp 20000 --min-patch-identity 0.95 --min-candidate-score 10"),
    ("patch70k_id95_mcs10", "--min-chain-both-anchors 1 --max-patch-gap-bp 70000 --min-patch-identity 0.95 --min-candidate-score 10"),
    ("strict_combo_patch20k", "--min-chain-both-anchors 1 --max-chain-gap 500 --min-chain-anchors 30 --min-candidate-score 10 --primary-margin 10 --max-patch-gap-bp 20000 --min-patch-identity 0.95"),
]


def run(cmd, log):
    start = time.time()
    with open(log, "w") as handle:
        proc = subprocess.run(cmd, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
    elapsed = time.time() - start
    return proc.returncode, elapsed


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
    out_root = ROOT / "results/minboth1_param_sweep_20260703"
    out_root.mkdir(parents=True, exist_ok=True)
    summary = out_root / "summary.tsv"
    fields = [
        "name", "extra", "construct_seconds", "benchmark_seconds",
        "tp", "fp", "fn", "precision", "recall", "f1",
        "gfa", "benchmark_dir", "construct_log", "benchmark_log",
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
        for name, extra in CONDITIONS:
            if name in done:
                continue
            run_name = f"vallescope2_minboth1_sweep_{name}_20260703"
            construct_log = out_root / f"{name}.construct.log"
            benchmark_log = out_root / f"{name}.benchmark.log"
            path_value = f"{ROOT / 'build'}:{ENV_BIN}:{os.environ.get('PATH', '')}"
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
            if rc != 0:
                writer.writerow({
                    "name": name, "extra": extra,
                    "construct_seconds": f"{construct_seconds:.3f}",
                    "benchmark_seconds": ".", "tp": ".", "fp": ".", "fn": ".",
                    "precision": ".", "recall": ".", "f1": ".",
                    "gfa": ".", "benchmark_dir": ".",
                    "construct_log": construct_log, "benchmark_log": benchmark_log,
                })
                handle.flush()
                continue

            gfa = ROOT / "results/graphs" / f"{run_name}.vallescope2.gfa"
            benchmark_dir = ROOT / "results/simulation" / f"{run_name}_type_harmonized"
            benchmark_cmd = [
                str(ENV_BIN / "vsg-benchmark"),
                "simulate", "benchmark-graph-vcf",
                "--gfa", str(gfa),
                "--reference-path", "chr1_124500000-125500000",
                "--reference-fasta", str(TRUTH_ROOT / "reference/visor_reference.fa"),
            ]
            for vcf in TRUTH_VCFS:
                benchmark_cmd.extend(["--truth-vcf", str(vcf)])
            benchmark_cmd.extend([
                "--output", str(benchmark_dir),
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
            rc, benchmark_seconds = run(benchmark_cmd, benchmark_log)
            if rc == 0:
                tp, fp, fn, precision, recall, f1 = aggregate(
                    benchmark_dir / "benchmark_summary.tsv")
            else:
                tp = fp = fn = "."
                precision = recall = f1 = "."
            writer.writerow({
                "name": name,
                "extra": extra,
                "construct_seconds": f"{construct_seconds:.3f}",
                "benchmark_seconds": f"{benchmark_seconds:.3f}",
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": f"{precision:.6f}" if isinstance(precision, float) else precision,
                "recall": f"{recall:.6f}" if isinstance(recall, float) else recall,
                "f1": f"{f1:.6f}" if isinstance(f1, float) else f1,
                "gfa": gfa,
                "benchmark_dir": benchmark_dir,
                "construct_log": construct_log,
                "benchmark_log": benchmark_log,
            })
            handle.flush()
            print(name, precision, recall, f1)


if __name__ == "__main__":
    main()
