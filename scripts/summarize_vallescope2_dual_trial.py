#!/usr/bin/env python3
"""Combine alpha-satellite coverage/gaps and simulation benchmark metrics."""

import argparse
import csv
from pathlib import Path


FIELDS = [
    "dataset", "sample", "blocks", "query_coverage", "reference_coverage",
    "query_largest_uncovered_bp", "reference_largest_uncovered_bp",
    "query_gap_count", "query_gap_count_ge_1kb", "query_gap_count_ge_10kb",
    "reference_gap_count", "reference_gap_count_ge_1kb",
    "reference_gap_count_ge_10kb", "tp", "fp", "fn", "precision", "recall",
    "f1",
]


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def coverage_values(row):
    return {
        "blocks": row["blocks"],
        "query_coverage": row["query_coverage"],
        "reference_coverage": row["reference_coverage"],
        "query_largest_uncovered_bp": row["query_largest_uncovered_bp"],
        "reference_largest_uncovered_bp": row["reference_largest_uncovered_bp"],
        "query_gap_count": row["query_gap_count"],
        "query_gap_count_ge_1kb": row["query_gap_count_ge_1kb"],
        "query_gap_count_ge_10kb": row["query_gap_count_ge_10kb"],
        "reference_gap_count": row["reference_gap_count"],
        "reference_gap_count_ge_1kb": row["reference_gap_count_ge_1kb"],
        "reference_gap_count_ge_10kb": row["reference_gap_count_ge_10kb"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha-coverage", type=Path, required=True)
    parser.add_argument("--simulation-coverage", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--raw-benchmark", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    for coverage in read_tsv(args.alpha_coverage):
        rows.append({
            **{field: "." for field in FIELDS},
            "dataset": "alpha_sat", "sample": coverage["sample"],
            **coverage_values(coverage),
        })

    simulation_coverage = {
        row["sample"]: row for row in read_tsv(args.simulation_coverage)
    }
    for dataset, benchmark_path in (
        ("simulation", args.benchmark),
        ("simulation_raw_seqwish", args.raw_benchmark),
    ):
        benchmark = read_tsv(benchmark_path)
        tp = fp = fn = 0
        for result in benchmark:
            sample = result["sample"]
            coverage = simulation_coverage.get(sample)
            row = {field: "." for field in FIELDS}
            row.update({
                "dataset": dataset, "sample": sample,
                "tp": result["tp_call"], "fp": result["fp"], "fn": result["fn"],
                "precision": result["precision"], "recall": result["recall"],
                "f1": result["f1"],
            })
            if coverage:
                row.update(coverage_values(coverage))
            rows.append(row)
            tp += int(result["tp_call"])
            fp += int(result["fp"])
            fn += int(result["fn"])

        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        aggregate = {field: "." for field in FIELDS}
        aggregate.update({
            "dataset": dataset, "sample": "aggregate", "tp": tp, "fp": fp,
            "fn": fn, "precision": f"{precision:.6f}",
            "recall": f"{recall:.6f}", "f1": f"{f1:.6f}",
        })
        rows.append(aggregate)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
