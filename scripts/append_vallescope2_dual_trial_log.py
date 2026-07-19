#!/usr/bin/env python3
"""Append one completed dual-dataset trial to the shared comparison table."""

import argparse
import csv
import datetime as dt
import fcntl
from pathlib import Path


ALPHA_SAMPLES = ["HG01114_hap1", "NA18534_hap2", "NA18974_hap2"]
SIMULATION_SAMPLES = ["sample01", "sample02", "sample03"]
ALPHA_METRICS = [
    "query_coverage", "reference_coverage", "query_gap_count",
    "query_gap_count_ge_1kb", "query_gap_count_ge_10kb",
    "reference_gap_count", "reference_gap_count_ge_1kb",
    "reference_gap_count_ge_10kb", "query_largest_uncovered_bp",
    "reference_largest_uncovered_bp",
]


def read_key_values(path):
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["metric"]: row["value"] for row in rows}


def field_names():
    fields = [
        "timestamp", "output_name", "git_commit", "options",
        "alpha_construct_seconds", "simulation_total_seconds",
    ]
    for sample in ALPHA_SAMPLES:
        fields.extend(f"alpha_{sample}_{metric}" for metric in ALPHA_METRICS)
    fields.extend(f"simulation_{sample}_f1" for sample in SIMULATION_SAMPLES)
    fields.extend([
        "simulation_tp", "simulation_fp", "simulation_fn",
        "simulation_precision", "simulation_recall", "simulation_f1",
    ])
    fields.extend(f"raw_seqwish_{sample}_f1" for sample in SIMULATION_SAMPLES)
    fields.extend([
        "raw_seqwish_tp", "raw_seqwish_fp", "raw_seqwish_fn",
        "raw_seqwish_precision", "raw_seqwish_recall", "raw_seqwish_f1",
    ])
    return fields


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trial-dir", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()

    summary_path = args.trial_dir / "trial_summary.tsv"
    with summary_path.open() as handle:
        summary = list(csv.DictReader(handle, delimiter="\t"))
    indexed = {(row["dataset"], row["sample"]): row for row in summary}
    metadata = read_key_values(args.trial_dir / "run_metadata.tsv")

    row = {field: "." for field in field_names()}
    row.update({
        "timestamp": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "output_name": args.trial_dir.name,
        "git_commit": (args.trial_dir / "git_commit.txt").read_text().strip(),
        "options": metadata["options"],
        "alpha_construct_seconds": metadata["alpha_construct_seconds"],
        "simulation_total_seconds": metadata["simulation_total_seconds"],
    })
    for sample in ALPHA_SAMPLES:
        result = indexed[("alpha_sat", sample)]
        for metric in ALPHA_METRICS:
            row[f"alpha_{sample}_{metric}"] = result[metric]
    for sample in SIMULATION_SAMPLES:
        row[f"simulation_{sample}_f1"] = indexed[("simulation", sample)]["f1"]
    aggregate = indexed[("simulation", "aggregate")]
    for metric in ("tp", "fp", "fn", "precision", "recall", "f1"):
        row[f"simulation_{metric}"] = aggregate[metric]
    for sample in SIMULATION_SAMPLES:
        row[f"raw_seqwish_{sample}_f1"] = indexed[
            ("simulation_raw_seqwish", sample)]["f1"]
    raw_aggregate = indexed[("simulation_raw_seqwish", "aggregate")]
    for metric in ("tp", "fp", "fn", "precision", "recall", "f1"):
        row[f"raw_seqwish_{metric}"] = raw_aggregate[metric]

    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.touch(exist_ok=True)
    with args.log.open("r+", newline="") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter="\t")
        existing_header = reader.fieldnames or []
        existing_rows = list(reader)
        deprecated_fields = {"output_dir"}
        unknown_fields = set(existing_header) - set(field_names()) - deprecated_fields
        if unknown_fields:
            raise RuntimeError(
                f"existing sweep log has unknown columns: {sorted(unknown_fields)}")
        existing_rows = [
            {field: existing.get(field, ".") or "." for field in field_names()}
            for existing in existing_rows
        ]
        replaced = False
        for index, existing in enumerate(existing_rows):
            if existing["output_name"] != row["output_name"]:
                continue
            row["timestamp"] = existing["timestamp"]
            existing_rows[index] = row
            replaced = True
            break
        if not replaced:
            existing_rows.append(row)
        handle.seek(0)
        handle.truncate()
        writer = csv.DictWriter(
            handle, fieldnames=field_names(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(existing_rows)
        handle.flush()
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


if __name__ == "__main__":
    main()
