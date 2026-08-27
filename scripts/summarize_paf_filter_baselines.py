#!/usr/bin/env python3
"""Summarize PAF-filter, seqwish-graph, snarl, and SV callset baselines."""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def merge_bp(intervals):
    total = 0
    end = -1
    for start, stop in sorted(intervals):
        if start > end:
            total += stop - start
            end = stop
        elif stop > end:
            total += stop - end
            end = stop
    return total


def read_paf(path):
    records = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            records.append(
                (
                    fields[0],
                    int(fields[2]),
                    int(fields[3]),
                    fields[5],
                    int(fields[7]),
                    int(fields[8]),
                    fields[4],
                )
            )
    return records


def paf_key(record):
    return record


def duplication_candidates(records, overlap_fraction=0.8, distance=10_000):
    """Return records participating in strong query overlap at distinct loci."""
    by_query = defaultdict(list)
    for record in records:
        by_query[record[0]].append(record)
    candidates = set()
    for query_records in by_query.values():
        for index, left in enumerate(query_records):
            left_len = left[2] - left[1]
            for right in query_records[index + 1 :]:
                overlap = min(left[2], right[2]) - max(left[1], right[1])
                if overlap <= 0 or overlap / min(left_len, right[2] - right[1]) < overlap_fraction:
                    continue
                distinct = left[3] != right[3] or abs(left[4] - right[4]) >= distance
                if distinct:
                    candidates.add(paf_key(left))
                    candidates.add(paf_key(right))
    return candidates


def paf_metrics(records, reference, original_duplication_candidates):
    query_intervals = defaultdict(list)
    target_intervals = defaultdict(list)
    reference_intervals = []
    for query, q_start, q_end, target, t_start, t_end, _strand in records:
        query_intervals[query].append((q_start, q_end))
        target_intervals[target].append((t_start, t_end))
        if query == reference:
            reference_intervals.append((q_start, q_end))
        if target == reference:
            reference_intervals.append((t_start, t_end))
    retained = {paf_key(record) for record in records}
    return {
        "paf_records": len(records),
        "paf_query_span_sum": sum(record[2] - record[1] for record in records),
        "paf_target_span_sum": sum(record[5] - record[4] for record in records),
        "paf_query_union_bp": sum(merge_bp(value) for value in query_intervals.values()),
        "paf_target_union_bp": sum(merge_bp(value) for value in target_intervals.values()),
        "reference_direct_union_bp": merge_bp(reference_intervals),
        "dup_candidate_records_retained": len(original_duplication_candidates & retained),
    }


def graph_metrics(path):
    metrics = {"graph_nodes": 0, "graph_edges": 0, "graph_paths": 0, "graph_node_bp": 0}
    with path.open() as handle:
        for line in handle:
            if line.startswith("S\t"):
                metrics["graph_nodes"] += 1
                metrics["graph_node_bp"] += len(line.split("\t", 3)[2])
            elif line.startswith("L\t"):
                metrics["graph_edges"] += 1
            elif line.startswith("P\t"):
                metrics["graph_paths"] += 1
    return metrics


def snarl_metrics(path):
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    return {
        "snarls": len(rows),
        "max_snarl_shallow_bp": max((int(row["Shallow-bases"]) for row in rows), default=0),
        "max_snarl_deep_bp": max((int(row["Deep-Bases"]) for row in rows), default=0),
        "max_snarl_deep_nodes": max((int(row["Deep-Nodes"]) for row in rows), default=0),
    }


def vcf_metrics(path):
    opener = gzip.open if path.suffix == ".gz" else open
    counts = {"sv50": 0, "sv50_ins": 0, "sv50_del": 0, "max_sv_bp": 0}
    ac_values = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            info = {
                key: value
                for item in fields[7].split(";")
                if "=" in item
                for key, value in [item.split("=", 1)]
            }
            size = abs(len(fields[3]) - len(fields[4]))
            kind = "INS" if len(fields[4]) > len(fields[3]) else "DEL"
            counts["sv50"] += 1
            counts[f"sv50_{kind.lower()}"] += 1
            counts["max_sv_bp"] = max(counts["max_sv_bp"], size)
            ac_values.append(int(info.get("AC", "0").split(",")[0]))
    counts["ac1"] = sum(value == 1 for value in ac_values)
    counts["ac2_3"] = sum(2 <= value <= 3 for value in ac_values)
    counts["ac4_plus"] = sum(value >= 4 for value in ac_values)
    return counts


def parse_baseline(value):
    fields = value.split("=", 1)
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("expected LABEL=DIR")
    return fields[0], Path(fields[1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-paf", type=Path, required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--baseline", action="append", type=parse_baseline, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    original = read_paf(args.original_paf)
    duplication = duplication_candidates(original)
    rows = []
    for label, directory in args.baseline:
        row = {"baseline": label}
        row.update(paf_metrics(read_paf(directory / "filtered.paf"), args.reference, duplication))
        graph = directory / "graph.k0.gfa"
        if not graph.exists():
            graph = directory / "graph.gfa"
        row.update(graph_metrics(graph))
        row.update(snarl_metrics(directory / "snarl_stats.tsv"))
        row.update(vcf_metrics(directory / "deconstruct.sv50.ac_positive.vcf.gz"))
        rows.append(row)

    for row in rows:
        row["paf_record_retention"] = row["paf_records"] / len(original)
        row["dup_candidate_retention"] = (
            row["dup_candidate_records_retained"] / len(duplication) if duplication else 0.0
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
