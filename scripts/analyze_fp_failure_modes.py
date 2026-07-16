#!/usr/bin/env python3
"""Classify Truvari false positives by graph representation failure mode."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Record:
    sample: str
    pos: int
    record_id: str
    ref: str
    alt: str
    info: dict[str, str]

    @property
    def start(self) -> int:
        return self.pos - 1

    @property
    def end(self) -> int:
        return self.start + len(self.ref)


def parse_info(value: str) -> dict[str, str]:
    result = {}
    for field in value.split(";"):
        if "=" in field:
            key, item = field.split("=", 1)
            result[key] = item
        elif field:
            result[field] = "1"
    return result


def read_vcf(path: Path, sample: str) -> list[Record]:
    records = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            records.append(
                Record(
                    sample=sample,
                    pos=int(fields[1]),
                    record_id=fields[2],
                    ref=fields[3],
                    alt=fields[4],
                    info=parse_info(fields[7]),
                )
            )
    return records


def comparison_match_id(record: Record) -> str:
    values = record.info.get("MatchId", "").split(",")
    return values[1] if len(values) > 1 else ""


def as_float(record: Record, key: str) -> float | None:
    value = record.info.get(key)
    return float(value) if value not in (None, "", ".") else None


def graph_kind(record: Record) -> str:
    delta = len(record.alt) - len(record.ref)
    if len(record.ref) == 1 and delta > 0:
        return "INS"
    if len(record.alt) == 1 and delta < 0:
        return "DEL"
    if delta == 0:
        return "COMPLEX_EQUAL_LENGTH"
    return "COMPLEX_NET_INS" if delta > 0 else "COMPLEX_NET_DEL"


def truth_interval(record: Record, truth_map: dict[str, tuple[int, int]]) -> tuple[int, int]:
    if record.record_id in truth_map:
        return truth_map[record.record_id]
    end = int(record.info.get("END", record.end))
    return record.start, max(record.start + 1, end)


def interval_distance(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def classify(
    fp: Record,
    direct_fns: list[Record],
    nearby_fns: list[Record],
    nearby_tps: list[Record],
) -> str:
    linked = {record.record_id: record for record in direct_fns + nearby_fns}
    linked_types = {record.info.get("SVTYPE", "") for record in linked.values()}
    if "DUP" in linked_types and len(linked) >= 2:
        return "dup_collapsed_with_adjacent_event"
    if len(direct_fns) >= 2 or len(linked) >= 2:
        return "multi_event_bubble"
    if linked_types == {"INV"}:
        return "inversion_representation_fragmentation"
    if linked_types == {"DUP"}:
        return "duplication_representation_mismatch"
    if direct_fns:
        size_similarity = as_float(fp, "PctSizeSimilarity")
        sequence_similarity = as_float(fp, "PctSeqSimilarity")
        start_distance = abs(int(fp.info.get("StartDistance", "0")))
        end_distance = abs(int(fp.info.get("EndDistance", "0")))
        if size_similarity is not None and size_similarity < 0.7:
            return "single_event_size_mismatch"
        if sequence_similarity is not None and sequence_similarity < 0.7:
            return "single_event_sequence_mismatch"
        if max(start_distance, end_distance) > 500:
            return "single_event_breakpoint_mismatch"
        return "single_event_threshold_mismatch"
    if nearby_tps:
        return "redundant_call_near_matched_event"
    graph_size = max(abs(len(fp.alt) - len(fp.ref)), len(fp.ref), len(fp.alt))
    if int(fp.info.get("LV", "0")) > 0:
        return "unmatched_nested_snarl_call"
    if graph_size <= 500:
        return "unmatched_small_call"
    return "unmatched_graph_call"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", type=Path)
    parser.add_argument("--truth-map", type=Path)
    parser.add_argument("--raw-benchmark", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nearby-bp", type=int, default=2500)
    args = parser.parse_args()

    truth_map: dict[str, tuple[int, int]] = {}
    if args.truth_map:
        with args.truth_map.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                truth_map[row["id"]] = (int(row["ref_start"]), int(row["ref_end"]))

    raw_fp_keys = set()
    if args.raw_benchmark:
        for sample_dir in sorted((args.raw_benchmark / "samples").glob("sample*")):
            sample = sample_dir.name
            for fp in read_vcf(sample_dir / "truvari" / "fp.vcf.gz", sample):
                raw_fp_keys.add((sample, fp.pos, fp.ref, fp.alt))

    rows = []
    category_counts = Counter()
    stage_counts = Counter()
    sample_category_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for sample_dir in sorted((args.benchmark / "samples").glob("sample*")):
        sample = sample_dir.name
        truvari = sample_dir / "truvari"
        fps = read_vcf(truvari / "fp.vcf.gz", sample)
        fns = read_vcf(truvari / "fn.vcf.gz", sample)
        tps = read_vcf(truvari / "tp-base.vcf.gz", sample)
        fns_by_match: dict[str, list[Record]] = defaultdict(list)
        for fn in fns:
            match_id = comparison_match_id(fn)
            if match_id:
                fns_by_match[match_id].append(fn)

        for fp in fps:
            match_id = comparison_match_id(fp)
            direct_fns = fns_by_match.get(match_id, []) if match_id else []
            nearby_fns = []
            nearby_tps = []
            for truth, destination in ((fns, nearby_fns), (tps, nearby_tps)):
                for record in truth:
                    start, end = truth_interval(record, truth_map)
                    distance = interval_distance(fp.start, fp.end, start, end)
                    if distance <= args.nearby_bp:
                        destination.append(record)

            category = classify(fp, direct_fns, nearby_fns, nearby_tps)
            stage_origin = "not_compared"
            if args.raw_benchmark:
                key = (sample, fp.pos, fp.ref, fp.alt)
                stage_origin = (
                    "present_in_raw_seqwish"
                    if key in raw_fp_keys
                    else "introduced_after_seqwish_smoothing"
                )
            category_counts[category] += 1
            stage_counts[stage_origin] += 1
            sample_category_counts[sample][category] += 1
            linked_fns = {r.record_id: r for r in direct_fns + nearby_fns}
            rows.append(
                {
                    "sample": sample,
                    "position": fp.pos,
                    "graph_id": fp.record_id,
                    "graph_kind": graph_kind(fp),
                    "ref_length": len(fp.ref),
                    "alt_length": len(fp.alt),
                    "net_size": len(fp.alt) - len(fp.ref),
                    "snarl_level": fp.info.get("LV", "0"),
                    "parent_snarl": fp.info.get("PS", ""),
                    "match_id": match_id,
                    "pct_sequence": fp.info.get("PctSeqSimilarity", ""),
                    "pct_size": fp.info.get("PctSizeSimilarity", ""),
                    "start_distance": fp.info.get("StartDistance", ""),
                    "end_distance": fp.info.get("EndDistance", ""),
                    "direct_fn_ids": ",".join(r.record_id for r in direct_fns),
                    "nearby_fn_ids": ",".join(linked_fns),
                    "nearby_tp_ids": ",".join(r.record_id for r in nearby_tps),
                    "stage_origin": stage_origin,
                    "category": category,
                }
            )

    args.output.mkdir(parents=True, exist_ok=True)
    detail_path = args.output / "fp_failure_modes.tsv"
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary_path = args.output / "fp_failure_mode_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "sample01", "sample02", "sample03", "total"])
        for category, total in category_counts.most_common():
            writer.writerow(
                [category]
                + [sample_category_counts[s][category] for s in ("sample01", "sample02", "sample03")]
                + [total]
            )

    stage_path = args.output / "fp_stage_summary.tsv"
    with stage_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["stage_origin", "count"])
        writer.writerows(stage_counts.most_common())

    print(f"FP records: {len(rows)}")
    for category, count in category_counts.most_common():
        print(f"{category}\t{count}")
    print(detail_path)
    print(summary_path)
    print(stage_path)


if __name__ == "__main__":
    main()
