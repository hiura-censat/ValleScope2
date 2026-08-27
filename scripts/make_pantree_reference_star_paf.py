#!/usr/bin/env python3
"""Create a reference-star PAF with non-overlapping sample intervals."""

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List


def alignment_score(fields: List[str]) -> int:
    for tag in fields[12:]:
        if tag.startswith("as:i:"):
            return int(tag[5:])
    return int(fields[9])


def swap_cigar(tag: str) -> str:
    if not tag.startswith("cg:Z:"):
        return tag
    return "cg:Z:" + re.sub(
        "[ID]", lambda match: "D" if match.group() == "I" else "I", tag[5:]
    )


def reference_first(fields: List[str], reference: str) -> List[str]:
    if fields[0] == reference:
        return fields
    if fields[5] != reference:
        return []
    return (
        fields[5:9]
        + [fields[4]]
        + fields[0:4]
        + [fields[9], fields[10], fields[11]]
        + [swap_cigar(tag) for tag in fields[12:]]
    )


def overlaps(start: int, end: int, intervals: List[List[int]]) -> bool:
    return any(start < kept_end and kept_start < end for kept_start, kept_end in intervals)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-paf", required=True, type=Path)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output-paf", required=True, type=Path)
    parser.add_argument("--summary-tsv", type=Path)
    args = parser.parse_args()

    by_sample: Dict[str, List[List[str]]] = defaultdict(list)
    input_records = 0
    with args.input_paf.open() as handle:
        for line in handle:
            input_records += 1
            fields = reference_first(line.rstrip().split("\t"), args.reference)
            if fields:
                by_sample[fields[5]].append(fields)

    selected = []
    for sample, records in by_sample.items():
        kept_intervals: List[List[int]] = []
        ordered = sorted(
            records,
            key=lambda fields: (
                alignment_score(fields),
                int(fields[8]) - int(fields[7]),
                int(fields[9]),
            ),
            reverse=True,
        )
        for fields in ordered:
            start, end = int(fields[7]), int(fields[8])
            if overlaps(start, end, kept_intervals):
                continue
            kept_intervals.append([start, end])
            selected.append(fields)

    selected.sort(key=lambda fields: (fields[5], int(fields[7]), int(fields[8])))
    args.output_paf.parent.mkdir(parents=True, exist_ok=True)
    with args.output_paf.open("w") as output:
        for fields in selected:
            output.write("\t".join(fields) + "\n")

    if args.summary_tsv:
        args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_tsv.open("w") as summary:
            summary.write(
                "input_records\treference_star_records\toutput_records\t"
                "removed_query_overlaps\tsamples\n"
            )
            star_records = sum(len(records) for records in by_sample.values())
            summary.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    input_records,
                    star_records,
                    len(selected),
                    star_records - len(selected),
                    len(by_sample),
                )
            )


if __name__ == "__main__":
    main()
