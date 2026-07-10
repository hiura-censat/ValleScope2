#!/usr/bin/env python3
"""Compare truth-like PAF intervals with the best overlapping query PAF row."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


def parse_paf(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for field in fields[12:]:
                parts = field.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            rows.append(
                {
                    "line": line_number,
                    "q": fields[0],
                    "qs": int(fields[2]),
                    "qe": int(fields[3]),
                    "strand": fields[4],
                    "t": fields[5],
                    "ts": int(fields[7]),
                    "te": int(fields[8]),
                    "tv": tags.get("tv", "."),
                    "sv": tags.get("sv", "."),
                    "cigar": tags.get("cg", ""),
                }
            )
    return rows


def split_at_long_indels(row: dict[str, object], minimum_indel: int) -> list[dict[str, object]]:
    if minimum_indel <= 0 or not row["cigar"]:
        return [row]
    operations = [(int(length), op) for length, op in re.findall(r"(\d+)([=XMID])", row["cigar"])]
    if not operations:
        return [row]

    ref_cursor = row["ts"]
    query_cursor = row["qs"] if row["strand"] == "+" else row["qe"]
    block_ref_start = ref_cursor
    block_query_start = query_cursor
    blocks: list[dict[str, object]] = []

    def finish_block() -> None:
        nonlocal block_ref_start, block_query_start
        if ref_cursor == block_ref_start or query_cursor == block_query_start:
            return
        block = dict(row)
        block["ts"], block["te"] = sorted((block_ref_start, ref_cursor))
        block["qs"], block["qe"] = sorted((block_query_start, query_cursor))
        block["cigar"] = ""
        blocks.append(block)

    for length, op in operations:
        if op in "ID" and length >= minimum_indel:
            finish_block()
            if op == "D":
                ref_cursor += length
            elif row["strand"] == "+":
                query_cursor += length
            else:
                query_cursor -= length
            block_ref_start = ref_cursor
            block_query_start = query_cursor
            continue
        if op in "=XM":
            ref_cursor += length
            query_cursor += length if row["strand"] == "+" else -length
        elif op == "D":
            ref_cursor += length
        elif row["strand"] == "+":
            query_cursor += length
        else:
            query_cursor -= length
    finish_block()
    return blocks or [row]


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def boundary_deltas(truth: dict[str, object], call: dict[str, object]) -> tuple[int, int]:
    if truth["strand"] == "+":
        start = max(abs(truth["ts"] - call["ts"]), abs(truth["qs"] - call["qs"]))
        end = max(abs(truth["te"] - call["te"]), abs(truth["qe"] - call["qe"]))
    else:
        start = max(abs(truth["ts"] - call["ts"]), abs(truth["qe"] - call["qe"]))
        end = max(abs(truth["te"] - call["te"]), abs(truth["qs"] - call["qs"]))
    return start, end


def percentile(values: list[int], fraction: float) -> int:
    ordered = sorted(values)
    return ordered[math.floor(fraction * (len(ordered) - 1))]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth-paf", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--details", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument(
        "--split-cigar-indel-bp", type=int, default=0,
        help="Split call PAF rows into homology blocks at CIGAR I/D of at least this size",
    )
    args = parser.parse_args()

    truth_rows = parse_paf(args.truth_paf)
    call_rows = [
        block
        for row in parse_paf(args.paf)
        for block in split_at_long_indels(row, args.split_cigar_indel_bp)
    ]
    by_track: dict[tuple[object, object, object], list[dict[str, object]]] = defaultdict(list)
    for row in call_rows:
        by_track[(row["q"], row["strand"], row["t"])].append(row)

    fields = [
        "line", "q", "strand", "t", "ts", "te", "qs", "qe", "tv", "sv",
        "best_line", "ref_cov", "query_cov", "max_boundary_delta", "start_delta",
        "end_delta",
    ]
    counts = defaultdict(int)
    boundary_values: list[int] = []
    args.details.parent.mkdir(parents=True, exist_ok=True)
    with args.details.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for truth in truth_rows:
            best = None
            best_values = None
            for call in by_track[(truth["q"], truth["strand"], truth["t"])]:
                ref_cov = overlap(truth["ts"], truth["te"], call["ts"], call["te"]) / (
                    truth["te"] - truth["ts"]
                )
                query_cov = overlap(truth["qs"], truth["qe"], call["qs"], call["qe"]) / (
                    truth["qe"] - truth["qs"]
                )
                if ref_cov == 0 or query_cov == 0:
                    continue
                start_delta, end_delta = boundary_deltas(truth, call)
                max_delta = max(start_delta, end_delta)
                rank = (min(ref_cov, query_cov), ref_cov + query_cov, -max_delta)
                if best_values is None or rank > best_values[0]:
                    best = call
                    best_values = (rank, ref_cov, query_cov, max_delta, start_delta, end_delta)

            kind = "flank" if truth["tv"] == "flank" else "sv"
            counts[f"{kind}:total"] += 1
            output = {field: truth.get(field, ".") for field in fields[:10]}
            if best is None or best_values is None:
                counts[f"{kind}:no_overlap"] += 1
                output.update(
                    best_line=".", ref_cov="0", query_cov="0", max_boundary_delta=".",
                    start_delta=".", end_delta="."
                )
            else:
                _, ref_cov, query_cov, max_delta, start_delta, end_delta = best_values
                if ref_cov >= 0.95 and query_cov >= 0.95:
                    counts[f"{kind}:cover95"] += 1
                for threshold in (10, 50, 100, 500, 1000):
                    if max_delta <= threshold:
                        counts[f"{kind}:within_{threshold}bp"] += 1
                boundary_values.append(max_delta)
                output.update(
                    best_line=best["line"], ref_cov=f"{ref_cov:.6f}",
                    query_cov=f"{query_cov:.6f}", max_boundary_delta=max_delta,
                    start_delta=start_delta, end_delta=end_delta,
                )
            writer.writerow(output)

    summary_order = []
    for kind in ("flank", "sv"):
        summary_order.extend(
            [
                f"{kind}:cover95", f"{kind}:no_overlap", f"{kind}:total",
                f"{kind}:within_1000bp", f"{kind}:within_100bp",
                f"{kind}:within_10bp", f"{kind}:within_500bp", f"{kind}:within_50bp",
            ]
        )
    with args.summary.open("w") as handle:
        handle.write("metric\tvalue\n")
        for metric in summary_order:
            handle.write(f"{metric}\t{counts[metric]}\n")
        if boundary_values:
            handle.write(f"all:median_max_boundary_delta\t{percentile(boundary_values, 0.5)}\n")
            handle.write(f"all:p90_max_boundary_delta\t{percentile(boundary_values, 0.9)}\n")


if __name__ == "__main__":
    main()
