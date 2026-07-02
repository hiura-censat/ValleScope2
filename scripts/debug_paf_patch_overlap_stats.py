#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([=XIDM])")


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def cigar_stats(cigar: str) -> dict[str, int]:
    stats = {"=": 0, "X": 0, "I": 0, "D": 0, "M": 0}
    for length, op in CIGAR_RE.findall(cigar):
        stats[op] += int(length)
    return stats


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if a_end <= b_start or b_end <= a_start:
        return 0
    return min(a_end, b_end) - max(a_start, b_start)


def span(start: int, end: int) -> int:
    return max(1, end - start)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paf")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--overlap-threshold", type=float, default=0.01)
    args = parser.parse_args()

    records: list[dict] = []
    patch_counter: Counter[str] = Counter()
    method_counter: Counter[str] = Counter()
    total_cigar = Counter()
    total_pg = 0

    with open(args.paf) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            tags = parse_tags(fields)
            cigar = tags.get("cg", "")
            cstats = cigar_stats(cigar)
            total_cigar.update(cstats)
            bt = tags.get("bt", ".")
            am = tags.get("am", ".")
            pg = int(tags.get("pg", "0"))
            total_pg += pg
            patch_counter[bt] += 1
            method_counter[am] += 1
            records.append(
                {
                    "line": line_number,
                    "qname": fields[0],
                    "qlen": int(fields[1]),
                    "qstart": int(fields[2]),
                    "qend": int(fields[3]),
                    "strand": fields[4],
                    "tname": fields[5],
                    "tlen": int(fields[6]),
                    "tstart": int(fields[7]),
                    "tend": int(fields[8]),
                    "matches": int(fields[9]),
                    "block": int(fields[10]),
                    "bt": bt,
                    "pg": pg,
                    "am": am,
                    "ch": tags.get("ch", "."),
                    "ins": cstats["I"],
                    "del": cstats["D"],
                    "sub": cstats["X"],
                    "eq": cstats["="] + cstats["M"],
                }
            )

    by_track: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for rec in records:
        by_track[(rec["qname"], rec["tname"], rec["strand"])].append(rec)

    overlap_rows: list[dict] = []
    overlap_summary = Counter()
    max_ref_ratio = 0.0
    max_query_ratio = 0.0
    for track, recs in by_track.items():
        recs.sort(key=lambda r: (r["tstart"], r["qstart"], r["line"]))
        for i in range(len(recs)):
            a = recs[i]
            for j in range(i + 1, len(recs)):
                b = recs[j]
                if b["tstart"] >= a["tend"] and b["qstart"] >= a["qend"]:
                    break
                ref_ov = overlap(a["tstart"], a["tend"], b["tstart"], b["tend"])
                query_ov = overlap(a["qstart"], a["qend"], b["qstart"], b["qend"])
                if ref_ov == 0 and query_ov == 0:
                    continue
                ref_ratio = ref_ov / min(span(a["tstart"], a["tend"]), span(b["tstart"], b["tend"]))
                query_ratio = query_ov / min(span(a["qstart"], a["qend"]), span(b["qstart"], b["qend"]))
                max_ref_ratio = max(max_ref_ratio, ref_ratio)
                max_query_ratio = max(max_query_ratio, query_ratio)
                if ref_ratio >= args.overlap_threshold and query_ratio >= args.overlap_threshold:
                    overlap_summary["both_ge_threshold"] += 1
                elif ref_ratio >= args.overlap_threshold:
                    overlap_summary["ref_only_ge_threshold"] += 1
                elif query_ratio >= args.overlap_threshold:
                    overlap_summary["query_only_ge_threshold"] += 1
                else:
                    overlap_summary["below_threshold"] += 1
                overlap_rows.append(
                    {
                        "qname": track[0],
                        "tname": track[1],
                        "strand": track[2],
                        "line_a": a["line"],
                        "line_b": b["line"],
                        "chain_a": a["ch"],
                        "chain_b": b["ch"],
                        "bt_a": a["bt"],
                        "bt_b": b["bt"],
                        "ref_overlap": ref_ov,
                        "query_overlap": query_ov,
                        "ref_overlap_ratio": f"{ref_ratio:.6f}",
                        "query_overlap_ratio": f"{query_ratio:.6f}",
                    }
                )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{out_prefix}.overlaps.tsv", "w", newline="") as handle:
        fieldnames = [
            "qname",
            "tname",
            "strand",
            "line_a",
            "line_b",
            "chain_a",
            "chain_b",
            "bt_a",
            "bt_b",
            "ref_overlap",
            "query_overlap",
            "ref_overlap_ratio",
            "query_overlap_ratio",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(overlap_rows)

    patched_records = patch_counter["patched"]
    total_records = len(records)
    patch_gap_records = sum(1 for rec in records if rec["pg"] > 0)
    with open(f"{out_prefix}.summary.tsv", "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"records\t{total_records}\n")
        out.write(f"patched_bt_records\t{patched_records}\n")
        out.write(f"patched_bt_fraction\t{patched_records / total_records if total_records else 0:.6f}\n")
        out.write(f"records_with_pg_gt0\t{patch_gap_records}\n")
        out.write(f"records_with_pg_gt0_fraction\t{patch_gap_records / total_records if total_records else 0:.6f}\n")
        out.write(f"total_pg_tag_sum\t{total_pg}\n")
        for key, value in sorted(patch_counter.items()):
            out.write(f"bt:{key}\t{value}\n")
        for key, value in sorted(method_counter.items()):
            out.write(f"am:{key}\t{value}\n")
        for key in ["=", "M", "X", "I", "D"]:
            out.write(f"cigar_{key}\t{total_cigar[key]}\n")
        out.write(f"overlap_pairs_total\t{len(overlap_rows)}\n")
        for key, value in sorted(overlap_summary.items()):
            out.write(f"overlap_{key}\t{value}\n")
        out.write(f"max_ref_overlap_ratio\t{max_ref_ratio:.6f}\n")
        out.write(f"max_query_overlap_ratio\t{max_query_ratio:.6f}\n")

    print(f"{out_prefix}.summary.tsv")
    print(f"{out_prefix}.overlaps.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
