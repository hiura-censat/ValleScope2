#!/usr/bin/env python3
"""Trace gaps between final PAF records back to ValleScope2 debug decisions."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


class DisjointSet:
    def __init__(self) -> None:
        self.parent: dict[int, int] = {}

    def find(self, value: int) -> int:
        self.parent.setdefault(value, value)
        if self.parent[value] != value:
            self.parent[value] = self.find(self.parent[value])
        return self.parent[value]

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--chain-extensions", required=True, type=Path)
    parser.add_argument("--patch-intervals", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top-per-sample", type=int, default=20)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_paf(path: Path) -> list[dict[str, object]]:
    records = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{line_number}: fewer than 12 PAF columns")
            tags = {}
            for field in fields[12:]:
                parts = field.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            records.append(
                {
                    "query": fields[0],
                    "query_length": int(fields[1]),
                    "query_start": int(fields[2]),
                    "query_end": int(fields[3]),
                    "strand": fields[4],
                    "reference": fields[5],
                    "reference_length": int(fields[6]),
                    "reference_start": int(fields[7]),
                    "reference_end": int(fields[8]),
                    "chain_id": int(tags.get("ch", "-1")),
                    "block_type": tags.get("bt", "."),
                }
            )
    return records


def integer(row: dict[str, str] | None, key: str, default: int = 0) -> int:
    if not row or row.get(key, "") in {"", "."}:
        return default
    return int(row[key])


def classify(
    patch: dict[str, str] | None, extension: dict[str, str] | None
) -> tuple[str, str]:
    if patch:
        status = patch["status"]
        if status == "patch_interval_low_identity":
            return status, "WFA completed, but the alignment did not meet patch acceptance"
        if status == "patch_zdrop_reject":
            return status, "long-I/D rescue was considered, but Z-drop rejected the noisy alignment"
        if patch.get("merge") == "1":
            return "accepted_patch_boundary_remains", "patch merged; the final gap is after trimming or chain lineage changes"
        return status, "patch was evaluated but not merged"
    if extension:
        if extension["stop_reason"] == "gap_too_large":
            return "not_patch_eligible_gap_too_large", "at least one original gap exceeded the 70 kb extension/patch limit"
        return "no_patch_record_after_extension", f"extension stopped with {extension['stop_reason']}; no matching patch decision was recorded"
    return "merged_lineage_or_nonadjacent_chain", "final neighboring records were not direct neighbors in the pre-merge chain graph"


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    extensions = read_tsv(args.chain_extensions)
    patches = read_tsv(args.patch_intervals)
    extension_by_edge = {
        (int(row["left_chain_id"]), int(row["right_chain_id"])): row
        for row in extensions
    }
    patch_by_edge = {
        (int(row["left_chain_id"]), int(row["right_chain_id"])): row
        for row in patches
    }
    merged = DisjointSet()
    for row in patches:
        if row.get("merge") == "1":
            merged.union(int(row["left_chain_id"]), int(row["right_chain_id"]))

    extensions_by_component: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
    patches_by_component: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
    for row in extensions:
        key = (merged.find(int(row["left_chain_id"])), merged.find(int(row["right_chain_id"])))
        if key[0] != key[1]:
            extensions_by_component[key].append(row)
    for row in patches:
        key = (merged.find(int(row["left_chain_id"])), merged.find(int(row["right_chain_id"])))
        if key[0] != key[1]:
            patches_by_component[key].append(row)

    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for record in parse_paf(args.paf):
        grouped[(str(record["query"]), str(record["reference"]), str(record["strand"]))].append(record)

    gaps: list[dict[str, object]] = []
    for (sample, reference, strand), records in grouped.items():
        records.sort(key=lambda row: (int(row["reference_start"]), int(row["reference_end"])))
        for left, right in zip(records, records[1:]):
            ref_gap = int(right["reference_start"]) - int(left["reference_end"])
            if strand == "+":
                query_gap = int(right["query_start"]) - int(left["query_end"])
            else:
                query_gap = int(left["query_start"]) - int(right["query_end"])
            if ref_gap <= 0 and query_gap <= 0:
                continue

            edge = (int(left["chain_id"]), int(right["chain_id"]))
            extension = extension_by_edge.get(edge)
            patch = patch_by_edge.get(edge)
            decision_match = "exact" if extension or patch else "."
            component_edge = (merged.find(edge[0]), merged.find(edge[1]))

            def gap_distance(row: dict[str, str]) -> int:
                return abs(integer(row, "ref_gap", ref_gap) - ref_gap) + abs(
                    integer(row, "query_gap", query_gap) - query_gap
                )

            if extension is None:
                candidates = extensions_by_component.get(component_edge, [])
                if candidates:
                    extension = min(candidates, key=gap_distance)
                    decision_match = "merged_component"
            if extension is None:
                # A competing chain at the same locus can replace the chain used by
                # extension. Recover that decision by requiring the same left merged
                # component and near-identical right-chain start coordinates.
                coordinate_candidates = []
                for row in extensions:
                    if merged.find(int(row["left_chain_id"])) != component_edge[0]:
                        continue
                    if not row["sample_b"].startswith(sample):
                        continue
                    ref_delta = abs(integer(row, "right_ref_start") - int(right["reference_start"]))
                    query_delta = abs(integer(row, "right_query_start") - int(right["query_start"]))
                    if ref_delta <= 5000 and query_delta <= 5000:
                        coordinate_candidates.append((ref_delta + query_delta, row))
                if coordinate_candidates:
                    extension = min(coordinate_candidates, key=lambda item: item[0])[1]
                    decision_match = "competing_chain_coordinate_proxy"
            if patch is None and extension is not None:
                extension_edge = (int(extension["left_chain_id"]), int(extension["right_chain_id"]))
                patch = patch_by_edge.get(extension_edge)
            if patch is None:
                candidates = patches_by_component.get(component_edge, [])
                if candidates:
                    patch = min(candidates, key=gap_distance)
            cause, explanation = classify(patch, extension)
            gaps.append(
                {
                    "sample": sample,
                    "reference": reference,
                    "strand": strand,
                    "left_chain_id": edge[0],
                    "right_chain_id": edge[1],
                    "decision_left_chain_id": extension["left_chain_id"] if extension else (patch["left_chain_id"] if patch else "."),
                    "decision_right_chain_id": extension["right_chain_id"] if extension else (patch["right_chain_id"] if patch else "."),
                    "decision_match": decision_match,
                    "extension_left_chain_id": extension["left_chain_id"] if extension else ".",
                    "extension_right_chain_id": extension["right_chain_id"] if extension else ".",
                    "patch_left_chain_id": patch["left_chain_id"] if patch else ".",
                    "patch_right_chain_id": patch["right_chain_id"] if patch else ".",
                    "left_block_type": left["block_type"],
                    "right_block_type": right["block_type"],
                    "reference_left_end": left["reference_end"],
                    "reference_right_start": right["reference_start"],
                    "query_left_boundary": left["query_end"] if strand == "+" else left["query_start"],
                    "query_right_boundary": right["query_start"] if strand == "+" else right["query_end"],
                    "final_ref_gap": ref_gap,
                    "final_query_gap": query_gap,
                    "max_final_gap": max(ref_gap, query_gap),
                    "extension_classification": extension["classification"] if extension else ".",
                    "extension_stop_reason": extension["stop_reason"] if extension else ".",
                    "original_ref_gap": integer(extension, "ref_gap", integer(patch, "ref_gap", -1)),
                    "original_query_gap": integer(extension, "query_gap", integer(patch, "query_gap", -1)),
                    "patch_status": patch["status"] if patch else ".",
                    "patch_identity": patch["identity"] if patch else ".",
                    "patch_merge": patch["merge"] if patch else ".",
                    "zdrop_fired": patch["zdrop_fired"] if patch else ".",
                    "cause": cause,
                    "explanation": explanation,
                }
            )

    gaps.sort(key=lambda row: (str(row["sample"]), int(row["reference_left_end"])))
    fields = list(gaps[0]) if gaps else []
    with (args.output_dir / "all_final_paf_gaps.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(gaps)

    summary_rows = []
    samples = sorted({str(row["sample"]) for row in gaps})
    causes = sorted({str(row["cause"]) for row in gaps})
    for sample in samples:
        sample_rows = [row for row in gaps if row["sample"] == sample]
        for cause in causes:
            rows = [row for row in sample_rows if row["cause"] == cause]
            if not rows:
                continue
            summary_rows.append(
                {
                    "sample": sample,
                    "cause": cause,
                    "gap_count": len(rows),
                    "gap_count_ge_1kb": sum(int(row["max_final_gap"]) >= 1000 for row in rows),
                    "max_ref_gap": max(int(row["final_ref_gap"]) for row in rows),
                    "max_query_gap": max(int(row["final_query_gap"]) for row in rows),
                }
            )
    summary_fields = ["sample", "cause", "gap_count", "gap_count_ge_1kb", "max_ref_gap", "max_query_gap"]
    with (args.output_dir / "gap_cause_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    top_rows = []
    for sample in samples:
        rows = [row for row in gaps if row["sample"] == sample]
        top_rows.extend(sorted(rows, key=lambda row: int(row["max_final_gap"]), reverse=True)[: args.top_per_sample])
    with (args.output_dir / "top_gaps.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(top_rows)

    print(f"gaps={len(gaps)}")
    print("causes=" + ", ".join(f"{key}:{value}" for key, value in sorted(Counter(row['cause'] for row in gaps).items())))


if __name__ == "__main__":
    main()
