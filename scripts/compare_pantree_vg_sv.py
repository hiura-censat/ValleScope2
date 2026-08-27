#!/usr/bin/env python3
"""Compare Pantree and vg deconstruct SV VCFs with one-to-one event matching."""

import argparse
import csv
import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Event:
    pos: int
    kind: str
    size: int
    ref: str
    alt: str
    ac: int
    an: int

    @property
    def end(self) -> int:
        return self.pos if self.kind == "INS" else self.pos + self.size

    @property
    def normalized_key(self) -> Tuple[int, str, str]:
        pos, ref, alt = normalize_alleles(self.pos, self.ref, self.alt)
        return pos, ref, alt


def parse_info(value: str) -> Dict[str, str]:
    return {
        key: item_value
        for item in value.split(";")
        if "=" in item
        for key, item_value in [item.split("=", 1)]
    }


def allele_length(allele: str) -> int:
    return 0 if allele == "." else len(allele)


def normalize_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    ref = ref.upper()
    alt = alt.upper()
    if ref == "." or alt == ".":
        return pos, ref, alt
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return pos, ref, alt


def read_events(path: Path) -> List[Event]:
    opener = gzip.open if path.suffix == ".gz" else open
    events = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            info = parse_info(fields[7])
            ref_length = allele_length(fields[3])
            alt_length = allele_length(fields[4])
            kind = info.get("VT", info.get("SVTYPE", ""))
            if kind not in ("INS", "DEL"):
                kind = "INS" if alt_length > ref_length else "DEL"
            events.append(
                Event(
                    pos=int(fields[1]),
                    kind=kind,
                    size=abs(alt_length - ref_length),
                    ref=fields[3],
                    alt=fields[4],
                    ac=int(info.get("AC", "0").split(",")[0]),
                    an=int(info.get("AN", "0")),
                )
            )
    return events


def compatible(left: Event, right: Event, distance: int, size_ratio: float) -> bool:
    if left.kind != right.kind or max(left.size, right.size) == 0:
        return False
    return (
        abs(left.pos - right.pos) <= distance
        and abs(left.end - right.end) <= distance
        and min(left.size, right.size) / max(left.size, right.size) >= size_ratio
    )


def maximum_match(
    left: List[Event], right: List[Event], distance: int, size_ratio: float
) -> List[int]:
    edges = [
        [
            index
            for index, candidate in enumerate(right)
            if compatible(event, candidate, distance, size_ratio)
        ]
        for event in left
    ]
    right_match = [-1] * len(right)
    left_match = [-1] * len(left)

    def augment(left_index: int, seen: List[bool]) -> bool:
        ordered = sorted(
            edges[left_index],
            key=lambda right_index: (
                abs(left[left_index].pos - right[right_index].pos)
                + abs(left[left_index].end - right[right_index].end),
                abs(left[left_index].size - right[right_index].size),
            ),
        )
        for right_index in ordered:
            if seen[right_index]:
                continue
            seen[right_index] = True
            previous = right_match[right_index]
            if previous < 0 or augment(previous, seen):
                left_match[left_index] = right_index
                right_match[right_index] = left_index
                return True
        return False

    for left_index in sorted(range(len(left)), key=lambda index: len(edges[index])):
        augment(left_index, [False] * len(right))
    return left_match


def exact_match(left: List[Event], right: List[Event]) -> List[int]:
    buckets = {}
    for right_index, event in enumerate(right):
        buckets.setdefault(event.normalized_key, []).append(right_index)
    matches = [-1] * len(left)
    for left_index, event in enumerate(left):
        candidates = buckets.get(event.normalized_key, [])
        if candidates:
            matches[left_index] = candidates.pop()
    return matches


def deduplicate(events: List[Event]) -> List[Event]:
    unique = {}
    for event in events:
        key = (event.kind,) + event.normalized_key
        previous = unique.get(key)
        if previous is None or event.ac > previous.ac:
            unique[key] = event
    return list(unique.values())


def summarize_matches(
    method: str,
    criterion: str,
    pantree: List[Event],
    vg: List[Event],
    matches: List[int],
) -> Dict[str, object]:
    pairs = [(pantree[index], vg[right]) for index, right in enumerate(matches) if right >= 0]
    count = len(pairs)
    return {
        "method": method,
        "criterion": criterion,
        "pantree_n": len(pantree),
        "vg_n": len(vg),
        "matched": count,
        "pantree_match_rate": count / len(pantree) if pantree else 0.0,
        "vg_match_rate": count / len(vg) if vg else 0.0,
        "overlap_f1": 2 * count / (len(pantree) + len(vg))
        if pantree or vg
        else 0.0,
        "jaccard_like": count / (len(pantree) + len(vg) - count)
        if pantree or vg
        else 0.0,
        "matched_ins": sum(first.kind == "INS" for first, _ in pairs),
        "matched_del": sum(first.kind == "DEL" for first, _ in pairs),
        "ac_exact_rate": sum(first.ac == second.ac for first, second in pairs) / count
        if count
        else 0.0,
        "mean_ac_abs_diff": sum(abs(first.ac - second.ac) for first, second in pairs)
        / count
        if count
        else 0.0,
    }


def write_details(
    path: Path, pantree: List[Event], vg: List[Event], matches: List[int]
) -> None:
    with path.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "type",
                "pantree_pos",
                "pantree_end",
                "pantree_size",
                "pantree_AC",
                "pantree_AN",
                "vg_pos",
                "vg_end",
                "vg_size",
                "vg_AC",
                "vg_AN",
                "start_diff",
                "end_diff",
                "size_ratio",
            ]
        )
        for index, right_index in enumerate(matches):
            if right_index < 0:
                continue
            first, second = pantree[index], vg[right_index]
            writer.writerow(
                [
                    first.kind,
                    first.pos,
                    first.end,
                    first.size,
                    first.ac,
                    first.an,
                    second.pos,
                    second.end,
                    second.size,
                    second.ac,
                    second.an,
                    first.pos - second.pos,
                    first.end - second.end,
                    min(first.size, second.size) / max(first.size, second.size),
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", action="append", required=True)
    parser.add_argument("--pantree", action="append", required=True, type=Path)
    parser.add_argument("--vg", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--size-ratio", type=float, default=0.7)
    args = parser.parse_args()
    if not (len(args.method) == len(args.pantree) == len(args.vg)):
        parser.error("--method, --pantree, and --vg must have equal counts")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    details = {}
    counts = {}
    events_by_method = {}
    for method, pantree_path, vg_path in zip(args.method, args.pantree, args.vg):
        pantree = read_events(pantree_path)
        vg = read_events(vg_path)
        events_by_method[method] = (pantree, vg)
        counts[method] = (
            Counter(event.kind for event in pantree),
            Counter(event.kind for event in vg),
        )
        rows.append(
            summarize_matches(method, "normalized_exact", pantree, vg, exact_match(pantree, vg))
        )
        for distance in (50, 500, 1000):
            matches = maximum_match(pantree, vg, distance, args.size_ratio)
            rows.append(
                summarize_matches(
                    method, "bp{}_size0.70".format(distance), pantree, vg, matches
                )
            )
            if distance == 500:
                details[method] = (pantree, vg, matches)
                write_details(
                    args.output_dir / "{}.bp500.matches.tsv".format(method.lower()),
                    pantree,
                    vg,
                    matches,
                )

    with (args.output_dir / "bp500_unique_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for method in args.method:
            pantree, vg = events_by_method[method]
            pantree = deduplicate(pantree)
            vg = deduplicate(vg)
            writer.writerow(
                summarize_matches(
                    method,
                    "bp500_size0.70_unique",
                    pantree,
                    vg,
                    maximum_match(pantree, vg, 500, args.size_ratio),
                )
            )

    fields = list(rows[0])
    with (args.output_dir / "comparison_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    with (args.output_dir / "bp500_type_summary.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "method",
                "type",
                "pantree_n",
                "vg_n",
                "matched",
                "pantree_match_rate",
                "vg_match_rate",
            ]
        )
        for method in args.method:
            pantree, vg, matches = details[method]
            for kind in ("INS", "DEL"):
                matched = sum(
                    right >= 0 and pantree[index].kind == kind
                    for index, right in enumerate(matches)
                )
                pantree_n = counts[method][0][kind]
                vg_n = counts[method][1][kind]
                writer.writerow(
                    [
                        method,
                        kind,
                        pantree_n,
                        vg_n,
                        matched,
                        matched / pantree_n if pantree_n else 0.0,
                        matched / vg_n if vg_n else 0.0,
                    ]
                )

    criteria = ["normalized_exact", "bp50_size0.70", "bp500_size0.70", "bp1000_size0.70"]
    x = list(range(len(args.method)))
    width = 0.19
    fig, ax = plt.subplots(figsize=(10, 5.8), dpi=180)
    colors = ["#555555", "#2F6B9A", "#3A8D68", "#D97732"]
    for offset, (criterion, color) in enumerate(zip(criteria, colors)):
        values = [
            next(
                row["pantree_match_rate"]
                for row in rows
                if row["method"] == method and row["criterion"] == criterion
            )
            for method in args.method
        ]
        positions = [value + (offset - 1.5) * width for value in x]
        ax.bar(positions, values, width, label=criterion, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(args.method)
    ax.set_ylabel("Fraction of Pantree SVs matched to vg")
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(args.output_dir / "pantree_vs_vg_match_rates.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
