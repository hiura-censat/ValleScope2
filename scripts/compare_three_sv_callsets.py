#!/usr/bin/env python3
"""Summarize exact and breakpoint-compatible overlap among three SV VCFs."""

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle


@dataclass(frozen=True)
class Event:
    chrom: str
    pos: int
    ref: str
    alt: str
    record_id: str

    @property
    def kind(self):
        return "INS" if len(self.alt) > len(self.ref) else "DEL"

    @property
    def size(self):
        return abs(len(self.alt) - len(self.ref))

    @property
    def end(self):
        return self.pos if self.kind == "INS" else self.pos + self.size

    @property
    def key(self):
        return self.chrom, self.pos, self.ref.upper(), self.alt.upper()


def read_vcf(path):
    opener = gzip.open if path.suffix == ".gz" else open
    events = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            events.append(Event(fields[0], int(fields[1]), fields[3], fields[4], fields[2]))
    return events


def maximum_match(left, right, distance, size_ratio):
    edges = [[] for _ in left]
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            if first.chrom != second.chrom or first.kind != second.kind:
                continue
            if abs(first.pos - second.pos) > distance or abs(first.end - second.end) > distance:
                continue
            if min(first.size, second.size) / max(first.size, second.size) < size_ratio:
                continue
            edges[i].append(j)
    right_match = [-1] * len(right)
    left_match = [-1] * len(left)

    def augment(i, seen):
        for j in edges[i]:
            if seen[j]:
                continue
            seen[j] = True
            if right_match[j] < 0 or augment(right_match[j], seen):
                left_match[i] = j
                right_match[j] = i
                return True
        return False

    for i in sorted(range(len(left)), key=lambda index: len(edges[index])):
        augment(i, [False] * len(right))
    return left_match


def draw_venn(counts, labels, output):
    fig, ax = plt.subplots(figsize=(9, 7), dpi=180)
    circles = [
        Circle((0.39, 0.57), 0.29, color="#2F6B9A", alpha=0.35),
        Circle((0.61, 0.57), 0.29, color="#D97732", alpha=0.35),
        Circle((0.50, 0.37), 0.29, color="#3A8D68", alpha=0.35),
    ]
    for circle in circles:
        ax.add_patch(circle)
    positions = {
        "100": (0.25, 0.64), "010": (0.75, 0.64), "001": (0.50, 0.18),
        "110": (0.50, 0.69), "101": (0.36, 0.40), "011": (0.64, 0.40),
        "111": (0.50, 0.50),
    }
    for region, position in positions.items():
        ax.text(*position, str(counts[region]), ha="center", va="center", fontsize=14, weight="bold")
    ax.text(0.18, 0.91, labels[0], ha="center", fontsize=12, weight="bold", color="#17466B")
    ax.text(0.82, 0.91, labels[1], ha="center", fontsize=12, weight="bold", color="#914615")
    ax.text(0.50, 0.03, labels[2], ha="center", fontsize=12, weight="bold", color="#1F6045")
    ax.set_title("Exact normalized SV overlap (net length difference >= 50 bp, AC > 0)", fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", action="append", type=Path, required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--distance", type=int, default=500)
    parser.add_argument("--minimum-size-ratio", type=float, default=0.7)
    args = parser.parse_args()
    if len(args.vcf) != 3 or len(args.label) != 3:
        parser.error("exactly three --vcf and --label arguments are required")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    event_lists = [read_vcf(path) for path in args.vcf]
    key_sets = [{event.key for event in events} for events in event_lists]
    all_keys = set.union(*key_sets)
    counts = {f"{a}{b}{c}": 0 for a in (0, 1) for b in (0, 1) for c in (0, 1) if a or b or c}
    membership_path = args.output_dir / "exact_sv_membership.tsv"
    with membership_path.open("w") as handle:
        handle.write("chrom\tpos\tref\talt\t" + "\t".join(args.label) + "\n")
        for key in sorted(all_keys, key=lambda item: (item[0], item[1], item[2], item[3])):
            membership = tuple(int(key in keys) for keys in key_sets)
            region = "".join(map(str, membership))
            counts[region] += 1
            handle.write("\t".join(map(str, key + membership)) + "\n")

    with (args.output_dir / "exact_venn_summary.tsv").open("w") as handle:
        handle.write("region\tcount\n")
        for region in ("100", "010", "001", "110", "101", "011", "111"):
            handle.write(f"{region}\t{counts[region]}\n")
        for label, keys in zip(args.label, key_sets):
            handle.write(f"total:{label}\t{len(keys)}\n")
    draw_venn(counts, args.label, args.output_dir / "exact_sv_venn.png")

    match_01 = maximum_match(event_lists[0], event_lists[1], args.distance, args.minimum_size_ratio)
    match_02 = maximum_match(event_lists[0], event_lists[2], args.distance, args.minimum_size_ratio)
    match_12 = maximum_match(event_lists[1], event_lists[2], args.distance, args.minimum_size_ratio)
    shared_all = [i for i in range(len(event_lists[0])) if match_01[i] >= 0 and match_02[i] >= 0]
    with (args.output_dir / "fuzzy_shared_summary.tsv").open("w") as handle:
        handle.write("comparison\tmatched\tdistance_bp\tminimum_size_ratio\n")
        handle.write(f"{args.label[0]}__{args.label[1]}\t{sum(i >= 0 for i in match_01)}\t{args.distance}\t{args.minimum_size_ratio}\n")
        handle.write(f"{args.label[0]}__{args.label[2]}\t{sum(i >= 0 for i in match_02)}\t{args.distance}\t{args.minimum_size_ratio}\n")
        handle.write(f"{args.label[1]}__{args.label[2]}\t{sum(i >= 0 for i in match_12)}\t{args.distance}\t{args.minimum_size_ratio}\n")
        handle.write(f"all_three_{args.label[0]}_anchored\t{len(shared_all)}\t{args.distance}\t{args.minimum_size_ratio}\n")
    with (args.output_dir / "fuzzy_all_three_matches.tsv").open("w") as handle:
        handle.write("anchor_pos\tanchor_type\tanchor_size\tsecond_pos\tthird_pos\n")
        for i in shared_all:
            handle.write(
                f"{event_lists[0][i].pos}\t{event_lists[0][i].kind}\t{event_lists[0][i].size}\t"
                f"{event_lists[1][match_01[i]].pos}\t{event_lists[2][match_02[i]].pos}\n"
            )


if __name__ == "__main__":
    main()
