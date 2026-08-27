#!/usr/bin/env python3
"""Compare three Pantree SV callsets and summarize frequency/fragmentation metrics."""

import argparse
import gzip
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle


@dataclass(frozen=True)
class Event:
    chrom: str
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
    def exact_key(self) -> Tuple[str, int, str, str]:
        return self.chrom, self.pos, self.ref.upper(), self.alt.upper()


def allele_length(allele: str) -> int:
    return 0 if allele == "." else len(allele)


def parse_info(value: str) -> Dict[str, str]:
    return {
        key: item_value
        for item in value.split(";")
        if "=" in item
        for key, item_value in [item.split("=", 1)]
    }


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
            kind = info.get("VT", "INS" if alt_length > ref_length else "DEL")
            events.append(
                Event(
                    fields[0],
                    int(fields[1]),
                    kind,
                    abs(alt_length - ref_length),
                    fields[3],
                    fields[4],
                    int(info.get("AC", "0").split(",")[0]),
                    int(info.get("AN", "0")),
                )
            )
    return events


def compatible(first: Event, second: Event, distance: int, size_ratio: float) -> bool:
    return (
        first.chrom == second.chrom
        and first.kind == second.kind
        and abs(first.pos - second.pos) <= distance
        and abs(first.end - second.end) <= distance
        and min(first.size, second.size) / max(first.size, second.size) >= size_ratio
    )


def maximum_match(
    left: List[Event], right: List[Event], distance: int, size_ratio: float
) -> List[int]:
    edges = [
        [
            right_index
            for right_index, second in enumerate(right)
            if compatible(first, second, distance, size_ratio)
        ]
        for first in left
    ]
    right_match = [-1] * len(right)
    left_match = [-1] * len(left)

    def augment(left_index: int, seen: List[bool]) -> bool:
        candidates = sorted(
            edges[left_index],
            key=lambda right_index: (
                abs(left[left_index].pos - right[right_index].pos)
                + abs(left[left_index].end - right[right_index].end),
                abs(left[left_index].size - right[right_index].size),
            ),
        )
        for right_index in candidates:
            if seen[right_index]:
                continue
            seen[right_index] = True
            if right_match[right_index] < 0 or augment(right_match[right_index], seen):
                left_match[left_index] = right_index
                right_match[right_index] = left_index
                return True
        return False

    for left_index in sorted(range(len(left)), key=lambda index: len(edges[index])):
        augment(left_index, [False] * len(right))
    return left_match


def three_way_matches(
    first: List[Event],
    second: List[Event],
    third: List[Event],
    distance: int,
    size_ratio: float,
) -> List[Tuple[int, int, int]]:
    candidates = []
    for i, event_a in enumerate(first):
        for j, event_b in enumerate(second):
            if not compatible(event_a, event_b, distance, size_ratio):
                continue
            for k, event_c in enumerate(third):
                if not compatible(event_a, event_c, distance, size_ratio):
                    continue
                if not compatible(event_b, event_c, distance, size_ratio):
                    continue
                cost = (
                    abs(event_a.pos - event_b.pos)
                    + abs(event_a.pos - event_c.pos)
                    + abs(event_b.pos - event_c.pos)
                    + abs(event_a.size - event_b.size)
                    + abs(event_a.size - event_c.size)
                )
                candidates.append((cost, i, j, k))
    used_a, used_b, used_c = set(), set(), set()
    matches = []
    for _, i, j, k in sorted(candidates):
        if i in used_a or j in used_b or k in used_c:
            continue
        used_a.add(i)
        used_b.add(j)
        used_c.add(k)
        matches.append((i, j, k))
    return matches


def percentile(values: List[int], fraction: float) -> float:
    return float(np.percentile(np.asarray(values), fraction * 100))


def adjacent_fraction(events: List[Event], distance: int) -> float:
    ordered = sorted(events, key=lambda event: (event.chrom, event.pos))
    close = set()
    for index in range(len(ordered) - 1):
        if (
            ordered[index].chrom == ordered[index + 1].chrom
            and ordered[index + 1].pos - ordered[index].pos <= distance
        ):
            close.update((index, index + 1))
    return len(close) / len(events) if events else 0.0


def draw_exact_venn(counts, labels, output):
    fig, ax = plt.subplots(figsize=(9, 7), dpi=180)
    for center, color in [
        ((0.39, 0.57), "#2F6B9A"),
        ((0.61, 0.57), "#D97732"),
        ((0.50, 0.37), "#3A8D68"),
    ]:
        ax.add_patch(Circle(center, 0.29, color=color, alpha=0.35))
    positions = {
        "100": (0.25, 0.64),
        "010": (0.75, 0.64),
        "001": (0.50, 0.18),
        "110": (0.50, 0.69),
        "101": (0.36, 0.40),
        "011": (0.64, 0.40),
        "111": (0.50, 0.50),
    }
    for region, position in positions.items():
        ax.text(*position, str(counts[region]), ha="center", va="center", fontsize=14)
    ax.text(0.18, 0.91, labels[0], ha="center", weight="bold")
    ax.text(0.82, 0.91, labels[1], ha="center", weight="bold")
    ax.text(0.50, 0.03, labels[2], ha="center", weight="bold")
    ax.set_title("Exact Pantree SV allele overlap")
    ax.set(xlim=(0, 1), ylim=(0, 1), aspect="equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", action="append", required=True, type=Path)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if len(args.vcf) != 3 or len(args.label) != 3:
        parser.error("exactly three --vcf and --label arguments are required")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    event_lists = [read_events(path) for path in args.vcf]
    with (args.output_dir / "callset_metrics.tsv").open("w") as handle:
        handle.write(
            "method\trecords\tunique_alleles\tduplicate_records\tINS\tDEL\t"
            "median_bp\tp90_bp\tmax_bp\tmean_AC\t"
            "mean_AN\tmean_AF\tAC1_fraction\tAF_ge_0.5_fraction\tAN_lt_10_fraction\t"
            "adjacent_100bp_fraction\n"
        )
        for label, events in zip(args.label, event_lists):
            sizes = [event.size for event in events]
            af = [event.ac / event.an for event in events if event.an]
            unique_alleles = len({event.exact_key for event in events})
            handle.write(
                f"{label}\t{len(events)}\t{unique_alleles}\t"
                f"{len(events) - unique_alleles}\t"
                f"{sum(event.kind == 'INS' for event in events)}\t"
                f"{sum(event.kind == 'DEL' for event in events)}\t"
                f"{percentile(sizes, 0.5):.1f}\t{percentile(sizes, 0.9):.1f}\t"
                f"{max(sizes)}\t{np.mean([event.ac for event in events]):.4f}\t"
                f"{np.mean([event.an for event in events]):.4f}\t{np.mean(af):.4f}\t"
                f"{np.mean([event.ac == 1 for event in events]):.4f}\t"
                f"{np.mean([value >= 0.5 for value in af]):.4f}\t"
                f"{np.mean([event.an < 10 for event in events]):.4f}\t"
                f"{adjacent_fraction(events, 100):.4f}\n"
            )

    exact_sets = [{event.exact_key for event in events} for events in event_lists]
    counts = Counter()
    all_exact_keys = set.union(*exact_sets)
    with (args.output_dir / "exact_membership.tsv").open("w") as handle:
        handle.write("chrom\tpos\tref\talt\t" + "\t".join(args.label) + "\n")
        for key in sorted(all_exact_keys, key=lambda value: (value[0], value[1], value[2], value[3])):
            membership = "".join(str(int(key in keys)) for keys in exact_sets)
            counts[membership] += 1
            handle.write("\t".join(map(str, key)) + "\t" + "\t".join(membership) + "\n")
    with (args.output_dir / "exact_overlap.tsv").open("w") as handle:
        handle.write("region\tcount\n")
        for region in ("100", "010", "001", "110", "101", "011", "111"):
            handle.write(f"{region}\t{counts[region]}\n")
    draw_exact_venn(counts, args.label, args.output_dir / "exact_overlap_venn.png")

    with (args.output_dir / "fuzzy_overlap.tsv").open("w") as handle:
        handle.write(
            "distance_bp\tminimum_size_ratio\tcomparison\tmatched\t"
            "left_total\tright_total\tleft_fraction\tright_fraction\n"
        )
        for distance in (50, 500, 1000):
            for left_index, right_index in ((0, 1), (0, 2), (1, 2)):
                matching = maximum_match(
                    event_lists[left_index],
                    event_lists[right_index],
                    distance,
                    0.7,
                )
                matched = sum(index >= 0 for index in matching)
                handle.write(
                    f"{distance}\t0.7\t{args.label[left_index]}__"
                    f"{args.label[right_index]}\t{matched}\t"
                    f"{len(event_lists[left_index])}\t{len(event_lists[right_index])}\t"
                    f"{matched / len(event_lists[left_index]):.6f}\t"
                    f"{matched / len(event_lists[right_index]):.6f}\n"
                )
            triples = three_way_matches(
                event_lists[0], event_lists[1], event_lists[2], distance, 0.7
            )
            if distance == 500:
                with (args.output_dir / "fuzzy_all_three_500bp.tsv").open("w") as triples_out:
                    triples_out.write(
                        "type\tValleScope2_pos\tValleScope2_size\tValleScope2_AC\t"
                        "ValleScope2_AN\tUnialigner_pos\tUnialigner_size\tUnialigner_AC\t"
                        "Unialigner_AN\tCentrolign_pos\tCentrolign_size\tCentrolign_AC\t"
                        "Centrolign_AN\n"
                    )
                    for first_index, second_index, third_index in triples:
                        matched_events = (
                            event_lists[0][first_index],
                            event_lists[1][second_index],
                            event_lists[2][third_index],
                        )
                        triples_out.write(
                            f"{matched_events[0].kind}\t"
                            + "\t".join(
                                str(value)
                                for event in matched_events
                                for value in (event.pos, event.size, event.ac, event.an)
                            )
                            + "\n"
                        )
            handle.write(
                f"{distance}\t0.7\tall_three\t{len(triples)}\t"
                f"{len(event_lists[0])}\t{len(event_lists[1])}\t"
                f"{len(triples) / len(event_lists[0]):.6f}\t"
                f"{len(triples) / len(event_lists[1]):.6f}\n"
            )

    with (args.output_dir / "ac_an_distribution.tsv").open("w") as handle:
        handle.write("method\tmetric\tvalue\tcount\tfraction\n")
        for label, events in zip(args.label, event_lists):
            for metric, values in (
                ("AC", [event.ac for event in events]),
                ("AN", [event.an for event in events]),
            ):
                distribution = Counter(values)
                for value in sorted(distribution):
                    handle.write(
                        f"{label}\t{metric}\t{value}\t{distribution[value]}\t"
                        f"{distribution[value] / len(events):.6f}\n"
                    )

    colors = ["#2F6B9A", "#D97732", "#3A8D68"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), dpi=180)
    x = np.arange(1, 11)
    width = 0.24
    for index, (label, events) in enumerate(zip(args.label, event_lists)):
        ac = Counter(event.ac for event in events)
        an = Counter(event.an for event in events)
        offset = (index - 1) * width
        axes[0].bar(x + offset, [ac[value] / len(events) for value in x], width,
                    label=label, color=colors[index])
        axes[1].bar(x + offset, [an[value] / len(events) for value in x], width,
                    label=label, color=colors[index])
        axes[2].hist(
            [event.ac / event.an for event in events if event.an],
            bins=np.linspace(0, 1, 11),
            histtype="step",
            linewidth=2,
            label=label,
            color=colors[index],
        )
    for axis, title, xlabel in (
        (axes[0], "Alternative allele count", "AC"),
        (axes[1], "Called allele count", "AN"),
        (axes[2], "Alternative allele frequency", "AF = AC / AN"),
    ):
        axis.set_title(title)
        axis.set_xlabel(xlabel)
        axis.set_ylabel("Fraction" if axis is not axes[2] else "SV records")
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False, fontsize=8)
    axes[0].set_xticks(x)
    axes[1].set_xticks(x)
    fig.tight_layout()
    fig.savefig(args.output_dir / "ac_an_af_distribution.png", bbox_inches="tight")
    plt.close(fig)

    size_edges = [50, 100, 500, 1000, 5000, 10000, 50000, math.inf]
    size_labels = ["50-99", "100-499", "500-999", "1k-4.9k", "5k-9.9k", "10k-49k", ">=50k"]
    fig, ax = plt.subplots(figsize=(11, 5), dpi=180)
    x = np.arange(len(size_labels))
    for index, (label, events) in enumerate(zip(args.label, event_lists)):
        values = Counter(
            next(i for i in range(len(size_edges) - 1)
                 if size_edges[i] <= event.size < size_edges[i + 1])
            for event in events
        )
        ax.bar(
            x + (index - 1) * width,
            [values[i] / len(events) for i in x],
            width,
            label=label,
            color=colors[index],
        )
    ax.set_xticks(x, size_labels)
    ax.set_ylabel("Fraction of SV records")
    ax.set_xlabel("SV size (bp)")
    ax.set_title("Pantree SV size distribution")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.output_dir / "sv_size_distribution.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
