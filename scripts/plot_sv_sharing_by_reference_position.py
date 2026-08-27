#!/usr/bin/env python3
"""Plot cross-method SV sharing counts and rates along a reference."""

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_vcf(path):
    opener = gzip.open if path.suffix == ".gz" else open
    events = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            ref_len = len(fields[3])
            alt_len = len(fields[4])
            size = abs(alt_len - ref_len)
            kind = "INS" if alt_len > ref_len else "DEL"
            pos = int(fields[1])
            end = pos if kind == "INS" else pos + size
            events.append((pos, (kind, pos, end, size)))
    return events


def read_matched_keys(path):
    keys = Counter()
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            keys[
                (
                    row["type"],
                    int(row["pantree_pos"]),
                    int(row["pantree_end"]),
                    int(row["pantree_size"]),
                )
            ] += 1
    return keys


def histogram(events, selected, edges):
    positions = [pos for (pos, _key), keep in zip(events, selected) if keep]
    return np.histogram(positions, bins=edges)[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--centrolign-matches", type=Path, required=True)
    parser.add_argument("--unialigner-matches", type=Path, required=True)
    parser.add_argument("--reference-length", type=int, required=True)
    parser.add_argument("--bin-bp", type=int, default=50_000)
    parser.add_argument("--output-png", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    events = read_vcf(args.vcf)
    centrolign = read_matched_keys(args.centrolign_matches)
    unialigner = read_matched_keys(args.unialigner_matches)
    centrolign_remaining = centrolign.copy()
    unialigner_remaining = unialigner.copy()
    selected = {"centrolign": [], "unialigner": [], "either": [], "both": []}
    for _pos, key in events:
        in_centrolign = centrolign_remaining[key] > 0
        in_unialigner = unialigner_remaining[key] > 0
        if in_centrolign:
            centrolign_remaining[key] -= 1
        if in_unialigner:
            unialigner_remaining[key] -= 1
        selected["centrolign"].append(in_centrolign)
        selected["unialigner"].append(in_unialigner)
        selected["either"].append(in_centrolign or in_unialigner)
        selected["both"].append(in_centrolign and in_unialigner)
    edges = np.arange(0, args.reference_length + args.bin_bp, args.bin_bp)
    total = histogram(events, [True] * len(events), edges)
    counts = {
        name: histogram(events, selected[name], edges)
        for name in ("centrolign", "unialigner", "either", "both")
    }
    rates = {
        name: np.divide(values, total, out=np.zeros_like(values, dtype=float), where=total > 0)
        for name, values in counts.items()
    }

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "bin_start",
                "bin_end",
                "paffy_total",
                "shared_centrolign",
                "shared_unialigner",
                "shared_either",
                "shared_both",
                "rate_centrolign",
                "rate_unialigner",
                "rate_either",
                "rate_both",
            ]
        )
        for index in range(len(total)):
            writer.writerow(
                [
                    edges[index],
                    min(edges[index + 1], args.reference_length),
                    total[index],
                    counts["centrolign"][index],
                    counts["unialigner"][index],
                    counts["either"][index],
                    counts["both"][index],
                    *(rates[name][index] for name in ("centrolign", "unialigner", "either", "both")),
                ]
            )

    x = edges[:-1] / 1_000_000
    width = args.bin_bp / 1_000_000 * 0.88
    colors = {
        "centrolign": "#2F6B9A",
        "unialigner": "#D97732",
        "either": "#3A8D68",
        "both": "#805A9F",
    }
    fig, axes = plt.subplots(2, 1, figsize=(16, 8.5), sharex=True)
    axes[0].bar(x, total, width=width, align="edge", color="#D5D7DA", label="All paffy SVs")
    axes[0].step(x, counts["either"], where="post", color=colors["either"], linewidth=1.8,
                 label=f"Shared with either ({counts['either'].sum()})")
    axes[0].step(x, counts["both"], where="post", color=colors["both"], linewidth=1.8,
                 label=f"Shared with both ({counts['both'].sum()})")
    axes[0].set_ylabel("SV records / 50 kb")
    axes[0].set_title("Paffy chain+tile SV sharing along CHM13 chr8 alpha-satellite")
    axes[0].legend(frameon=False, ncol=3)

    for name, label in (
        ("centrolign", "Centrolign"),
        ("unialigner", "Unialigner"),
        ("either", "Either"),
        ("both", "Both"),
    ):
        axes[1].step(x, rates[name], where="post", linewidth=1.7, color=colors[name], label=label)
    axes[1].set_xlabel("CHM13 reference position (Mb)")
    axes[1].set_ylabel("Fraction of paffy SVs shared")
    axes[1].set_ylim(0, 1.02)
    axes[1].legend(frameon=False, ncol=4)
    for axis in axes:
        axis.set_xlim(0, args.reference_length / 1_000_000)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=180)


if __name__ == "__main__":
    main()
