#!/usr/bin/env python3
"""Plot reference-star PAF coverage against Pantree and vg SV density."""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_sv_positions(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [
            int(line.split("\t", 2)[1])
            for line in handle
            if not line.startswith("#")
        ]


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--reference-length", type=int, required=True)
    parser.add_argument("--pantree-vcf", type=Path, required=True)
    parser.add_argument("--vg-vcf", type=Path, required=True)
    parser.add_argument("--bin-bp", type=int, default=50_000)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-png", type=Path, required=True)
    args = parser.parse_args()

    by_sample = defaultdict(list)
    with args.paf.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields[0] == args.reference:
                by_sample[fields[5]].append((int(fields[2]), int(fields[3])))
            elif fields[5] == args.reference:
                by_sample[fields[0]].append((int(fields[7]), int(fields[8])))

    edges = np.arange(
        0, args.reference_length + args.bin_bp, args.bin_bp, dtype=np.int64
    )
    covered_bp = np.zeros(len(edges) - 1, dtype=np.int64)
    covered_samples = np.zeros(len(edges) - 1, dtype=np.int64)
    for intervals in by_sample.values():
        sample_covered = np.zeros(len(edges) - 1, dtype=np.int64)
        for start, end in merge_intervals(intervals):
            first = start // args.bin_bp
            last = min(len(sample_covered) - 1, (end - 1) // args.bin_bp)
            for index in range(first, last + 1):
                overlap = max(
                    0, min(end, edges[index + 1]) - max(start, edges[index])
                )
                sample_covered[index] += overlap
        sample_covered = np.minimum(sample_covered, np.diff(edges))
        covered_bp += sample_covered
        covered_samples += sample_covered > 0

    pantree = np.histogram(read_sv_positions(args.pantree_vcf), bins=edges)[0]
    vg = np.histogram(read_sv_positions(args.vg_vcf), bins=edges)[0]
    mean_depth = covered_bp / np.diff(edges)

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "bin_start",
                "bin_end",
                "direct_ref_aligned_sample_count",
                "mean_direct_ref_alignment_depth",
                "pantree_sv",
                "vg_sv",
            ]
        )
        for index in range(len(mean_depth)):
            writer.writerow(
                [
                    edges[index],
                    min(edges[index + 1], args.reference_length),
                    covered_samples[index],
                    mean_depth[index],
                    pantree[index],
                    vg[index],
                ]
            )

    x = edges[:-1] / 1_000_000
    fig, axes = plt.subplots(2, 1, figsize=(15, 7), sharex=True, dpi=180)
    axes[0].step(
        x,
        mean_depth,
        where="post",
        color="#2F6B9A",
        label="Mean direct CHM13 alignment depth",
    )
    axes[0].step(
        x,
        covered_samples,
        where="post",
        color="#3A8D68",
        label="Samples with direct CHM13 alignment",
    )
    axes[0].set_ylabel("Samples / 50 kb bin")
    axes[0].legend(frameon=False, ncol=2)
    axes[1].step(x, pantree, where="post", color="#D97732", label="Pantree")
    axes[1].step(x, vg, where="post", color="#8B5A9F", label="vg deconstruct")
    axes[1].set_ylabel("SV records / 50 kb bin")
    axes[1].set_xlabel("CHM13 reference position (Mb)")
    axes[1].legend(frameon=False, ncol=2)
    for axis in axes:
        axis.grid(axis="y", alpha=0.25)
        axis.spines[["top", "right"]].set_visible(False)
    axes[1].set_xlim(0, args.reference_length / 1_000_000)
    fig.tight_layout()
    fig.savefig(args.output_png)
    plt.close(fig)


if __name__ == "__main__":
    main()
