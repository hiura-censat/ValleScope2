#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_events(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    events = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            position = int(fields[1])
            ref_length = len(fields[3])
            alt_length = len(fields[4])
            event_type = "INS" if alt_length > ref_length else "DEL"
            events.append((position, event_type))
    return events


def bin_events(events, edges):
    positions = np.array([event[0] for event in events], dtype=np.int64)
    insertions = np.array(
        [event[0] for event in events if event[1] == "INS"], dtype=np.int64
    )
    deletions = np.array(
        [event[0] for event in events if event[1] == "DEL"], dtype=np.int64
    )
    return {
        "all": np.histogram(positions, bins=edges)[0],
        "ins": np.histogram(insertions, bins=edges)[0],
        "del": np.histogram(deletions, bins=edges)[0],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Plot chm13-relative SV start-position histograms."
    )
    parser.add_argument("--raw-vcf", required=True, type=Path)
    parser.add_argument("--final-vcf", required=True, type=Path)
    parser.add_argument("--raw-label", default="Raw")
    parser.add_argument("--final-label", default="Final")
    parser.add_argument("--reference-length", required=True, type=int)
    parser.add_argument("--bin-bp", type=int, default=50_000)
    parser.add_argument("--output-png", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    args = parser.parse_args()

    edges = np.arange(0, args.reference_length + args.bin_bp, args.bin_bp)
    if edges[-1] < args.reference_length + 1:
        edges = np.append(edges, args.reference_length + 1)

    raw_events = read_events(args.raw_vcf)
    final_events = read_events(args.final_vcf)
    raw = bin_events(raw_events, edges)
    final = bin_events(final_events, edges)
    centers_mb = (edges[:-1] + edges[1:]) / 2_000_000
    widths_mb = np.diff(edges) / 1_000_000

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w") as handle:
        handle.write(
            "bin_start\tbin_end\traw_total\traw_ins\traw_del\t"
            "final_total\tfinal_ins\tfinal_del\n"
        )
        for index in range(len(edges) - 1):
            handle.write(
                f"{edges[index]}\t{min(edges[index + 1], args.reference_length)}\t"
                f"{raw['all'][index]}\t{raw['ins'][index]}\t{raw['del'][index]}\t"
                f"{final['all'][index]}\t{final['ins'][index]}\t"
                f"{final['del'][index]}\n"
            )

    colors = {"raw": "#3B6FB6", "final": "#D05A47", "ins": "#2A9D8F", "del": "#E9A23B"}
    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True)

    axes[0].step(edges[:-1] / 1_000_000, raw["all"], where="post", color=colors["raw"], linewidth=1.8, label=f"{args.raw_label} ({len(raw_events)})")
    axes[0].step(edges[:-1] / 1_000_000, final["all"], where="post", color=colors["final"], linewidth=1.8, label=f"{args.final_label} ({len(final_events)})")
    axes[0].set_ylabel("SV records / bin")
    axes[0].set_title(f"chm13 chr8 alpha-satellite SV start positions ({args.bin_bp // 1000} kb bins)")
    axes[0].legend(frameon=False, ncol=2)

    for axis, label, counts in ((axes[1], args.raw_label, raw), (axes[2], args.final_label, final)):
        axis.bar(centers_mb, counts["ins"], width=widths_mb * 0.9, color=colors["ins"], label="INS")
        axis.bar(centers_mb, counts["del"], width=widths_mb * 0.9, bottom=counts["ins"], color=colors["del"], label="DEL")
        axis.set_ylabel(f"{label} records")
        axis.legend(frameon=False, ncol=2, loc="upper right")

    axes[2].set_xlabel("chm13 reference position (Mb)")
    axes[2].set_xlim(0, args.reference_length / 1_000_000)
    for axis in axes:
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    fig.tight_layout()
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_png, dpi=180)


if __name__ == "__main__":
    main()
