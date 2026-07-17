#!/usr/bin/env python3
"""Plot reference-position histograms for multiple SV VCFs."""

import argparse
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_events(path):
    opener = gzip.open if path.suffix == ".gz" else open
    events = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            kind = "INS" if len(fields[4]) > len(fields[3]) else "DEL"
            events.append((int(fields[1]), kind))
    return events


def bin_events(events, edges):
    return {
        key: np.histogram([pos for pos, kind in events if key == "total" or kind == key], bins=edges)[0]
        for key in ("total", "INS", "DEL")
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", action="append", type=Path, required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--reference-length", type=int, required=True)
    parser.add_argument("--bin-bp", type=int, default=50_000)
    parser.add_argument("--output-png", type=Path, required=True)
    parser.add_argument("--output-normalized-png", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()
    if len(args.vcf) != len(args.label):
        parser.error("--vcf and --label counts must match")

    edges = np.arange(0, args.reference_length + args.bin_bp, args.bin_bp)
    events = [read_events(path) for path in args.vcf]
    counts = [bin_events(items, edges) for items in events]
    centers = (edges[:-1] + edges[1:]) / 2_000_000
    widths = np.diff(edges) / 1_000_000 * 0.9
    line_colors = ["#2F6B9A", "#D97732", "#3A8D68", "#8B5A9F"]

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w") as handle:
        columns = ["bin_start", "bin_end"]
        for label in args.label:
            slug = label.lower().replace(" ", "_")
            columns.extend((f"{slug}_total", f"{slug}_ins", f"{slug}_del"))
        handle.write("\t".join(columns) + "\n")
        for index in range(len(edges) - 1):
            row = [edges[index], min(edges[index + 1], args.reference_length)]
            for item in counts:
                row.extend((item["total"][index], item["INS"][index], item["DEL"][index]))
            handle.write("\t".join(map(str, row)) + "\n")

    fig, axes = plt.subplots(len(events) + 1, 1, figsize=(15, 3 + 2.4 * len(events)), sharex=True)
    for index, (label, item, event_list) in enumerate(zip(args.label, counts, events)):
        axes[0].step(edges[:-1] / 1_000_000, item["total"], where="post", linewidth=1.7,
                     color=line_colors[index], label=f"{label} ({len(event_list)})")
    axes[0].set_title(f"CHM13 chr8 alpha-satellite SV start positions ({args.bin_bp // 1000} kb bins)")
    axes[0].set_ylabel("SV records / bin")
    axes[0].legend(frameon=False, ncol=len(events), fontsize=9)
    for axis, label, item in zip(axes[1:], args.label, counts):
        axis.bar(centers, item["INS"], width=widths, color="#2A9D8F", label="INS")
        axis.bar(centers, item["DEL"], width=widths, bottom=item["INS"], color="#E9A23B", label="DEL")
        axis.set_ylabel(f"{label}\nrecords")
        axis.legend(frameon=False, ncol=2, loc="upper right")
    axes[-1].set_xlabel("CHM13 reference position (Mb)")
    axes[-1].set_xlim(0, args.reference_length / 1_000_000)
    for axis in axes:
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(15, 4.8))
    for index, (label, item) in enumerate(zip(args.label, counts)):
        fractions = item["total"] / max(1, item["total"].sum())
        ax.step(edges[:-1] / 1_000_000, fractions, where="post", linewidth=1.8,
                color=line_colors[index], label=label)
    ax.set_title(f"Normalized CHM13 SV start-position distribution ({args.bin_bp // 1000} kb bins)")
    ax.set_xlabel("CHM13 reference position (Mb)")
    ax.set_ylabel("Fraction of method's SVs / bin")
    ax.set_xlim(0, args.reference_length / 1_000_000)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=len(events))
    fig.tight_layout()
    fig.savefig(args.output_normalized_png, dpi=180)


if __name__ == "__main__":
    main()
