#!/usr/bin/env python3
"""Plot graph SV positions after projecting concatenated reference coordinates."""

import argparse
import gzip
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


CIGAR = re.compile(r"(\d+)([=XMID])")


def paf_tags(fields):
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def load_coordinate_map(path):
    segments = []
    target_spans = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            tags = paf_tags(fields)
            cigar = tags.get("cg")
            if not cigar:
                continue
            strand = fields[4]
            qpos = int(fields[2]) if strand == "+" else int(fields[3])
            tpos = int(fields[7])
            target_spans.append((int(fields[7]), int(fields[8])))
            for length_text, operation in CIGAR.findall(cigar):
                length = int(length_text)
                if operation in "=XM":
                    if strand == "+":
                        segments.append((qpos, qpos + length, tpos, tpos + length, "+"))
                        qpos += length
                    else:
                        segments.append((qpos - length, qpos, tpos, tpos + length, "-"))
                        qpos -= length
                    tpos += length
                elif operation == "I":
                    if strand == "+":
                        segments.append((qpos, qpos + length, tpos, tpos, "+"))
                        qpos += length
                    else:
                        segments.append((qpos - length, qpos, tpos, tpos, "-"))
                        qpos -= length
                elif operation == "D":
                    tpos += length
    return sorted(segments), min(x[0] for x in target_spans), max(x[1] for x in target_spans)


def project(position, segments):
    position -= 1
    for qstart, qend, tstart, tend, strand in segments:
        if qstart <= position < qend:
            if tstart == tend:
                return tstart
            if strand == "+":
                return tstart + position - qstart
            return tstart + qend - 1 - position
    return None


def read_events(path, segments):
    opener = gzip.open if path.suffix == ".gz" else open
    events = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            position = project(int(fields[1]), segments)
            if position is None:
                continue
            ref_length = len(fields[3])
            deltas = [len(alt) - ref_length for alt in fields[4].split(",")
                      if not alt.startswith("<")]
            if not deltas:
                continue
            kind = "INS" if max(deltas, key=abs) > 0 else "DEL"
            events.append((position, kind))
    return events


def bin_events(events, edges):
    return {
        key: np.histogram(
            [pos for pos, kind in events if key == "total" or kind == key],
            bins=edges,
        )[0]
        for key in ("total", "INS", "DEL")
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", action="append", type=Path, required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--coordinate-map-paf", type=Path, required=True)
    parser.add_argument("--cdr-bed", type=Path, required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--bin-bp", type=int, default=50_000)
    parser.add_argument("--output-png", type=Path, required=True)
    parser.add_argument("--output-normalized-png", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    segments, plot_start, plot_end = load_coordinate_map(args.coordinate_map_paf)
    plot_start = plot_start // args.bin_bp * args.bin_bp
    plot_end = ((plot_end + args.bin_bp - 1) // args.bin_bp) * args.bin_bp
    edges = np.arange(plot_start, plot_end + args.bin_bp, args.bin_bp)
    events = [read_events(path, segments) for path in args.vcf]
    counts = [bin_events(items, edges) for items in events]
    highlights = []
    with args.cdr_bed.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields[0] == args.chromosome:
                highlights.append((int(fields[1]), int(fields[2]), fields[3]))

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w") as handle:
        columns = ["genomic_bin_start", "genomic_bin_end"]
        for label in args.label:
            slug = label.lower().replace(" ", "_")
            columns.extend((f"{slug}_total", f"{slug}_ins", f"{slug}_del"))
        handle.write("\t".join(columns) + "\n")
        for index in range(len(edges) - 1):
            row = [edges[index], edges[index + 1]]
            for item in counts:
                row.extend((item["total"][index], item["INS"][index], item["DEL"][index]))
            handle.write("\t".join(map(str, row)) + "\n")

    colors = ["#2F6B9A", "#D97732", "#3A8D68", "#8B5A9F"]
    centers = (edges[:-1] + edges[1:]) / 2_000_000
    widths = np.diff(edges) / 1_000_000 * 0.9
    fig, axes = plt.subplots(len(events) + 1, 1, figsize=(15, 3 + 2.4 * len(events)),
                             sharex=True)
    for axis in axes:
        for start, end, _ in highlights:
            axis.axvspan(start / 1_000_000, end / 1_000_000,
                         color="#8E6BBE", alpha=0.18, linewidth=0)
    for index, (label, item, records) in enumerate(zip(args.label, counts, events)):
        axes[0].step(edges[:-1] / 1_000_000, item["total"], where="post",
                     linewidth=1.7, color=colors[index],
                     label=f"{label} ({len(records)})")
    handles = axes[0].get_legend_handles_labels()[0]
    handles.append(Patch(facecolor="#8E6BBE", alpha=0.28, label="CDR"))
    axes[0].legend(handles=handles, frameon=False, ncol=len(events) + 1, fontsize=9)
    axes[0].set_title(
        f"CHM13 {args.chromosome} genomic SV start positions "
        f"({args.bin_bp // 1000} kb bins)"
    )
    axes[0].set_ylabel("SV records / bin")
    for axis, label, item in zip(axes[1:], args.label, counts):
        axis.bar(centers, item["INS"], width=widths, color="#2A9D8F", label="INS")
        axis.bar(centers, item["DEL"], width=widths, bottom=item["INS"],
                 color="#E9A23B", label="DEL")
        axis.set_ylabel(f"{label}\nrecords")
        axis.legend(frameon=False, ncol=2, loc="upper right")
    axes[-1].set_xlabel(f"CHM13 {args.chromosome} genomic position (Mb)")
    axes[-1].set_xlim(plot_start / 1_000_000, plot_end / 1_000_000)
    for axis in axes:
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(15, 4.8))
    for start, end, _ in highlights:
        axis.axvspan(start / 1_000_000, end / 1_000_000,
                     color="#8E6BBE", alpha=0.18, linewidth=0)
    for index, (label, item) in enumerate(zip(args.label, counts)):
        fractions = item["total"] / max(1, item["total"].sum())
        axis.step(edges[:-1] / 1_000_000, fractions, where="post",
                  linewidth=1.8, color=colors[index], label=label)
    handles = axis.get_legend_handles_labels()[0]
    handles.append(Patch(facecolor="#8E6BBE", alpha=0.28, label="CDR"))
    axis.legend(handles=handles, frameon=False, ncol=len(events) + 1)
    axis.set_title(
        f"Normalized CHM13 {args.chromosome} genomic SV distribution "
        f"({args.bin_bp // 1000} kb bins)"
    )
    axis.set_xlabel(f"CHM13 {args.chromosome} genomic position (Mb)")
    axis.set_ylabel("Fraction of method's SVs / bin")
    axis.set_xlim(plot_start / 1_000_000, plot_end / 1_000_000)
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output_normalized_png, dpi=180)


if __name__ == "__main__":
    main()
