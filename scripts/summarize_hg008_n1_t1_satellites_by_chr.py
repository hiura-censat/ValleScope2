#!/usr/bin/env python3
"""Summarize per-chromosome HG008 N1/T1 satellite alignments and SV sharing."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = ("minimap2", "winnowmap2", "wfmash", "unialigner", "vallescope2")


def chrom_key(chrom):
    value = chrom.removeprefix("chr")
    return 23 if value == "X" else int(value)


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prep-root", required=True, type=Path)
    parser.add_argument("--alignment-root", required=True, type=Path)
    args = parser.parse_args()

    analysis = args.alignment_root / "analysis"
    plots = analysis / "plots_by_chr"
    plots.mkdir(parents=True, exist_ok=True)
    metrics = read_tsv(analysis / "alignment_metrics.tsv")
    events = read_tsv(analysis / "sv_candidates.tsv")
    clusters = read_tsv(analysis / "sv_clusters.tsv")
    intervals = read_tsv(args.prep_root / "extracted_intervals.tsv")
    selected = read_tsv(args.prep_root / "selected_contigs.tsv")

    interval = {(r["sample"], r["chrom"]): r for r in intervals}
    selection = {(r["sample"], r["chrom"]): r for r in selected}
    chroms = sorted(
        {r["pair"].split("_", 1)[0] for r in metrics}, key=chrom_key
    )

    summary = []
    for row in metrics:
        chrom = row["pair"].split("_", 1)[0]
        n1 = interval.get(("N1", chrom), {})
        t1 = interval.get(("T1", chrom), {})
        t1_selection = selection.get(("T1", chrom), {})
        summary.append({
            "chrom": chrom,
            "aligner": row["aligner"],
            "status": row["status"],
            "n1_alpha_bp": n1.get("length_bp", 0),
            "t1_alpha_bp": t1.get("length_bp", 0),
            "t1_chr_primary_aligned_bp": t1_selection.get("primary_aligned_bp", 0),
            "t1_cdr_mapping_distance_bp": t1.get("cdr_mapping_distance_bp", 0),
            "t1_cdr_mapping_mapq": t1.get("cdr_mapping_mapq", 0),
            "query_coverage": row["query_coverage"],
            "target_coverage": row["target_coverage"],
            "weighted_identity": row["weighted_identity"],
            "query_gap_count": row["query_gap_count"],
            "query_largest_gap_bp": row["query_largest_gap_bp"],
            "target_gap_count": row["target_gap_count"],
            "target_largest_gap_bp": row["target_largest_gap_bp"],
            "sv_count": row["sv_count"],
        })
    write_tsv(analysis / "chromosome_method_summary.tsv", summary)

    sharing = []
    for chrom in chroms:
        rows = [r for r in clusters if r["pair"].startswith(chrom + "_")]
        counts = Counter(int(r["aligner_count"]) for r in rows)
        sharing.append({
            "chrom": chrom,
            "clusters": len(rows),
            "single_method": counts[1],
            "shared_2_methods": counts[2],
            "shared_3_methods": counts[3],
            "shared_4_methods": counts[4],
            "shared_5_methods": counts[5],
            "shared_ge_2": sum(v for k, v in counts.items() if k >= 2),
            "shared_fraction": (
                sum(v for k, v in counts.items() if k >= 2) / len(rows)
                if rows else 0
            ),
        })
    write_tsv(analysis / "chromosome_sv_sharing.tsv", sharing)

    metric_lookup = {(r["pair"].split("_", 1)[0], r["aligner"]): r for r in metrics}
    coverage = np.array([
        [float(metric_lookup[(c, m)]["query_coverage"]) for m in METHODS]
        for c in chroms
    ])
    sv_counts = np.array([
        [int(metric_lookup[(c, m)]["sv_count"]) for m in METHODS]
        for c in chroms
    ])
    fig, axes = plt.subplots(1, 2, figsize=(13, 9))
    for axis, data, title, label in (
        (axes[0], coverage, "T1 alpha-satellite coverage", "fraction"),
        (axes[1], np.log10(sv_counts + 1), "SV count", "log10(count + 1)"),
    ):
        image = axis.imshow(data, aspect="auto", cmap="viridis")
        axis.set_xticks(range(len(METHODS)), METHODS, rotation=35, ha="right")
        axis.set_yticks(range(len(chroms)), chroms)
        axis.set_title(title)
        fig.colorbar(image, ax=axis, label=label, fraction=0.046)
    fig.tight_layout()
    fig.savefig(plots / "coverage_and_sv_count_heatmaps.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(7, 3, figsize=(15, 20))
    colors = dict(zip(METHODS, ("#3366A8", "#D97732", "#4C9A68", "#8B6BB1", "#C34A4A")))
    for axis, chrom in zip(axes.flat, chroms):
        n1_length = int(interval[("N1", chrom)]["length_bp"])
        bins = np.linspace(0, n1_length, 31)
        for method in METHODS:
            positions = [
                int(r["target_pos"]) for r in events
                if r["pair"].startswith(chrom + "_") and r["aligner"] == method
            ]
            if positions:
                axis.hist(positions, bins=bins, histtype="step",
                          linewidth=1.2, color=colors[method], label=method)
        axis.set_title(chrom)
        axis.set_xlim(0, n1_length)
        axis.tick_params(labelsize=8)
    for axis in axes.flat[len(chroms):]:
        axis.axis("off")
    handles = [
        plt.Line2D([0], [0], color=colors[m], label=m) for m in METHODS
    ]
    fig.legend(handles=handles, loc="upper center", ncol=5)
    fig.supxlabel("N1 alpha-satellite coordinate (bp)")
    fig.supylabel("SV candidates per bin")
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(plots / "sv_position_histograms_by_chromosome.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(13, 5))
    x = np.arange(len(chroms))
    shared = [int(r["shared_ge_2"]) for r in sharing]
    single = [int(r["single_method"]) for r in sharing]
    axis.bar(x, shared, label="shared by >=2 methods", color="#397A55")
    axis.bar(x, single, bottom=shared, label="single method", color="#A7B0B5")
    axis.set_xticks(x, chroms, rotation=45, ha="right")
    axis.set_ylabel("SV clusters")
    axis.legend()
    fig.tight_layout()
    fig.savefig(plots / "sv_sharing_by_chromosome.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
