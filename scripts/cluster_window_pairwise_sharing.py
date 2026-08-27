#!/usr/bin/env python3
"""Cluster haplotypes by bp-weighted Jaccard in CHM13-anchored windows."""

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from analyze_gfa_path_sharing import SAMPLE_ORDER, parse_gfa, sample_name


def load_window_jaccard(method, gfa, window_bp):
    lengths, paths = parse_gfa(gfa)
    sample_paths = {}
    for path_name, steps in paths.items():
        sample = sample_name(path_name)
        if sample in SAMPLE_ORDER:
            sample_paths[sample] = steps
    missing = [sample for sample in SAMPLE_ORDER if sample not in sample_paths]
    if missing:
        raise ValueError(f"{method}: missing paths: {','.join(missing)}")
    counts = {sample: Counter(sample_paths[sample]) for sample in SAMPLE_ORDER}
    chm_nodes = sample_paths["chm13v2.0"]
    path_length = sum(lengths[node] for node in chm_nodes)
    n_windows = (path_length + window_bp - 1) // window_bp
    intersections = np.zeros((n_windows, len(SAMPLE_ORDER), len(SAMPLE_ORDER)), dtype=float)
    unions = np.zeros_like(intersections)
    pos = 0
    for node in chm_nodes:
        end = pos + lengths[node]
        presence = np.array([counts[sample][node] > 0 for sample in SAMPLE_ORDER], dtype=bool)
        both = np.logical_and.outer(presence, presence)
        either = np.logical_or.outer(presence, presence)
        first = pos // window_bp
        last = (end - 1) // window_bp
        for index in range(first, last + 1):
            overlap = max(0, min(end, (index + 1) * window_bp) - max(pos, index * window_bp))
            intersections[index] += overlap * both
            unions[index] += overlap * either
        pos = end
    matrices = np.ones_like(intersections)
    np.divide(intersections, unions, out=matrices, where=unions > 0)
    return path_length, matrices


def canonical_partition(labels):
    groups = []
    for label in sorted(set(labels)):
        groups.append(tuple(index for index, value in enumerate(labels) if value == label))
    return tuple(sorted(groups))


def cluster_matrix(matrix, minimum_jaccard):
    distance = np.clip(1.0 - matrix, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)
    tree = linkage(condensed, method="average")
    labels = fcluster(tree, t=1.0 - minimum_jaccard, criterion="distance")
    return canonical_partition(labels)


def partition_text(partition):
    return "|".join(",".join(SAMPLE_ORDER[index] for index in group) for group in partition)


def group_signature(group):
    return "".join("1" if index in group else "0" for index in range(len(SAMPLE_ORDER)))


def merge_partitions(partitions, window_bp, path_length):
    blocks = []
    for index, partition in enumerate(partitions):
        start = index * window_bp
        end = min((index + 1) * window_bp, path_length)
        if blocks and blocks[-1]["partition"] == partition:
            blocks[-1]["end"] = end
            blocks[-1]["windows"] += 1
        else:
            blocks.append({"start": start, "end": end, "windows": 1, "partition": partition})
    return blocks


def write_pairwise(method, matrices, window_bp, path_length, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "window_start", "window_end", "haplotype_1",
                         "haplotype_2", "bp_weighted_jaccard"])
        for window, matrix in enumerate(matrices):
            start = window * window_bp
            end = min((window + 1) * window_bp, path_length)
            for i in range(len(SAMPLE_ORDER)):
                for j in range(i + 1, len(SAMPLE_ORDER)):
                    writer.writerow([method, start, end, SAMPLE_ORDER[i], SAMPLE_ORDER[j], matrix[i, j]])


def write_window_clusters(method, partitions, window_bp, path_length, output):
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "start", "end", "number_of_groups", "partition",
                         *[f"{sample}_group" for sample in SAMPLE_ORDER]])
        for index, partition in enumerate(partitions):
            labels = {}
            for group_number, group in enumerate(partition, 1):
                signature = group_signature(group)
                for sample_index in group:
                    labels[sample_index] = signature
            writer.writerow([method, index * window_bp,
                             min((index + 1) * window_bp, path_length), len(partition),
                             partition_text(partition), *[labels[i] for i in range(len(SAMPLE_ORDER))]])


def write_blocks(method, blocks, output):
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["block_id", "method", "start", "end", "length", "windows",
                         "number_of_groups", "group_id", "group_signature", "haplotypes"])
        for block_index, block in enumerate(blocks, 1):
            block_id = f"{method}.SharingBlock_{block_index:04d}"
            for group_index, group in enumerate(block["partition"], 1):
                writer.writerow([
                    block_id, method, block["start"], block["end"],
                    block["end"] - block["start"], block["windows"],
                    len(block["partition"]), f"{block_id}.Group_{group_index}",
                    group_signature(group), ",".join(SAMPLE_ORDER[i] for i in group),
                ])


def signature_color(signature):
    value = int(signature, 2)
    return plt.get_cmap("turbo")((value * 2654435761 % 997) / 996.0)


def partition_image(partitions):
    image = np.zeros((len(SAMPLE_ORDER), len(partitions), 4), dtype=float)
    for column, partition in enumerate(partitions):
        for group in partition:
            color = signature_color(group_signature(group))
            for row in group:
                image[row, column] = color
    return image


def plot_group_map(method_partitions, output, window_bp, path_length, title, start=0, end=None):
    if end is None:
        end = path_length
    first = start // window_bp
    last = min(len(next(iter(method_partitions.values()))), (end + window_bp - 1) // window_bp)
    methods = list(method_partitions)
    fig, axes = plt.subplots(len(methods), 1, figsize=(18, 2.3 * len(methods)), sharex=True)
    for ax, method in zip(axes, methods):
        subset = method_partitions[method][first:last]
        image = partition_image(subset)
        ax.imshow(image, aspect="auto", interpolation="nearest",
                  extent=[first * window_bp / 1e6, min(last * window_bp, path_length) / 1e6,
                          len(SAMPLE_ORDER) - 0.5, -0.5])
        for index in range(1, len(subset)):
            if subset[index] != subset[index - 1]:
                ax.axvline((first + index) * window_bp / 1e6,
                           color="white", linewidth=0.4, alpha=0.8)
        ax.set_yticks(range(len(SAMPLE_ORDER)), SAMPLE_ORDER, fontsize=7)
        ax.set_ylabel(method, rotation=0, ha="right", va="center", labelpad=10, fontsize=9)
    axes[-1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    fig.suptitle(title)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_group_counts(method_partitions, output, window_bp, path_length):
    fig, ax = plt.subplots(figsize=(16, 4.5))
    for method, partitions in method_partitions.items():
        x = np.arange(len(partitions)) * window_bp / 1e6
        ax.step(x, [len(partition) for partition in partitions], where="post", label=method)
    ax.set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    ax.set_ylabel("haplotype groups per 10 kb window")
    ax.set_xlim(0, path_length / 1e6)
    ax.set_ylim(0.8, len(SAMPLE_ORDER) + 0.2)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", action="append", required=True, help="METHOD=GFA")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--window-bp", type=int, default=10_000)
    parser.add_argument("--thresholds", default="0.7,0.8,0.9")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = [float(value) for value in args.thresholds.split(",")]
    matrices_by_method = {}
    path_length = None
    for item in args.graph:
        method, gfa = item.split("=", 1)
        length, matrices = load_window_jaccard(method, Path(gfa), args.window_bp)
        if path_length is None:
            path_length = length
        elif length != path_length:
            raise ValueError(f"CHM13 path length mismatch for {method}")
        matrices_by_method[method] = matrices
        write_pairwise(method, matrices, args.window_bp, path_length,
                       args.output_dir / method / "pairwise_jaccard_10kb.tsv")

    summary_rows = []
    representative_partitions = {}
    for threshold in thresholds:
        threshold_dir = args.output_dir / f"jaccard_{threshold:.1f}"
        threshold_dir.mkdir(parents=True, exist_ok=True)
        for method, matrices in matrices_by_method.items():
            partitions = [cluster_matrix(matrix, threshold) for matrix in matrices]
            blocks = merge_partitions(partitions, args.window_bp, path_length)
            write_window_clusters(method, partitions, args.window_bp, path_length,
                                  threshold_dir / f"{method}.window_clusters.tsv")
            write_blocks(method, blocks, threshold_dir / f"{method}.sharing_blocks.tsv")
            group_counts = np.array([len(partition) for partition in partitions])
            lengths = np.array([block["end"] - block["start"] for block in blocks])
            summary_rows.append({
                "method": method, "minimum_jaccard": threshold, "windows": len(partitions),
                "sharing_blocks": len(blocks), "median_block_bp": float(np.median(lengths)),
                "max_block_bp": int(lengths.max()), "mean_groups_per_window": group_counts.mean(),
                "max_groups_per_window": int(group_counts.max()),
                "single_group_windows": int((group_counts == 1).sum()),
            })
            if threshold == 0.8:
                representative_partitions[method] = partitions

    fields = list(summary_rows[0])
    with open(args.output_dir / "clustering_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    with open(args.output_dir / "README.txt", "w") as handle:
        handle.write("Jaccard is computed over CHM13-path node bp in each window.\n")
        handle.write("It measures shared retention of CHM13-anchored graph sequence and does not prove that paths use the same non-reference branch.\n")
        handle.write("If both paths visit none of the CHM13 nodes in a window, their empty-set Jaccard is defined as 1.0.\n")
    plot_group_map(representative_partitions, args.output_dir / "haplotype_group_sharing_map_jaccard0.8.png",
                   args.window_bp, path_length,
                   "10 kb pairwise sharing groups (bp-weighted Jaccard, average linkage, threshold 0.8)")
    plot_group_map(representative_partitions, args.output_dir / "haplotype_group_sharing_map_zoom_1p9_2p2Mb.png",
                   args.window_bp, path_length,
                   "Pairwise sharing groups around 1.9-2.2 Mb", 1_900_000, 2_200_000)
    plot_group_counts(representative_partitions, args.output_dir / "haplotype_group_count_jaccard0.8.png",
                      args.window_bp, path_length)


if __name__ == "__main__":
    main()
