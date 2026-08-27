#!/usr/bin/env python3
"""Window and segment CHM13-anchored haplotype-sharing profiles."""

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_gfa_path_sharing import SAMPLE_ORDER, parse_gfa, sample_name


def load_profiles(method, gfa, window_sizes):
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
    path_length = sum(lengths[node] for node in sample_paths["chm13v2.0"])
    profiles = {
        size: np.zeros(((path_length + size - 1) // size, len(SAMPLE_ORDER)), dtype=float)
        for size in window_sizes
    }
    denominators = {size: np.zeros(matrix.shape[0], dtype=float) for size, matrix in profiles.items()}
    pos = 0
    for node in sample_paths["chm13v2.0"]:
        end = pos + lengths[node]
        signature = np.array([counts[sample][node] > 0 for sample in SAMPLE_ORDER], dtype=float)
        for size, matrix in profiles.items():
            first = pos // size
            last = (end - 1) // size
            for index in range(first, last + 1):
                overlap = max(0, min(end, (index + 1) * size) - max(pos, index * size))
                matrix[index] += overlap * signature
                denominators[size][index] += overlap
        pos = end
    for size in window_sizes:
        profiles[size] /= denominators[size][:, None]
    return path_length, profiles


def write_windows(method, matrix, window_bp, path_length, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "start", "end", "mean_support", *SAMPLE_ORDER])
        for index, vector in enumerate(matrix):
            writer.writerow([method, index * window_bp, min((index + 1) * window_bp, path_length),
                             vector.mean(), *vector])


def cosine(a, b):
    denominator = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denominator) if denominator else 0.0


def cosine_blocks(matrix, threshold):
    blocks = []
    start = 0
    total = matrix[0].copy()
    count = 1
    for index in range(1, len(matrix)):
        centroid = total / count
        if cosine(centroid, matrix[index]) >= threshold:
            total += matrix[index]
            count += 1
        else:
            blocks.append((start, index, total / count))
            start = index
            total = matrix[index].copy()
            count = 1
    blocks.append((start, len(matrix), total / count))
    return blocks


def prefix_stats(matrix):
    sums = np.vstack([np.zeros(matrix.shape[1]), np.cumsum(matrix, axis=0)])
    squares = np.concatenate([[0.0], np.cumsum(np.sum(matrix * matrix, axis=1))])
    return sums, squares


def segment_cost(sums, squares, start, end):
    n = end - start
    vector_sum = sums[end] - sums[start]
    return float(squares[end] - squares[start] - np.dot(vector_sum, vector_sum) / n)


def binary_segmentation(matrix, penalty_scale, min_windows=2):
    sums, squares = prefix_stats(matrix)
    differences = np.diff(matrix, axis=0)
    noise = float(np.median(np.sum(differences * differences, axis=1)) / 2.0)
    noise = max(noise, 1e-6)
    penalty = penalty_scale * noise * matrix.shape[1] * math.log(len(matrix))
    segments = [(0, len(matrix))]
    changed = True
    while changed:
        changed = False
        output = []
        for start, end in segments:
            if end - start < 2 * min_windows:
                output.append((start, end))
                continue
            unsplit = segment_cost(sums, squares, start, end)
            best_gain = -1.0
            best_split = None
            for split in range(start + min_windows, end - min_windows + 1):
                gain = unsplit - segment_cost(sums, squares, start, split) - segment_cost(sums, squares, split, end)
                if gain > best_gain:
                    best_gain = gain
                    best_split = split
            if best_split is not None and best_gain > penalty:
                output.extend([(start, best_split), (best_split, end)])
                changed = True
            else:
                output.append((start, end))
        segments = output
    return [(start, end, matrix[start:end].mean(axis=0)) for start, end in segments], noise, penalty


def write_blocks(method, blocks, window_bp, path_length, output, metadata):
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["block_id", "method", "start", "end", "length", "windows",
                         "mean_support", *SAMPLE_ORDER])
        for index, (first, last, vector) in enumerate(blocks, 1):
            start = first * window_bp
            end = min(last * window_bp, path_length)
            writer.writerow([f"{method}.Block_{index:04d}", method, start, end, end - start,
                             last - first, vector.mean(), *vector])
    with open(output.with_suffix(".meta.tsv"), "w") as handle:
        for key, value in metadata.items():
            handle.write(f"{key}\t{value}\n")


def block_vector_track(blocks, n_windows):
    track = np.zeros((n_windows, len(SAMPLE_ORDER)), dtype=float)
    for start, end, vector in blocks:
        track[start:end] = vector
    return track


def plot_window_heatmap(method_matrices, output, window_bp, path_length, title):
    methods = list(method_matrices)
    fig, axes = plt.subplots(len(methods), 1, figsize=(18, 2.35 * len(methods)), sharex=True)
    for ax, method in zip(axes, methods):
        matrix = method_matrices[method]
        image = ax.imshow(matrix.T, aspect="auto", interpolation="nearest", vmin=0, vmax=1,
                          cmap="viridis", extent=[0, path_length / 1e6,
                                                  len(SAMPLE_ORDER) - 0.5, -0.5])
        ax.set_yticks(range(len(SAMPLE_ORDER)), SAMPLE_ORDER, fontsize=7)
        ax.set_ylabel(method, rotation=0, ha="right", va="center", labelpad=10, fontsize=9)
    axes[-1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    fig.colorbar(image, ax=axes, label="fraction of window bp shared with path", pad=0.01)
    fig.suptitle(title)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_segmentation(method_tracks, method_blocks, output, window_bp, path_length, title):
    methods = list(method_tracks)
    fig, axes = plt.subplots(len(methods), 1, figsize=(18, 2.35 * len(methods)), sharex=True)
    for ax, method in zip(axes, methods):
        track = method_tracks[method]
        image = ax.imshow(track.T, aspect="auto", interpolation="nearest", vmin=0, vmax=1,
                          cmap="viridis", extent=[0, path_length / 1e6,
                                                  len(SAMPLE_ORDER) - 0.5, -0.5])
        for _, end, _ in method_blocks[method][:-1]:
            ax.axvline(end * window_bp / 1e6, color="white", linewidth=0.35, alpha=0.8)
        ax.set_yticks(range(len(SAMPLE_ORDER)), SAMPLE_ORDER, fontsize=7)
        ax.set_ylabel(method, rotation=0, ha="right", va="center", labelpad=10, fontsize=9)
    axes[-1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    fig.colorbar(image, ax=axes, label="block mean sharing fraction", pad=0.01)
    fig.suptitle(title)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def summarize(method, approach, parameter, blocks, window_bp, path_length):
    lengths = np.array([min(end * window_bp, path_length) - start * window_bp
                        for start, end, _ in blocks])
    return {
        "method": method, "approach": approach, "parameter": parameter,
        "blocks": len(blocks), "median_block_bp": float(np.median(lengths)),
        "mean_block_bp": float(np.mean(lengths)), "max_block_bp": int(lengths.max()),
        "blocks_per_mb": len(blocks) / (path_length / 1e6),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", action="append", required=True, help="METHOD=GFA")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    window_sizes = [1000, 5000, 10000]
    profiles = {}
    path_length = None
    for item in args.graph:
        method, gfa = item.split("=", 1)
        length, method_profiles = load_profiles(method, Path(gfa), window_sizes)
        if path_length is None:
            path_length = length
        elif length != path_length:
            raise ValueError(f"CHM13 path length mismatch for {method}")
        profiles[method] = method_profiles
        for size, matrix in method_profiles.items():
            write_windows(method, matrix, size, path_length,
                          args.output_dir / method / "method_A" / f"windows_{size}bp.tsv")

    for size in window_sizes:
        plot_window_heatmap({method: data[size] for method, data in profiles.items()},
                            args.output_dir / f"method_A_windows_{size}bp.png", size, path_length,
                            f"Method A: {size // 1000} kb window sharing profiles")

    base_window = 5000
    cosine_thresholds = [0.95, 0.98, 0.99, 0.995]
    penalty_scales = [1.0, 3.0, 5.0]
    summaries = []
    b_representative = {}
    c_representative = {}
    b_tracks = {}
    c_tracks = {}
    for method, data in profiles.items():
        matrix = data[base_window]
        for threshold in cosine_thresholds:
            blocks = cosine_blocks(matrix, threshold)
            write_blocks(method, blocks, base_window, path_length,
                         args.output_dir / method / "method_B" / f"cosine_{threshold:.3f}.tsv",
                         {"approach": "adjacent_window_cosine", "threshold": threshold,
                          "window_bp": base_window})
            summaries.append(summarize(method, "B_cosine", threshold, blocks,
                                       base_window, path_length))
            if threshold == 0.99:
                b_representative[method] = blocks
                b_tracks[method] = block_vector_track(blocks, len(matrix))
        for scale in penalty_scales:
            blocks, noise, penalty = binary_segmentation(matrix, scale)
            write_blocks(method, blocks, base_window, path_length,
                         args.output_dir / method / "method_C" / f"binary_bic_{scale:g}.tsv",
                         {"approach": "multivariate_l2_binary_segmentation",
                          "penalty_scale": scale, "noise_estimate": noise,
                          "penalty": penalty, "window_bp": base_window, "min_windows": 2})
            summaries.append(summarize(method, "C_binary_BIC", scale, blocks,
                                       base_window, path_length))
            if scale == 3.0:
                c_representative[method] = blocks
                c_tracks[method] = block_vector_track(blocks, len(matrix))

    fields = list(summaries[0])
    with open(args.output_dir / "segmentation_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)
    plot_segmentation(b_tracks, b_representative, args.output_dir / "method_B_cosine_0.99.png",
                      base_window, path_length,
                      "Method B: adjacent 5 kb windows merged by cosine similarity >= 0.99")
    plot_segmentation(c_tracks, c_representative, args.output_dir / "method_C_binary_BIC3.png",
                      base_window, path_length,
                      "Method C: multivariate change-point segmentation (binary L2, BIC scale 3)")


if __name__ == "__main__":
    main()
