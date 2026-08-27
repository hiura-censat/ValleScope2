#!/usr/bin/env python3
"""Project full graph paths to CHM13 and build non-reference-aware sharing blocks."""

import argparse
import bisect
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_gfa_path_sharing import SAMPLE_ORDER, sample_name
from cluster_window_pairwise_sharing import (
    cluster_matrix, merge_partitions, partition_text, group_signature,
    plot_group_counts, plot_group_map,
)


def parse_oriented_gfa(path):
    lengths = {}
    paths = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S":
                if fields[2] != "*":
                    lengths[fields[1]] = len(fields[2])
                else:
                    tag = next(x for x in fields[3:] if x.startswith("LN:i:"))
                    lengths[fields[1]] = int(tag[5:])
            elif fields[0] == "P":
                paths[fields[1]] = [(token[:-1], token[-1]) for token in fields[2].split(",") if token]
    return lengths, paths


def cumulative_positions(steps, lengths):
    positions = np.zeros(len(steps) + 1, dtype=np.int64)
    positions[1:] = np.cumsum([lengths[node] for node, _ in steps], dtype=np.int64)
    return positions


def longest_increasing_anchors(pairs):
    tails = []
    tail_indices = []
    previous = [-1] * len(pairs)
    for index, (_, ref_index) in enumerate(pairs):
        place = bisect.bisect_left(tails, ref_index)
        if place == len(tails):
            tails.append(ref_index)
            tail_indices.append(index)
        else:
            tails[place] = ref_index
            tail_indices[place] = index
        if place:
            previous[index] = tail_indices[place - 1]
    if not tail_indices:
        return []
    selected = []
    index = tail_indices[-1]
    while index >= 0:
        selected.append(pairs[index])
        index = previous[index]
    return selected[::-1]


def build_anchor_chain(ref_steps, sample_steps, ref_positions, sample_positions):
    ref_count = Counter(ref_steps)
    sample_count = Counter(sample_steps)
    ref_unique = {step: index for index, step in enumerate(ref_steps) if ref_count[step] == 1}
    pairs = [(index, ref_unique[step]) for index, step in enumerate(sample_steps)
             if sample_count[step] == 1 and step in ref_unique]
    anchors = longest_increasing_anchors(pairs)
    coordinates = [(0.0, 0.0, -1, -1)]
    for sample_index, ref_index in anchors:
        sample_bp = float(sample_positions[sample_index])
        ref_bp = float(ref_positions[ref_index])
        if sample_bp > coordinates[-1][0] and ref_bp > coordinates[-1][1]:
            coordinates.append((sample_bp, ref_bp, sample_index, ref_index))
    sample_end = float(sample_positions[-1])
    ref_end = float(ref_positions[-1])
    if sample_end > coordinates[-1][0] and ref_end > coordinates[-1][1]:
        coordinates.append((sample_end, ref_end, len(sample_steps), len(ref_steps)))
    return coordinates, len(pairs)


def distribute_weight(windows, feature, weight, projected_start, projected_end, window_bp, path_length):
    projected_start = min(max(projected_start, 0.0), path_length)
    projected_end = min(max(projected_end, 0.0), path_length)
    if projected_end < projected_start:
        projected_start, projected_end = projected_end, projected_start
    if projected_end == projected_start:
        index = min(len(windows) - 1, int(projected_start // window_bp))
        windows[index][feature] += weight
        return
    span = projected_end - projected_start
    first = min(len(windows) - 1, int(projected_start // window_bp))
    last = min(len(windows) - 1, int((projected_end - 1e-9) // window_bp))
    assigned = 0.0
    for index in range(first, last + 1):
        overlap = max(0.0, min(projected_end, (index + 1) * window_bp) -
                      max(projected_start, index * window_bp))
        value = weight * overlap / span
        windows[index][feature] += value
        assigned += value
    if abs(assigned - weight) > max(1e-6, weight * 1e-8):
        raise ValueError("projected node weight was not conserved")


def project_reference(ref_steps, ref_positions, lengths, window_bp):
    path_length = int(ref_positions[-1])
    windows = [Counter() for _ in range((path_length + window_bp - 1) // window_bp)]
    for index, feature in enumerate(ref_steps):
        distribute_weight(windows, feature, lengths[feature[0]],
                          float(ref_positions[index]), float(ref_positions[index + 1]),
                          window_bp, path_length)
    return windows


def project_sample(sample, steps, positions, anchors, lengths, ref_node_ids,
                   window_bp, path_length, branch_writer, anchor_writer):
    windows = [Counter() for _ in range((path_length + window_bp - 1) // window_bp)]
    anchor_sample_bp = [anchor[0] for anchor in anchors]
    interval = 0
    low_confidence_bp = 0.0
    branch = None
    branch_index = 0
    for step_index, feature in enumerate(steps):
        midpoint = (positions[step_index] + positions[step_index + 1]) / 2.0
        while interval + 1 < len(anchors) - 1 and midpoint >= anchors[interval + 1][0]:
            interval += 1
        sample_left, ref_left, _, _ = anchors[interval]
        sample_right, ref_right, _, _ = anchors[interval + 1]
        sample_span = sample_right - sample_left
        ref_span = ref_right - ref_left
        if sample_span <= 0 or ref_span <= 0:
            raise ValueError(f"{sample}: non-positive anchor interval")
        scale = ref_span / sample_span
        projected_start = ref_left + (positions[step_index] - sample_left) * scale
        projected_end = ref_left + (positions[step_index + 1] - sample_left) * scale
        length = lengths[feature[0]]
        distribute_weight(windows, feature, length, projected_start, projected_end,
                          window_bp, path_length)
        confidence = "low" if ref_span > 100_000 else "high"
        if confidence == "low":
            low_confidence_bp += length
        is_nonreference = feature[0] not in ref_node_ids
        if is_nonreference:
            if branch is None:
                branch_index += 1
                branch = {
                    "start_step": step_index, "end_step": step_index + 1,
                    "sample_start": int(positions[step_index]),
                    "sample_end": int(positions[step_index + 1]),
                    "projected_start": projected_start, "projected_end": projected_end,
                    "bp": length, "nodes": {feature}, "confidence": confidence,
                    "anchor_ref_start": ref_left, "anchor_ref_end": ref_right,
                }
            else:
                branch["end_step"] = step_index + 1
                branch["sample_end"] = int(positions[step_index + 1])
                branch["projected_end"] = projected_end
                branch["bp"] += length
                branch["nodes"].add(feature)
                if confidence == "low":
                    branch["confidence"] = "low"
        elif branch is not None:
            write_branch(branch_writer, sample, branch_index, branch)
            branch = None
    if branch is not None:
        write_branch(branch_writer, sample, branch_index, branch)
    for index, anchor in enumerate(anchors):
        if index + 1 == len(anchors):
            break
        anchor_writer.writerow([
            sample, index, int(anchor[0]), int(anchors[index + 1][0]),
            int(anchor[1]), int(anchors[index + 1][1]),
            int(anchors[index + 1][0] - anchor[0]),
            int(anchors[index + 1][1] - anchor[1]),
            "low" if anchors[index + 1][1] - anchor[1] > 100_000 else "high",
        ])
    return windows, low_confidence_bp, branch_index


def write_branch(writer, sample, branch_index, branch):
    writer.writerow([
        f"{sample}.Branch_{branch_index:06d}", sample,
        branch["start_step"], branch["end_step"], branch["sample_start"], branch["sample_end"],
        branch["sample_end"] - branch["sample_start"],
        int(round(branch["projected_start"])), int(round(branch["projected_end"])),
        int(round(branch["anchor_ref_start"])), int(round(branch["anchor_ref_end"])),
        branch["bp"], len(branch["nodes"]), branch["confidence"],
    ])


def multiset_jaccard(a, b):
    keys = a.keys() | b.keys()
    intersection = sum(min(a[key], b[key]) for key in keys)
    union = sum(max(a[key], b[key]) for key in keys)
    return intersection / union if union else 1.0


def calculate_matrices(projected):
    samples = SAMPLE_ORDER
    n_windows = len(projected[samples[0]])
    matrices = np.ones((n_windows, len(samples), len(samples)), dtype=float)
    for window in range(n_windows):
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                value = multiset_jaccard(projected[samples[i]][window], projected[samples[j]][window])
                matrices[window, i, j] = matrices[window, j, i] = value
    return matrices


def write_pairwise(method, matrices, window_bp, path_length, output):
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "window_start", "window_end", "haplotype_1",
                         "haplotype_2", "full_graph_multiset_jaccard"])
        for window, matrix in enumerate(matrices):
            for i in range(len(SAMPLE_ORDER)):
                for j in range(i + 1, len(SAMPLE_ORDER)):
                    writer.writerow([method, window * window_bp,
                                     min((window + 1) * window_bp, path_length),
                                     SAMPLE_ORDER[i], SAMPLE_ORDER[j], matrix[i, j]])


def write_clusters(method, partitions, blocks, window_bp, path_length, directory):
    with open(directory / "full_window_clusters.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "start", "end", "number_of_groups", "partition"])
        for index, partition in enumerate(partitions):
            writer.writerow([method, index * window_bp, min((index + 1) * window_bp, path_length),
                             len(partition), partition_text(partition)])
    with open(directory / "full_sharing_blocks.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["block_id", "method", "start", "end", "length", "windows",
                         "number_of_groups", "group_id", "signature", "haplotypes"])
        for block_index, block in enumerate(blocks, 1):
            block_id = f"{method}.FullBlock_{block_index:04d}"
            for group_index, group in enumerate(block["partition"], 1):
                writer.writerow([
                    block_id, method, block["start"], block["end"], block["end"] - block["start"],
                    block["windows"], len(block["partition"]), f"{block_id}.Group_{group_index}",
                    group_signature(group), ",".join(SAMPLE_ORDER[i] for i in group),
                ])


def load_phase_a_partitions(path):
    partitions = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups = []
            for text in row["partition"].split("|"):
                members = tuple(SAMPLE_ORDER.index(sample) for sample in text.split(","))
                groups.append(members)
            partitions.append(tuple(sorted(groups)))
    return partitions


def coclustering_agreement(a, b):
    def group_map(partition):
        return {sample: group for group, members in enumerate(partition) for sample in members}
    left, right = group_map(a), group_map(b)
    agreements = 0
    total = 0
    for i in range(len(SAMPLE_ORDER)):
        for j in range(i + 1, len(SAMPLE_ORDER)):
            agreements += ((left[i] == left[j]) == (right[i] == right[j]))
            total += 1
    return agreements / total


def plot_phase_comparison(phase_a, phase_b, output, window_bp, path_length):
    methods = list(phase_b)
    fig, axes = plt.subplots(len(methods), 1, figsize=(16, 1.8 * len(methods)), sharex=True)
    for ax, method in zip(axes, methods):
        values = [coclustering_agreement(a, b) for a, b in zip(phase_a[method], phase_b[method])]
        x = np.arange(len(values)) * window_bp / 1e6
        ax.fill_between(x, values, 1, step="post", color="#2b8cbe", alpha=0.85)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel(method, rotation=0, ha="right", va="center")
    axes[-1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    fig.suptitle("Phase A vs Phase B pairwise co-clustering agreement")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_phase_comparison(phase_a, phase_b, output, window_bp, path_length):
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "window_start", "window_end",
                         "phase_a_phase_b_coclustering_agreement"])
        for method in phase_b:
            for index, (a, b) in enumerate(zip(phase_a[method], phase_b[method])):
                writer.writerow([
                    method,
                    index * window_bp,
                    min((index + 1) * window_bp, path_length),
                    coclustering_agreement(a, b),
                ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", action="append", required=True, help="METHOD=GFA")
    parser.add_argument("--phase-a-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--window-bp", type=int, default=10_000)
    parser.add_argument("--minimum-jaccard", type=float, default=0.8)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    phase_b_partitions = {}
    phase_a_partitions = {}
    summary = []
    path_length = None
    for item in args.graph:
        method, gfa = item.split("=", 1)
        method_dir = args.output_dir / method
        method_dir.mkdir(parents=True, exist_ok=True)
        lengths, paths = parse_oriented_gfa(Path(gfa))
        sample_paths = {sample_name(name): steps for name, steps in paths.items()
                        if sample_name(name) in SAMPLE_ORDER}
        ref_steps = sample_paths["chm13v2.0"]
        ref_positions = cumulative_positions(ref_steps, lengths)
        current_path_length = int(ref_positions[-1])
        if path_length is None:
            path_length = current_path_length
        elif current_path_length != path_length:
            raise ValueError(f"CHM13 length mismatch for {method}")
        projected = {"chm13v2.0": project_reference(ref_steps, ref_positions, lengths,
                                                    args.window_bp)}
        ref_node_ids = {node for node, _ in ref_steps}
        branch_path = method_dir / "nonreference_branches.tsv"
        anchor_path = method_dir / "anchor_intervals.tsv"
        with open(branch_path, "w") as branch_handle, open(anchor_path, "w") as anchor_handle:
            branch_writer = csv.writer(branch_handle, delimiter="\t", lineterminator="\n")
            anchor_writer = csv.writer(anchor_handle, delimiter="\t", lineterminator="\n")
            branch_writer.writerow(["branch_id", "sample", "start_step", "end_step",
                                    "sample_start", "sample_end", "sample_bp",
                                    "projected_start", "projected_end", "anchor_ref_start",
                                    "anchor_ref_end", "nonreference_bp", "unique_oriented_nodes",
                                    "confidence"])
            anchor_writer.writerow(["sample", "interval", "sample_start", "sample_end",
                                    "ref_start", "ref_end", "sample_span", "ref_span", "confidence"])
            for sample in SAMPLE_ORDER[1:]:
                steps = sample_paths[sample]
                positions = cumulative_positions(steps, lengths)
                anchors, candidate_count = build_anchor_chain(ref_steps, steps, ref_positions, positions)
                windows, low_bp, branch_count = project_sample(
                    sample, steps, positions, anchors, lengths, ref_node_ids,
                    args.window_bp, path_length, branch_writer, anchor_writer)
                projected[sample] = windows
                projected_bp = sum(sum(counter.values()) for counter in windows)
                summary.append({
                    "method": method, "sample": sample, "path_bp": int(positions[-1]),
                    "projected_bp": projected_bp, "candidate_unique_anchors": candidate_count,
                    "selected_monotonic_anchors": len(anchors) - 2,
                    "anchor_intervals": len(anchors) - 1, "nonreference_branches": branch_count,
                    "low_confidence_projected_bp": low_bp,
                })
        matrices = calculate_matrices(projected)
        write_pairwise(method, matrices, args.window_bp, path_length,
                       method_dir / "full_pairwise_jaccard_10kb.tsv")
        partitions = [cluster_matrix(matrix, args.minimum_jaccard) for matrix in matrices]
        blocks = merge_partitions(partitions, args.window_bp, path_length)
        write_clusters(method, partitions, blocks, args.window_bp, path_length, method_dir)
        phase_b_partitions[method] = partitions
        phase_a_path = args.phase_a_root / "jaccard_0.8" / f"{method}.window_clusters.tsv"
        phase_a_partitions[method] = load_phase_a_partitions(phase_a_path)

    fields = list(summary[0])
    with open(args.output_dir / "projection_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    with open(args.output_dir / "README.txt", "w") as handle:
        handle.write("Full sample path node-orientation multisets are projected to CHM13 with unique shared occurrence anchors.\n")
        handle.write("The longest increasing anchor chain defines the syntenic projection; intervals >100 kb are marked low confidence.\n")
        handle.write("Jaccard uses node bp with copy multiplicity: sum(min weights) / sum(max weights).\n")
    plot_group_map(phase_b_partitions, args.output_dir / "full_haplotype_group_map.png",
                   args.window_bp, path_length,
                   "Phase B full-graph haplotype groups (10 kb, multiset Jaccard >= 0.8)")
    plot_group_counts(phase_b_partitions, args.output_dir / "full_haplotype_group_counts.png",
                      args.window_bp, path_length)
    plot_phase_comparison(phase_a_partitions, phase_b_partitions,
                          args.output_dir / "phaseA_vs_phaseB_agreement.png",
                          args.window_bp, path_length)
    write_phase_comparison(phase_a_partitions, phase_b_partitions,
                           args.output_dir / "phaseA_vs_phaseB_agreement.tsv",
                           args.window_bp, path_length)


if __name__ == "__main__":
    main()
