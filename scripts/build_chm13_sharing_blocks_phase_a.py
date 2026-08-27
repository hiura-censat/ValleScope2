#!/usr/bin/env python3
"""Build CHM13-anchored haplotype-sharing blocks from path-embedded GFAs."""

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_gfa_path_sharing import SAMPLE_ORDER, parse_gfa, sample_name


CLASS_COLORS = {
    "core": "#1b7837",
    "common": "#5aae61",
    "shell": "#f1b65c",
    "private": "#b2182b",
}


def block_class(support):
    if support == len(SAMPLE_ORDER):
        return "core"
    if support >= 8:
        return "common"
    if support >= 2:
        return "shell"
    return "private"


def signature_text(signature):
    return "".join("1" if value else "0" for value in signature)


def load_chm13_nodes(method, gfa):
    lengths, paths = parse_gfa(gfa)
    sample_paths = {}
    for path_name, steps in paths.items():
        sample = sample_name(path_name)
        if sample in SAMPLE_ORDER:
            if sample in sample_paths:
                raise ValueError(f"{method}: duplicate path for {sample}")
            sample_paths[sample] = steps
    missing = [sample for sample in SAMPLE_ORDER if sample not in sample_paths]
    if missing:
        raise ValueError(f"{method}: missing paths: {','.join(missing)}")
    counts = {sample: Counter(sample_paths[sample]) for sample in SAMPLE_ORDER}
    nodes = []
    pos = 0
    for node_id in sample_paths["chm13v2.0"]:
        length = lengths[node_id]
        signature = tuple(counts[sample][node_id] > 0 for sample in SAMPLE_ORDER)
        traversals = tuple(counts[sample][node_id] for sample in SAMPLE_ORDER)
        nodes.append({
            "node": node_id, "start": pos, "end": pos + length,
            "length": length, "signature": signature, "traversals": traversals,
        })
        pos += length
    return nodes


def exact_blocks(nodes):
    blocks = []
    for node in nodes:
        if blocks and blocks[-1]["signature"] == node["signature"]:
            blocks[-1]["end"] = node["end"]
            blocks[-1]["length"] += node["length"]
            blocks[-1]["node_count"] += 1
        else:
            blocks.append({
                "start": node["start"], "end": node["end"], "length": node["length"],
                "signature": node["signature"], "node_count": 1,
            })
    return blocks


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def smooth_blocks(raw_blocks, threshold):
    blocks = [dict(block) for block in raw_blocks]
    changed = True
    merges = 0
    while changed:
        changed = False
        output = []
        i = 0
        while i < len(blocks):
            if i + 2 < len(blocks):
                left, middle, right = blocks[i:i + 3]
                if (left["signature"] == right["signature"] and
                        middle["length"] < threshold and
                        hamming(middle["signature"], left["signature"]) <= 1):
                    output.append({
                        "start": left["start"], "end": right["end"],
                        "length": right["end"] - left["start"],
                        "signature": left["signature"],
                        "node_count": left["node_count"] + middle["node_count"] + right["node_count"],
                    })
                    merges += 1
                    changed = True
                    i += 3
                    continue
            current = dict(blocks[i])
            if output and output[-1]["signature"] == current["signature"]:
                output[-1]["end"] = current["end"]
                output[-1]["length"] = output[-1]["end"] - output[-1]["start"]
                output[-1]["node_count"] += current["node_count"]
            else:
                output.append(current)
            i += 1
        blocks = output
    return blocks, merges


def write_nodes(nodes, output):
    with open(output, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["node_id", "start", "end", "length", "signature",
                         *[f"{sample}_presence" for sample in SAMPLE_ORDER],
                         *[f"{sample}_traversals" for sample in SAMPLE_ORDER]])
        for node in nodes:
            writer.writerow([node["node"], node["start"], node["end"], node["length"],
                             signature_text(node["signature"]),
                             *[int(x) for x in node["signature"]], *node["traversals"]])


def write_blocks(blocks, directory, method):
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / "sharing_blocks.tsv", "w") as tsv, open(directory / "sharing_blocks.bed", "w") as bed:
        writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")
        writer.writerow(["block_id", "method", "start", "end", "length", "signature",
                         "support_count", "class", "samples", "node_count"])
        for index, block in enumerate(blocks, 1):
            support = sum(block["signature"])
            cls = block_class(support)
            block_id = f"{method}.Block_{index:05d}"
            samples = ",".join(sample for sample, present in zip(SAMPLE_ORDER, block["signature"]) if present)
            writer.writerow([block_id, method, block["start"], block["end"], block["length"],
                             signature_text(block["signature"]), support, cls, samples,
                             block["node_count"]])
            bed.write(f"chm13_alphaSat\t{block['start']}\t{block['end']}\t{block_id}\t{support}\n")


def validate_blocks(nodes, blocks, expected_length):
    if not blocks or blocks[0]["start"] != 0 or blocks[-1]["end"] != expected_length:
        raise ValueError("blocks do not cover the full CHM13 path")
    for left, right in zip(blocks, blocks[1:]):
        if left["end"] != right["start"]:
            raise ValueError("block gap or overlap detected")
    if sum(block["length"] for block in blocks) != expected_length:
        raise ValueError("block lengths do not sum to the CHM13 path length")
    if sum(node["length"] for node in nodes) != expected_length:
        raise ValueError("node lengths do not sum to the CHM13 path length")


def summarize(method, label, blocks, path_length, merges):
    lengths = np.array([block["length"] for block in blocks], dtype=int)
    class_bp = Counter()
    class_count = Counter()
    for block in blocks:
        cls = block_class(sum(block["signature"]))
        class_bp[cls] += block["length"]
        class_count[cls] += 1
    row = {
        "method": method, "block_set": label, "blocks": len(blocks),
        "path_bp": path_length, "median_block_bp": float(np.median(lengths)),
        "mean_block_bp": float(np.mean(lengths)), "max_block_bp": int(lengths.max()),
        "blocks_per_mb": len(blocks) / (path_length / 1e6),
        "transitions_per_mb": (len(blocks) - 1) / (path_length / 1e6),
        "blocks_lt_1kb": int((lengths < 1000).sum()), "smoothing_merges": merges,
        "unique_signatures": len({block["signature"] for block in blocks}),
    }
    for cls in ["core", "common", "shell", "private"]:
        row[f"{cls}_blocks"] = class_count[cls]
        row[f"{cls}_bp"] = class_bp[cls]
    return row


def rasterize(blocks, width, start, end):
    support = np.zeros(width, dtype=float)
    membership = np.zeros((len(SAMPLE_ORDER), width), dtype=float)
    span = end - start
    for block in blocks:
        lo = max(start, block["start"])
        hi = min(end, block["end"])
        if lo >= hi:
            continue
        x0 = max(0, int((lo - start) / span * width))
        x1 = min(width, max(x0 + 1, int(np.ceil((hi - start) / span * width))))
        support[x0:x1] = sum(block["signature"])
        for row, present in enumerate(block["signature"]):
            membership[row, x0:x1] = int(present)
    return support, membership


def plot_membership(block_sets, output, start, end, title):
    methods = list(block_sets)
    width = 1800
    fig, axes = plt.subplots(len(methods), 1, figsize=(18, 2.45 * len(methods)), sharex=True)
    if len(methods) == 1:
        axes = [axes]
    cmap = plt.matplotlib.colors.ListedColormap(["#f2f2f2", "#202020"])
    for ax, method in zip(axes, methods):
        support, membership = rasterize(block_sets[method], width, start, end)
        ax.imshow(membership, aspect="auto", interpolation="nearest", cmap=cmap,
                  extent=[start / 1e6, end / 1e6, len(SAMPLE_ORDER) - 0.5, -0.5])
        ax.set_yticks(range(len(SAMPLE_ORDER)), SAMPLE_ORDER, fontsize=7)
        ax.set_ylabel(method, rotation=0, ha="right", va="center", labelpad=10, fontsize=9)
        twin = ax.inset_axes([0, 1.01, 1, 0.10])
        class_values = np.zeros((1, width, 3))
        for x, count in enumerate(support.astype(int)):
            color = CLASS_COLORS[block_class(count)]
            class_values[0, x] = tuple(int(color[i:i + 2], 16) / 255 for i in (1, 3, 5))
        twin.imshow(class_values, aspect="auto", interpolation="nearest")
        twin.set_axis_off()
    axes[-1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    fig.suptitle(title, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_block_classes(block_sets, output, path_length, title):
    methods = list(block_sets)
    fig, ax = plt.subplots(figsize=(17, 3.7))
    width = 2400
    image = np.zeros((len(methods), width, 3), dtype=float)
    for y, method in enumerate(methods):
        support, _ = rasterize(block_sets[method], width, 0, path_length)
        for x, count in enumerate(support.astype(int)):
            color = CLASS_COLORS[block_class(count)]
            image[y, x] = tuple(int(color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    ax.imshow(image, aspect="auto", interpolation="nearest",
              extent=[0, path_length / 1e6, len(methods) - 0.5, -0.5])
    ax.set_xlim(0, path_length / 1e6)
    ax.set_ylim(len(methods) - 0.5, -0.5)
    ax.set_yticks(range(len(methods)), methods)
    ax.set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    ax.set_title(title)
    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[x], label=x)
               for x in ["core", "common", "shell", "private"]]
    ax.legend(handles=handles, ncol=4, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.22))
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", action="append", required=True, help="METHOD=GFA")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--thresholds", default="100,500,1000,5000")
    args = parser.parse_args()
    thresholds = [int(x) for x in args.thresholds.split(",")]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    raw_by_method = {}
    smoothed_by_threshold = {threshold: {} for threshold in thresholds}
    summaries = []
    path_length = None
    for item in args.graph:
        method, gfa = item.split("=", 1)
        method_dir = args.output_dir / method
        method_dir.mkdir(parents=True, exist_ok=True)
        nodes = load_chm13_nodes(method, Path(gfa))
        current_length = nodes[-1]["end"]
        if path_length is None:
            path_length = current_length
        elif current_length != path_length:
            raise ValueError(f"CHM13 path length differs for {method}: {current_length} != {path_length}")
        write_nodes(nodes, method_dir / "chm13_node_membership.tsv")
        raw = exact_blocks(nodes)
        validate_blocks(nodes, raw, path_length)
        raw_by_method[method] = raw
        write_blocks(raw, method_dir / "raw", method)
        summaries.append(summarize(method, "raw", raw, path_length, 0))
        current_blocks = raw
        cumulative_merges = 0
        for threshold in sorted(thresholds):
            blocks, merges = smooth_blocks(current_blocks, threshold)
            cumulative_merges += merges
            validate_blocks(nodes, blocks, path_length)
            smoothed_by_threshold[threshold][method] = blocks
            write_blocks(blocks, method_dir / f"smoothed_{threshold}bp", method)
            summaries.append(summarize(method, f"smoothed_{threshold}bp", blocks,
                                       path_length, cumulative_merges))
            current_blocks = blocks

    fields = list(summaries[0])
    with open(args.output_dir / "block_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)

    plot_block_classes(raw_by_method, args.output_dir / "sharing_classes_raw.png", path_length,
                       "CHM13-anchored exact sharing blocks")
    plot_membership(raw_by_method, args.output_dir / "sharing_membership_raw.png", 0, path_length,
                    "CHM13-anchored exact haplotype membership")
    representative = smoothed_by_threshold[1000]
    plot_block_classes(representative, args.output_dir / "sharing_classes_smoothed_1kb.png",
                       path_length, "CHM13-anchored sharing blocks (A-B-A smoothing <1 kb)")
    plot_membership(representative, args.output_dir / "sharing_membership_smoothed_1kb.png",
                    0, path_length, "CHM13-anchored haplotype membership (smoothed <1 kb)")
    plot_membership(representative, args.output_dir / "sharing_membership_sparse_zoom.png",
                    1_900_000, 2_150_000,
                    "Haplotype membership around the ValleScope2 sparse interval (smoothed <1 kb)")


if __name__ == "__main__":
    main()
