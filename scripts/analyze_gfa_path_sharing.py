#!/usr/bin/env python3
"""Quantify sample sharing and copy-number structure from path-embedded GFAs."""

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SAMPLE_ORDER = [
    "chm13v2.0", "HG00146_hap1", "HG00358_hap2", "HG00544_hap2",
    "HG00741_hap1", "HG01114_hap1", "HG01255_hap2", "NA18534_hap2",
    "NA18974_hap2", "NA18982_hap1",
]


def sample_name(path_name):
    return path_name.split("|", 1)[0]


def parse_gfa(path):
    lengths = {}
    paths = {}
    with open(path) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                if fields[2] != "*":
                    lengths[fields[1]] = len(fields[2])
                else:
                    ln = next((x[5:] for x in fields[3:] if x.startswith("LN:i:")), None)
                    if ln is None:
                        raise ValueError(f"segment {fields[1]} has neither sequence nor LN tag")
                    lengths[fields[1]] = int(ln)
            elif fields[0] == "P":
                paths[fields[1]] = [x[:-1] for x in fields[2].split(",") if x]
    return lengths, paths


def support_class(n, total):
    if n == 0:
        return "unvisited"
    if n == 1:
        return "private"
    if n == total:
        return "core"
    return "shell"


def copy_bin(n):
    if n <= 1:
        return str(n)
    if n == 2:
        return "2"
    if n <= 5:
        return "3-5"
    if n <= 10:
        return "6-10"
    return ">10"


def analyze_graph(method, gfa, outdir, bin_bp):
    lengths, gfa_paths = parse_gfa(gfa)
    sample_paths = {}
    for name, steps in gfa_paths.items():
        sample = sample_name(name)
        if sample in SAMPLE_ORDER:
            if sample in sample_paths:
                raise ValueError(f"{method}: duplicate sample path {sample}")
            sample_paths[sample] = steps
    missing = [x for x in SAMPLE_ORDER if x not in sample_paths]
    if missing:
        raise ValueError(f"{method}: missing paths: {', '.join(missing)}")

    counts = {sample: Counter(sample_paths[sample]) for sample in SAMPLE_ORDER}
    visited = set().union(*(set(x) for x in counts.values()))
    graph_dir = outdir / method
    graph_dir.mkdir(parents=True, exist_ok=True)

    class_nodes = Counter()
    class_bp = Counter()
    copy_nodes = Counter()
    copy_bp = Counter()
    matrix_path = graph_dir / "node_sample_traversal_matrix.tsv.gz"
    with gzip.open(matrix_path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["node", "length", "sample_support", "class", *SAMPLE_ORDER])
        for node in sorted(visited, key=lambda x: (len(x), x)):
            values = [counts[s][node] for s in SAMPLE_ORDER]
            support = sum(x > 0 for x in values)
            cls = support_class(support, len(SAMPLE_ORDER))
            length = lengths[node]
            class_nodes[cls] += 1
            class_bp[cls] += length
            mb = copy_bin(max(values))
            copy_nodes[mb] += 1
            copy_bp[mb] += length
            writer.writerow([node, length, support, cls, *values])

    chm_steps = sample_paths["chm13v2.0"]
    chm_length = sum(lengths[x] for x in chm_steps)
    n_bins = (chm_length + bin_bp - 1) // bin_bp
    support_sum = np.zeros(n_bins, dtype=float)
    core_bp = np.zeros(n_bins, dtype=float)
    shell_bp = np.zeros(n_bins, dtype=float)
    private_bp = np.zeros(n_bins, dtype=float)
    covered_bp = np.zeros(n_bins, dtype=float)
    blocks = []
    pos = 0
    previous = None
    for node in chm_steps:
        length = lengths[node]
        node_counts = tuple(counts[s][node] for s in SAMPLE_ORDER)
        mask = tuple(i for i, value in enumerate(node_counts) if value > 0)
        cls = support_class(len(mask), len(SAMPLE_ORDER))
        end = pos + length
        if previous and previous[2] == mask:
            previous = (previous[0], end, mask, cls)
            blocks[-1] = previous
        else:
            previous = (pos, end, mask, cls)
            blocks.append(previous)
        left_bin = pos // bin_bp
        right_bin = (end - 1) // bin_bp
        for b in range(left_bin, right_bin + 1):
            overlap = max(0, min(end, (b + 1) * bin_bp) - max(pos, b * bin_bp))
            covered_bp[b] += overlap
            support_sum[b] += overlap * len(mask) / len(SAMPLE_ORDER)
            if cls == "core":
                core_bp[b] += overlap
            elif cls == "shell":
                shell_bp[b] += overlap
            elif cls == "private":
                private_bp[b] += overlap
        pos = end

    with gzip.open(graph_dir / "chm13_sharing_blocks.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["start", "end", "length", "class", "sample_support", "samples"])
        for start, end, mask, cls in blocks:
            writer.writerow([start, end, end - start, cls, len(mask),
                             ",".join(SAMPLE_ORDER[i] for i in mask)])

    rows = []
    with open(graph_dir / "chm13_support_10kb.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "start", "end", "mean_sample_fraction",
                         "core_fraction", "shell_fraction", "private_fraction"])
        for b in range(n_bins):
            denom = covered_bp[b]
            row = {
                "method": method, "start": b * bin_bp,
                "end": min((b + 1) * bin_bp, chm_length),
                "mean_sample_fraction": support_sum[b] / denom,
                "core_fraction": core_bp[b] / denom,
                "shell_fraction": shell_bp[b] / denom,
                "private_fraction": private_bp[b] / denom,
            }
            rows.append(row)
            writer.writerow([row[k] for k in ["method", "start", "end",
                                               "mean_sample_fraction", "core_fraction",
                                               "shell_fraction", "private_fraction"]])

    with open(graph_dir / "copy_number_distribution.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["max_traversal_count", "nodes", "node_bp"])
        for key in ["0", "1", "2", "3-5", "6-10", ">10"]:
            writer.writerow([key, copy_nodes[key], copy_bp[key]])

    summary = {
        "method": method, "gfa": str(gfa), "graph_bp": sum(lengths.values()),
        "nodes": len(lengths), "visited_nodes": len(visited), "paths": len(sample_paths),
        "chm13_path_bp": chm_length, "sharing_blocks": len(blocks),
    }
    for cls in ["core", "shell", "private", "unvisited"]:
        summary[f"{cls}_nodes"] = class_nodes[cls]
        summary[f"{cls}_bp"] = class_bp[cls]
    return summary, rows


def plot_support(method_rows, output, bin_bp):
    methods = list(method_rows)
    max_bins = max(len(x) for x in method_rows.values())
    matrix = np.full((len(methods), max_bins), np.nan)
    fig, axes = plt.subplots(2, 1, figsize=(16, 5.8), gridspec_kw={"height_ratios": [1, 1.6]})
    for i, method in enumerate(methods):
        values = [x["mean_sample_fraction"] for x in method_rows[method]]
        matrix[i, :len(values)] = values
        x = np.arange(len(values)) * bin_bp / 1e6
        axes[1].plot(x, values, linewidth=1.2, label=method)
    image = axes[0].imshow(matrix, aspect="auto", interpolation="nearest", vmin=0.1, vmax=1,
                           cmap="viridis", extent=[0, max_bins * bin_bp / 1e6,
                                                   len(methods) - 0.5, -0.5])
    axes[0].set_yticks(range(len(methods)), methods)
    axes[0].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    axes[0].set_title("Mean fraction of the 10 sample paths visiting each CHM13-path node")
    fig.colorbar(image, ax=axes[0], label="sample support fraction", pad=0.01)
    axes[1].set_ylim(0, 1.03)
    axes[1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    axes[1].set_ylabel("sample support fraction")
    axes[1].legend(ncol=2, frameon=False)
    axes[1].grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_cross_method_comparison(method_rows, outdir):
    fields = ["method", "mean_chm13_support", "bins", "bins_support_le_0.3",
              "bins_support_ge_0.9", "correlation_with_vallescope2"]
    baseline = np.array([x["mean_sample_fraction"] for x in method_rows["vallescope2"]])
    with open(outdir / "chm13_support_method_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for method, rows in method_rows.items():
            values = np.array([x["mean_sample_fraction"] for x in rows])
            writer.writerow({
                "method": method,
                "mean_chm13_support": values.mean(),
                "bins": len(values),
                "bins_support_le_0.3": int((values <= 0.3).sum()),
                "bins_support_ge_0.9": int((values >= 0.9).sum()),
                "correlation_with_vallescope2": np.corrcoef(baseline, values)[0, 1],
            })

    methods = list(method_rows)
    with open(outdir / "vallescope2_sparse_bins.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["start", "end", *methods])
        for i, base_row in enumerate(method_rows["vallescope2"]):
            if base_row["mean_sample_fraction"] <= 0.3:
                writer.writerow([base_row["start"], base_row["end"],
                                 *(method_rows[m][i]["mean_sample_fraction"] for m in methods)])

    comparators = [m for m in methods if m != "vallescope2"]
    with open(outdir / "vallescope2_sparse_rescued_regions.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["comparator", "start", "end", "bins",
                         "vallescope2_mean_support", "comparator_mean_support"])
        for comparator in comparators:
            active = []
            for i, base_row in enumerate(method_rows["vallescope2"]):
                other = method_rows[comparator][i]
                qualifies = (base_row["mean_sample_fraction"] <= 0.3 and
                             other["mean_sample_fraction"] >= 0.7)
                if qualifies:
                    active.append((base_row, other))
                if active and (not qualifies or i == len(method_rows["vallescope2"]) - 1):
                    writer.writerow([
                        comparator, active[0][0]["start"], active[-1][0]["end"], len(active),
                        np.mean([x[0]["mean_sample_fraction"] for x in active]),
                        np.mean([x[1]["mean_sample_fraction"] for x in active]),
                    ])
                    active = []


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", action="append", required=True,
                        help="METHOD=GFA; repeat for each graph")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bin-bp", type=int, default=10_000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    method_rows = {}
    for item in args.graph:
        method, path = item.split("=", 1)
        summary, rows = analyze_graph(method, Path(path), args.output_dir, args.bin_bp)
        summaries.append(summary)
        method_rows[method] = rows

    fields = list(summaries[0])
    with open(args.output_dir / "graph_sharing_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)

    with open(args.output_dir / "chm13_support_all_methods.tsv", "w") as handle:
        fields = ["method", "start", "end", "mean_sample_fraction",
                  "core_fraction", "shell_fraction", "private_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for rows in method_rows.values():
            writer.writerows(rows)
    plot_support(method_rows, args.output_dir / "chm13_sample_support_comparison.png", args.bin_bp)
    write_cross_method_comparison(method_rows, args.output_dir)


if __name__ == "__main__":
    main()
