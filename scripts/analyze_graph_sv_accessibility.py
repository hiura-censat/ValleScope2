#!/usr/bin/env python3
"""Analyze reference accessibility and graph support of Pantree SVs."""

import argparse
import csv
import gzip
import math
import re
import statistics
from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_pantree_vg_sv import maximum_match, read_events


Edge = Tuple[str, str]


def edge_key(left: str, right: str) -> Edge:
    return tuple(sorted((left, right)))


def parse_path_nodes(value: str) -> List[str]:
    return [token[:-1] for token in value.split(",") if token]


class UnionFind:
    def __init__(self) -> None:
        self.parent: Dict[str, str] = {}
        self.rank: Dict[str, int] = {}

    def find(self, item: str) -> str:
        if item not in self.parent:
            self.parent[item] = item
            self.rank[item] = 0
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def read_gfa(path: Path):
    node_lengths: Dict[str, int] = {}
    graph_edges: Set[Edge] = set()
    paths: Dict[str, List[str]] = {}
    union = UnionFind()
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields[0] == "S":
                node_lengths[fields[1]] = len(fields[2])
                union.find(fields[1])
            elif fields[0] == "L":
                edge = edge_key(fields[1], fields[3])
                graph_edges.add(edge)
                union.union(*edge)
            elif fields[0] == "P":
                paths[fields[1]] = parse_path_nodes(fields[2])

    node_paths: Dict[str, Set[str]] = defaultdict(set)
    node_positions: Dict[str, Dict[str, List[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    edge_paths: Dict[Edge, Set[str]] = defaultdict(set)
    for path_name, nodes in paths.items():
        position = 0
        for index, node in enumerate(nodes):
            node_paths[node].add(path_name)
            node_positions[node][path_name].append(position)
            position += node_lengths[node]
            if index:
                edge_paths[edge_key(nodes[index - 1], node)].add(path_name)

    adjacency: Dict[str, Set[str]] = defaultdict(set)
    for left, right in graph_edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    return (
        node_lengths,
        graph_edges,
        paths,
        node_paths,
        node_positions,
        edge_paths,
        adjacency,
        union,
    )


def reference_distances(adjacency, reference_nodes: Set[str]) -> Dict[str, int]:
    distances = {node: 0 for node in reference_nodes}
    queue = deque(reference_nodes)
    while queue:
        node = queue.popleft()
        for neighbor in adjacency.get(node, ()):
            if neighbor not in distances:
                distances[neighbor] = distances[node] + 1
                queue.append(neighbor)
    return distances


def parse_info(value: str) -> Dict[str, str]:
    return {
        key: item_value
        for item in value.split(";")
        if "=" in item
        for key, item_value in [item.split("=", 1)]
    }


def allele_length(allele: str) -> int:
    return 0 if allele == "." else len(allele)


def read_pantree_vcf(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    records = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            info = parse_info(fields[7])
            nodes = re.findall(r"[<>]([^<>]+)", fields[2])
            records.append(
                {
                    "pos": int(fields[1]),
                    "id": fields[2],
                    "nodes": nodes[:2],
                    "kind": info.get("VT", ""),
                    "size": abs(
                        allele_length(fields[4]) - allele_length(fields[3])
                    ),
                    "ac": int(info.get("AC", "0").split(",")[0]),
                    "an": int(info.get("AN", "0")),
                }
            )
    return records


def read_paf(path: Optional[Path]):
    by_pair = defaultdict(list)
    if path is None:
        return by_pair
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            first, second = fields[0], fields[5]
            intervals = {
                first: (int(fields[2]), int(fields[3])),
                second: (int(fields[7]), int(fields[8])),
            }
            by_pair[tuple(sorted((first, second)))].append(intervals)
    return by_pair


def locus_pair_support(
    nodes: List[str],
    node_positions: Dict[str, Dict[str, List[int]]],
    paf_by_pair,
) -> int:
    if not paf_by_pair:
        return -1
    path_positions: Dict[str, List[int]] = defaultdict(list)
    for node in nodes:
        for path, positions in node_positions.get(node, {}).items():
            path_positions[path].extend(positions)
    supported = 0
    for first, second in combinations(sorted(path_positions), 2):
        alignments = paf_by_pair.get(tuple(sorted((first, second))), ())
        if any(
            start_a <= pos_a <= end_a and start_b <= pos_b <= end_b
            for intervals in alignments
            for start_a, end_a in [intervals[first]]
            for start_b, end_b in [intervals[second]]
            for pos_a in path_positions[first]
            for pos_b in path_positions[second]
        ):
            supported += 1
    return supported


def an_bucket(an: int) -> str:
    if an <= 3:
        return "AN=2-3"
    if an <= 6:
        return "AN=4-6"
    if an <= 9:
        return "AN=7-9"
    return "AN=10"


def embedding_rows(label, node_lengths, paths, node_paths, adjacency, reference):
    reference_nodes = set(paths[reference])
    distances = reference_distances(adjacency, reference_nodes)
    rows = []
    for path, nodes in paths.items():
        counts = Counter()
        hop_bp = 0
        for node in nodes:
            length = node_lengths[node]
            if node in reference_nodes:
                counts["direct"] += length
            elif node in distances:
                counts["mediated"] += length
            else:
                counts["isolated"] += length
            if node in distances:
                hop_bp += distances[node] * length
            if node not in reference_nodes and len(node_paths[node]) >= 2:
                counts["other_shared"] += length
            elif node not in reference_nodes:
                counts["sample_unique"] += length
        total = sum(counts[key] for key in ("direct", "mediated", "isolated"))
        rows.append(
            {
                "method": label,
                "sample": path,
                "path_bp": total,
                "direct_bp": counts["direct"],
                "mediated_bp": counts["mediated"],
                "isolated_bp": counts["isolated"],
                "direct_rate": counts["direct"] / total,
                "mediated_rate": counts["mediated"] / total,
                "isolated_rate": counts["isolated"] / total,
                "mean_reference_hops": hop_bp
                / max(1, counts["direct"] + counts["mediated"]),
                "other_sample_shared_bp": counts["other_shared"],
                "sample_unique_bp": counts["sample_unique"],
            }
        )
    return rows


def merged_length(intervals: List[Tuple[int, int]]) -> int:
    total = 0
    end = -1
    for start, next_end in sorted(intervals):
        if next_end <= end:
            continue
        total += next_end - max(start, end)
        end = next_end
    return total


def paf_embedding_rows(label: str, paf: Optional[Path], reference: str):
    if paf is None:
        return []
    lengths = {}
    all_intervals = defaultdict(list)
    direct_intervals = defaultdict(list)
    with paf.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            first, second = fields[0], fields[5]
            lengths[first], lengths[second] = int(fields[1]), int(fields[6])
            all_intervals[first].append((int(fields[2]), int(fields[3])))
            all_intervals[second].append((int(fields[7]), int(fields[8])))
            if first == reference:
                direct_intervals[second].append((int(fields[7]), int(fields[8])))
            elif second == reference:
                direct_intervals[first].append((int(fields[2]), int(fields[3])))
    rows = []
    for sample, length in sorted(lengths.items()):
        if sample == reference:
            continue
        any_bp = merged_length(all_intervals[sample])
        direct_bp = merged_length(direct_intervals[sample])
        mediated_bp = max(0, any_bp - direct_bp)
        isolated_bp = max(0, length - any_bp)
        rows.append(
            {
                "method": label,
                "sample": sample,
                "sequence_bp": length,
                "direct_CHM13_alignment_bp": direct_bp,
                "other_sample_only_alignment_bp": mediated_bp,
                "unaligned_bp": isolated_bp,
                "direct_rate": direct_bp / length,
                "mediated_rate": mediated_bp / length,
                "unaligned_rate": isolated_bp / length,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--gfa", action="append", type=Path, required=True)
    parser.add_argument(
        "--embedding-gfa",
        action="append",
        type=Path,
        help="Optional graph used only for graph-embedding rates.",
    )
    parser.add_argument("--pantree-vcf", action="append", type=Path, required=True)
    parser.add_argument("--vg-vcf", action="append", type=Path, required=True)
    parser.add_argument("--reference", action="append", required=True)
    parser.add_argument("--paf", action="append", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    count = len(args.label)
    if not all(
        len(values) == count
        for values in (args.gfa, args.pantree_vcf, args.vg_vcf, args.reference)
    ):
        parser.error("all repeated arguments must have equal counts")
    paf_paths = args.paf or []
    paf_paths += [None] * (count - len(paf_paths))
    embedding_gfas = args.embedding_gfa or []
    embedding_gfas += [None] * (count - len(embedding_gfas))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_events = []
    all_embedding = []
    all_paf_embedding = []
    all_an = []
    for label, gfa, embedding_gfa, pantree_vcf, vg_vcf, reference, paf in zip(
        args.label,
        args.gfa,
        embedding_gfas,
        args.pantree_vcf,
        args.vg_vcf,
        args.reference,
        paf_paths,
    ):
        (
            node_lengths,
            graph_edges,
            paths,
            node_paths,
            node_positions,
            edge_paths,
            adjacency,
            union,
        ) = read_gfa(gfa)
        if reference not in paths:
            raise RuntimeError("{}: reference path is absent: {}".format(label, reference))
        reference_nodes = set(paths[reference])
        distances = reference_distances(adjacency, reference_nodes)
        reference_components = {union.find(node) for node in reference_nodes}
        paf_by_pair = read_paf(paf)
        all_paf_embedding.extend(paf_embedding_rows(label, paf, reference))

        records = read_pantree_vcf(pantree_vcf)
        pantree_events = read_events(pantree_vcf)
        vg_events = read_events(vg_vcf)
        matches = maximum_match(pantree_events, vg_events, 500, 0.7)
        for index, record in enumerate(records):
            nodes = record["nodes"]
            edge = edge_key(*nodes) if len(nodes) == 2 else ("", "")
            endpoint_on_reference = sum(node in reference_nodes for node in nodes)
            components = {union.find(node) for node in nodes if node in union.parent}
            if endpoint_on_reference:
                accessibility = "direct"
            elif components & reference_components:
                accessibility = "mediated"
            else:
                accessibility = "isolated"
            min_hops = min((distances.get(node, -1) for node in nodes), default=-1)
            edge_support = len(edge_paths.get(edge, ()))
            locus_paths = set()
            for node in nodes:
                locus_paths.update(node_paths.get(node, ()))
            support_entropy = (
                math.log(len(locus_paths)) / math.log(len(paths))
                if len(locus_paths) > 1 and len(paths) > 1
                else 0.0
            )
            all_events.append(
                {
                    "method": label,
                    **record,
                    "vg_matched_500bp": int(matches[index] >= 0),
                    "accessibility": accessibility,
                    "endpoint_on_reference_count": endpoint_on_reference,
                    "min_reference_hops": min_hops,
                    "variant_edge_in_gfa": int(edge in graph_edges),
                    "variant_edge_path_support": edge_support,
                    "locus_path_support": len(locus_paths),
                    "path_support_entropy": support_entropy,
                    "locus_spanning_paf_pairs": locus_pair_support(
                        nodes, node_positions, paf_by_pair
                    ),
                }
            )

        if embedding_gfa is None or embedding_gfa == gfa:
            embedding_data = (node_lengths, paths, node_paths, adjacency)
        else:
            (
                embedding_lengths,
                _,
                embedding_paths,
                embedding_node_paths,
                _,
                _,
                embedding_adjacency,
                _,
            ) = read_gfa(embedding_gfa)
            embedding_data = (
                embedding_lengths,
                embedding_paths,
                embedding_node_paths,
                embedding_adjacency,
            )
        all_embedding.extend(
            embedding_rows(label, *embedding_data, reference)
        )
        for bucket in ("AN=2-3", "AN=4-6", "AN=7-9", "AN=10"):
            indices = [
                index
                for index, record in enumerate(records)
                if an_bucket(record["an"]) == bucket
            ]
            matched = sum(matches[index] >= 0 for index in indices)
            all_an.append(
                {
                    "method": label,
                    "AN_bucket": bucket,
                    "pantree_SVs": len(indices),
                    "vg_matched": matched,
                    "match_rate": matched / len(indices) if indices else 0.0,
                }
            )

    event_fields = list(all_events[0])
    with (args.output_dir / "sv_graph_metrics.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=event_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(all_events)
    with (args.output_dir / "graph_embedding.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(all_embedding[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(all_embedding)
    with (args.output_dir / "an_stratified_match.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(all_an[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(all_an)
    if all_paf_embedding:
        with (args.output_dir / "paf_embedding.tsv").open("w") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(all_paf_embedding[0]),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(all_paf_embedding)

    summary_rows = []
    for label in args.label:
        items = [row for row in all_events if row["method"] == label]
        for matched_label, matched_value in (
            ("Pantree-only", 0),
            ("Pantree-vg matched", 1),
        ):
            subset = [
                row for row in items if row["vg_matched_500bp"] == matched_value
            ]
            for key in ("direct", "mediated", "isolated"):
                count_value = sum(row["accessibility"] == key for row in subset)
                summary_rows.append(
                    {
                        "method": label,
                        "caller_class": matched_label,
                        "metric": "accessibility_" + key,
                        "count": count_value,
                        "fraction": count_value / len(subset) if subset else 0.0,
                    }
                )
            for threshold, name in ((1, "edge_support_1"), (3, "edge_support_2_3")):
                if threshold == 1:
                    count_value = sum(
                        row["variant_edge_path_support"] <= 1 for row in subset
                    )
                else:
                    count_value = sum(
                        2 <= row["variant_edge_path_support"] <= 3 for row in subset
                    )
                summary_rows.append(
                    {
                        "method": label,
                        "caller_class": matched_label,
                        "metric": name,
                        "count": count_value,
                        "fraction": count_value / len(subset) if subset else 0.0,
                    }
                )
            count_value = sum(
                row["variant_edge_path_support"] >= 4 for row in subset
            )
            summary_rows.append(
                {
                    "method": label,
                    "caller_class": matched_label,
                    "metric": "edge_support_ge4",
                    "count": count_value,
                    "fraction": count_value / len(subset) if subset else 0.0,
                }
            )
    with (args.output_dir / "sv_graph_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    support_rows = []
    for label in args.label:
        for class_label, matched_value in (
            ("Pantree-only", 0),
            ("Pantree-vg matched", 1),
        ):
            subset = [
                row
                for row in all_events
                if row["method"] == label
                and row["vg_matched_500bp"] == matched_value
            ]
            locus_support = [row["locus_path_support"] for row in subset]
            paf_support = [
                row["locus_spanning_paf_pairs"]
                for row in subset
                if row["locus_spanning_paf_pairs"] >= 0
            ]
            support_rows.append(
                {
                    "method": label,
                    "caller_class": class_label,
                    "SVs": len(subset),
                    "edge_path_support_1_fraction": sum(
                        row["variant_edge_path_support"] <= 1 for row in subset
                    )
                    / len(subset),
                    "edge_path_support_2_3_fraction": sum(
                        2 <= row["variant_edge_path_support"] <= 3 for row in subset
                    )
                    / len(subset),
                    "edge_path_support_ge4_fraction": sum(
                        row["variant_edge_path_support"] >= 4 for row in subset
                    )
                    / len(subset),
                    "locus_path_support_median": statistics.median(locus_support),
                    "locus_path_support_ge4_fraction": sum(
                        value >= 4 for value in locus_support
                    )
                    / len(locus_support),
                    "mean_path_support_entropy": sum(
                        row["path_support_entropy"] for row in subset
                    )
                    / len(subset),
                    "locus_spanning_paf_pairs_median": statistics.median(paf_support)
                    if paf_support
                    else -1,
                    "locus_spanning_paf_pairs_0_1_fraction": sum(
                        value <= 1 for value in paf_support
                    )
                    / len(paf_support)
                    if paf_support
                    else -1,
                    "locus_spanning_paf_pairs_2_3_fraction": sum(
                        2 <= value <= 3 for value in paf_support
                    )
                    / len(paf_support)
                    if paf_support
                    else -1,
                    "locus_spanning_paf_pairs_ge4_fraction": sum(
                        value >= 4 for value in paf_support
                    )
                    / len(paf_support)
                    if paf_support
                    else -1,
                }
            )
    with (args.output_dir / "support_class_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(support_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(support_rows)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=180)
    buckets = ("AN=2-3", "AN=4-6", "AN=7-9", "AN=10")
    colors = ("#2F6B9A", "#D97732")
    x = list(range(len(buckets)))
    for method_index, label in enumerate(args.label):
        values = [
            next(
                row["match_rate"]
                for row in all_an
                if row["method"] == label and row["AN_bucket"] == bucket
            )
            for bucket in buckets
        ]
        axes[0].bar(
            [position + (method_index - 0.5) * 0.36 for position in x],
            values,
            0.36,
            label=label,
            color=colors[method_index],
        )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(buckets)
    axes[0].set_ylabel("Pantree SV fraction matched by vg")
    axes[0].set_ylim(0, 1)
    axes[0].legend(frameon=False)

    non_reference_embedding = [
        row for row in all_embedding if row["sample"] not in args.reference
    ]
    labels = args.label
    direct = [
        sum(row["direct_bp"] for row in non_reference_embedding if row["method"] == label)
        / sum(row["path_bp"] for row in non_reference_embedding if row["method"] == label)
        for label in labels
    ]
    mediated = [
        sum(
            row["mediated_bp"]
            for row in non_reference_embedding
            if row["method"] == label
        )
        / sum(row["path_bp"] for row in non_reference_embedding if row["method"] == label)
        for label in labels
    ]
    isolated = [1 - direct[index] - mediated[index] for index in range(len(labels))]
    axes[1].bar(labels, direct, color="#2A9D8F", label="On CHM13 nodes")
    axes[1].bar(
        labels,
        mediated,
        bottom=direct,
        color="#E9A23B",
        label="Same CHM13 component",
    )
    axes[1].bar(
        labels,
        isolated,
        bottom=[direct[index] + mediated[index] for index in range(len(labels))],
        color="#8C8C8C",
        label="Different component",
    )
    axes[1].set_ylabel("Fraction of non-reference path bp")
    axes[1].set_ylim(0, 1)
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.grid(axis="y", alpha=0.25)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output_dir / "an_match_and_embedding.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
