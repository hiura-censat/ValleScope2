#!/usr/bin/env python3
import argparse
import csv
import hashlib
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SAMPLES = [
    "chm13v2.0", "HG00146_hap1", "HG00358_hap2", "HG00544_hap2",
    "HG00741_hap1", "HG01114_hap1", "HG01255_hap2", "NA18534_hap2",
    "NA18974_hap2", "NA18982_hap1",
]


def sample_name(name):
    return name.split("|", 1)[0].split("#", 1)[0].split(".chr", 1)[0]


def read_fasta(path):
    chunks = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                chunks.append(line.strip().upper())
    return "".join(chunks)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def parse_gfa(path):
    sequences = {}
    paths = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("S\t"):
                fields = line.rstrip().split("\t")
                sequences[fields[1]] = fields[2]
            elif line.startswith("P\t"):
                fields = line.rstrip().split("\t")
                name = sample_name(fields[1])
                if name in SAMPLES:
                    paths[name] = [(token[:-1], token[-1]) for token in fields[2].split(",")]
    return sequences, paths


def reconstruct_path(steps, sequences):
    chunks = []
    positions = [0]
    for node, orientation in steps:
        sequence = sequences[node]
        chunks.append(sequence if orientation == "+" else reverse_complement(sequence))
        positions.append(positions[-1] + len(sequence))
    return "".join(chunks), positions


def load_rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_fasta(records, path):
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def run_flank_alignment(records, reference_fasta, output_paf, minimap2, threads):
    query_fasta = output_paf.with_suffix(".queries.fa")
    write_fasta(records, query_fasta)
    with open(output_paf, "w") as output:
        subprocess.run([
            minimap2, "-x", "asm20", "-c", "--secondary=yes", "-N", "50",
            "-t", str(threads), str(reference_fasta), str(query_fasta),
        ], check=True, stdout=output)


def parse_flank_alignments(path, expected):
    candidates = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            query = fields[0]
            matches = int(fields[9])
            block = int(fields[10])
            identity = matches / block if block else 0.0
            midpoint = (int(fields[7]) + int(fields[8])) / 2
            candidates[query].append((identity, midpoint, fields[4], int(fields[7]), int(fields[8])))
    result = {}
    for query, target in expected.items():
        hits = candidates.get(query, [])
        if not hits:
            result[query] = {"identity": 0.0, "distance": math.inf, "strand": ".", "hits": 0}
            continue
        best = min(hits, key=lambda hit: (abs(hit[1] - target), -hit[0]))
        result[query] = {
            "identity": best[0], "distance": abs(best[1] - target),
            "strand": best[2], "hits": len(hits),
        }
    return result


def interval_overlap(start, end, intervals):
    total = max(1, end - start)
    overlap = sum(max(0, min(end, right) - max(start, left)) for left, right, _ in intervals)
    labels = Counter()
    for left, right, label in intervals:
        amount = max(0, min(end, right) - max(start, left))
        labels[label] += amount
    return min(1.0, overlap / total), labels


def label_similarity(left, right):
    keys = left.keys() | right.keys()
    numerator = sum(min(left[key], right[key]) for key in keys)
    denominator = sum(max(left[key], right[key]) for key in keys)
    return numerator / denominator if denominator else 1.0


def confidence_score(ref_span, sample_span, left, right):
    gap_score = math.exp(-ref_span / 100000.0)
    ratio = sample_span / max(1.0, ref_span)
    stretch_score = math.exp(-abs(math.log(max(ratio, 1e-6))))
    flank_identity = (left["identity"] + right["identity"]) / 2
    distances = [value["distance"] for value in (left, right) if math.isfinite(value["distance"])]
    proximity = math.exp(-(sum(distances) / len(distances)) / 25000.0) if distances else 0.0
    flank_score = flank_identity * proximity
    score = 0.35 * gap_score + 0.20 * stretch_score + 0.45 * flank_score
    return score, gap_score, stretch_score, flank_score


def confidence_class(score):
    if score >= 0.75:
        return "high"
    if score >= 0.45:
        return "moderate"
    return "low"


def shannon_entropy(groups):
    total = sum(groups)
    return -sum((size / total) * math.log2(size / total) for size in groups if size)


def validate_method(method, gfa, phase_b, fastas, annotations, output, minimap2, threads,
                    high_controls):
    method_dir = output / method
    method_dir.mkdir(parents=True, exist_ok=True)
    sequences, paths = parse_gfa(gfa)
    path_sequences = {}
    path_positions = {}
    path_exact = {}
    for sample in SAMPLES:
        sequence, positions = reconstruct_path(paths[sample], sequences)
        path_sequences[sample] = sequence
        path_positions[sample] = positions
        path_exact[sample] = sequence == fastas[sample]

    branches = load_rows(phase_b / method / "nonreference_branches.tsv")
    selected = []
    controls_seen = Counter()
    for row in branches:
        if row["confidence"] == "low" or controls_seen[row["sample"]] < high_controls:
            selected.append(row)
            if row["confidence"] != "low":
                controls_seen[row["sample"]] += 1

    flank_records = []
    expected = {}
    flank_bp = 1000
    for row in selected:
        sample = row["sample"]
        start, end = int(row["sample_start"]), int(row["sample_end"])
        left_name, right_name = row["branch_id"] + "|L", row["branch_id"] + "|R"
        left_sequence = fastas[sample][max(0, start - flank_bp):start]
        right_sequence = fastas[sample][end:min(len(fastas[sample]), end + flank_bp)]
        if left_sequence:
            flank_records.append((left_name, left_sequence))
            expected[left_name] = float(row["projected_start"])
        if right_sequence:
            flank_records.append((right_name, right_sequence))
            expected[right_name] = float(row["projected_end"])
    paf = method_dir / "validation_flanks_to_chm13.paf"
    if not paf.exists():
        run_flank_alignment(flank_records, output / "chm13.reference.fa", paf, minimap2, threads)
    flank_hits = parse_flank_alignments(paf, expected)

    fields = [
        "method", "branch_id", "sample", "sample_start", "sample_end", "branch_bp",
        "projected_start", "projected_end", "anchor_ref_span", "anchor_sample_span",
        "path_matches_original_fasta", "branch_matches_original_fasta", "left_identity",
        "left_projection_error_bp", "right_identity", "right_projection_error_bp",
        "sample_alpha_fraction", "projected_alpha_fraction", "alpha_label_similarity",
        "gap_score", "stretch_score", "flank_score", "projection_confidence_score",
        "projection_confidence_class", "phase_b_confidence",
    ]
    validated = []
    selected_ids = {row["branch_id"] for row in selected}
    for row in branches:
        sample = row["sample"]
        start, end = int(row["sample_start"]), int(row["sample_end"])
        projected_start, projected_end = int(float(row["projected_start"])), int(float(row["projected_end"]))
        branch_sequence = path_sequences[sample][start:end]
        branch_exact = branch_sequence == fastas[sample][start:end]
        ref_span = int(row["anchor_ref_end"]) - int(row["anchor_ref_start"])
        sample_span = end - start
        left = flank_hits.get(row["branch_id"] + "|L", {"identity": 0.0, "distance": math.inf})
        right = flank_hits.get(row["branch_id"] + "|R", {"identity": 0.0, "distance": math.inf})
        if row["branch_id"] in selected_ids:
            score, gap_score, stretch_score, flank_score = confidence_score(ref_span, sample_span, left, right)
        else:
            gap_score = math.exp(-ref_span / 100000.0)
            stretch_score = math.exp(-abs(math.log(max(sample_span / max(1.0, ref_span), 1e-6))))
            flank_score = float("nan")
            score = 0.35 * gap_score + 0.20 * stretch_score + 0.45
        sample_alpha, sample_labels = interval_overlap(start, end, annotations[sample])
        ref_alpha, ref_labels = interval_overlap(projected_start, projected_end, annotations["chm13v2.0"])
        validated.append({
            "method": method, "branch_id": row["branch_id"], "sample": sample,
            "sample_start": start, "sample_end": end, "branch_bp": end - start,
            "projected_start": projected_start, "projected_end": projected_end,
            "anchor_ref_span": ref_span, "anchor_sample_span": sample_span,
            "path_matches_original_fasta": int(path_exact[sample]),
            "branch_matches_original_fasta": int(branch_exact),
            "left_identity": left.get("identity", float("nan")),
            "left_projection_error_bp": left.get("distance", float("nan")),
            "right_identity": right.get("identity", float("nan")),
            "right_projection_error_bp": right.get("distance", float("nan")),
            "sample_alpha_fraction": sample_alpha, "projected_alpha_fraction": ref_alpha,
            "alpha_label_similarity": label_similarity(sample_labels, ref_labels),
            "gap_score": gap_score, "stretch_score": stretch_score, "flank_score": flank_score,
            "projection_confidence_score": score,
            "projection_confidence_class": confidence_class(score),
            "phase_b_confidence": row["confidence"],
        })
    with open(method_dir / "branch_projection_validation.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(validated)
    return validated, path_exact


def analyze_blocks(method, phase_b, method_dir):
    rows = load_rows(phase_b / method / "full_sharing_blocks.tsv")
    by_block = defaultdict(list)
    signatures = defaultdict(list)
    for row in rows:
        by_block[row["block_id"]].append(row)
        signatures[row["signature"]].append(row)
    block_fields = ["method", "block_id", "start", "end", "length", "number_of_groups",
                    "group_entropy", "largest_group_fraction"]
    block_stats = []
    for block_id, group_rows in by_block.items():
        sizes = [len(row["haplotypes"].split(",")) for row in group_rows]
        first = group_rows[0]
        block_stats.append({
            "method": method, "block_id": block_id, "start": first["start"], "end": first["end"],
            "length": first["length"], "number_of_groups": len(sizes),
            "group_entropy": shannon_entropy(sizes), "largest_group_fraction": max(sizes) / len(SAMPLES),
        })
    with open(method_dir / "block_statistics.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=block_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(block_stats)
    recurrence = []
    for signature, group_rows in signatures.items():
        intervals = sorted({(int(row["start"]), int(row["end"])) for row in group_rows})
        recurrence.append({
            "method": method, "signature": signature,
            "haplotypes": group_rows[0]["haplotypes"], "occurrences": len(intervals),
            "total_bp": sum(end - start for start, end in intervals),
            "intervals": ",".join(f"{start}-{end}" for start, end in intervals),
        })
    recurrence.sort(key=lambda row: (-row["total_bp"], -row["occurrences"], row["signature"]))
    with open(method_dir / "group_signature_recurrence.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(recurrence[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(recurrence)
    return block_stats, recurrence


def plot_confidence(all_validated, output, path_length):
    methods = list(all_validated)
    colors = {"high": "#2ca25f", "moderate": "#fec44f", "low": "#de2d26"}
    fig, axes = plt.subplots(len(methods), 1, figsize=(16, 1.8 * len(methods)), sharex=True)
    for ax, method in zip(axes, methods):
        for row in all_validated[method]:
            start, end = row["projected_start"] / 1e6, row["projected_end"] / 1e6
            confidence = row["projection_confidence_class"]
            color = colors[confidence] if not math.isnan(row["flank_score"]) else "#bdbdbd"
            ax.plot([start, max(start + 0.001, end)], [0.5, 0.5], lw=4,
                    color=color, alpha=0.75)
        ax.set_xlim(0, path_length / 1e6); ax.set_yticks([]); ax.set_ylabel(method, rotation=0, ha="right")
    axes[-1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    fig.suptitle("Phase C projection confidence: high / moderate / low; gray = not flank-tested")
    fig.tight_layout(); fig.savefig(output, dpi=180); plt.close(fig)


def plot_recurrence(all_recurrence, output, path_length, top_n=12):
    methods = list(all_recurrence)
    fig, axes = plt.subplots(len(methods), 1, figsize=(16, 3.0 * len(methods)), sharex=True)
    for ax, method in zip(axes, methods):
        rows = [row for row in all_recurrence[method]
                if 1 < row["occurrences"] and 1 < len(row["haplotypes"].split(",")) < len(SAMPLES)][:top_n]
        for y, row in enumerate(rows):
            for interval in row["intervals"].split(","):
                start, end = map(int, interval.split("-"))
                ax.broken_barh([(start / 1e6, (end - start) / 1e6)], (y - 0.35, 0.7))
        ax.set_yticks(range(len(rows)))
        labels = []
        for index, row in enumerate(rows, 1):
            member_count = len(row["haplotypes"].split(","))
            labels.append(f"G{index:02d}: n={member_count}, occurrences={row['occurrences']}")
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_ylabel(method, rotation=0, ha="right")
        ax.set_xlim(0, path_length / 1e6)
    axes[-1].set_xlabel("CHM13 alphaSat path coordinate (Mb)")
    fig.suptitle("Recurrent haplotype-group signatures across the full array")
    fig.tight_layout(); fig.savefig(output, dpi=180); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", action="append", required=True, help="METHOD=GFA")
    parser.add_argument("--phase-b", type=Path, required=True)
    parser.add_argument("--fasta-list", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--high-controls-per-sample", type=int, default=20)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fasta_paths = [Path(line.strip()) for line in open(args.fasta_list) if line.strip()]
    fastas = {sample_name(path.name): read_fasta(path) for path in fasta_paths}
    write_fasta([("chm13v2.0", fastas["chm13v2.0"])], args.output_dir / "chm13.reference.fa")
    annotations = {}
    for sample in SAMPLES:
        annotation = args.annotation_dir / f"{sample}.alpha.bed"
        intervals = []
        with open(annotation) as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                intervals.append((int(fields[1]), int(fields[2]), fields[3] if len(fields) > 3 else "alpha"))
        annotations[sample] = intervals
    all_validated, all_blocks, all_recurrence = {}, {}, {}
    path_exact_rows = []
    for item in args.graph:
        method, graph = item.split("=", 1)
        validated, path_exact = validate_method(
            method, Path(graph), args.phase_b, fastas, annotations, args.output_dir,
            args.minimap2, args.threads, args.high_controls_per_sample)
        all_validated[method] = validated
        for sample, exact in path_exact.items():
            path_exact_rows.append({"method": method, "sample": sample, "path_matches_fasta": int(exact)})
        blocks, recurrence = analyze_blocks(method, args.phase_b, args.output_dir / method)
        all_blocks[method], all_recurrence[method] = blocks, recurrence
    with open(args.output_dir / "path_fasta_consistency.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(path_exact_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(path_exact_rows)
    summary = []
    for method, rows in all_validated.items():
        tested = [row for row in rows if not math.isnan(row["flank_score"])]
        counts = Counter(row["projection_confidence_class"] for row in tested)
        identities = [value for row in tested for value in (row["left_identity"], row["right_identity"])
                      if not math.isnan(value)]
        errors = [value for row in tested for value in
                  (row["left_projection_error_bp"], row["right_projection_error_bp"])
                  if math.isfinite(value)]
        summary.append({
            "method": method, "branches": len(rows), "flank_validated_branches": len(tested),
            "fasta_exact_branches": sum(row["branch_matches_original_fasta"] for row in rows),
            "high": counts["high"], "moderate": counts["moderate"], "low": counts["low"],
            "mean_confidence": np.mean([row["projection_confidence_score"] for row in tested]),
            "mean_flank_identity": np.mean(identities),
            "median_flank_projection_error_bp": np.median(errors),
            "mean_alpha_label_similarity": np.mean([row["alpha_label_similarity"] for row in tested]),
        })
    with open(args.output_dir / "projection_validation_summary.tsv", "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary)
    path_length = len(fastas["chm13v2.0"])
    plot_confidence(all_validated, args.output_dir / "projection_confidence_full_length.png", path_length)
    plot_recurrence(all_recurrence, args.output_dir / "recurrent_group_sharing_map.png", path_length)


if __name__ == "__main__":
    main()
