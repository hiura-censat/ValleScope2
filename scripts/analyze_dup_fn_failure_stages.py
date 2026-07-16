#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([=XIDM])")


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def load_anchors(path: Path) -> dict[str, tuple[str, int, int]]:
    anchors = {}
    for row in read_tsv(path):
        anchors[row["anchor_id"]] = (
            row["sequence_id"], int(row["start"]), int(row["end"])
        )
    return anchors


def load_pair_anchors(path: Path, anchors: dict[str, tuple[str, int, int]]):
    rows = []
    for row in read_tsv(path):
        left = anchors.get(row["anchor_a"])
        right = anchors.get(row["anchor_b"])
        if left and right:
            rows.append((row["sample_a"], row["sample_b"], left, right, row))
    return rows


def oriented_anchor_pair(item, sample: str):
    sample_a, sample_b, left, right, _ = item
    if sample_a == "reference" and sample_b == sample:
        return left, right
    if sample_b == "reference" and sample_a == sample:
        return right, left
    return None


def anchor_counts(items, sample: str, ref_interval, copies):
    counts = [0, 0]
    ref0, ref1 = ref_interval
    for item in items:
        pair = oriented_anchor_pair(item, sample)
        if not pair:
            continue
        ref_anchor, query_anchor = pair
        if not overlap(ref_anchor[1], ref_anchor[2], ref0, ref1):
            continue
        for index, (q0, q1) in enumerate(copies):
            if overlap(query_anchor[1], query_anchor[2], q0, q1):
                counts[index] += 1
    return counts


def read_fn_dup_records(benchmark: Path, samples: list[str]):
    result = {}
    for sample in samples:
        path = benchmark / "samples" / sample / "truvari" / "fn.vcf.gz"
        with open_text(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                info = parse_info(fields[7])
                if info.get("SVTYPE") == "DUP":
                    result[fields[2]] = info
    return result


def parse_info(text: str) -> dict[str, str]:
    info = {}
    for field in text.split(";"):
        key, sep, value = field.partition("=")
        info[key] = value if sep else "true"
    return info


def read_paf(path: Path):
    rows = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for value in fields[12:]:
                parts = value.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            rows.append({
                "qname": fields[0], "qstart": int(fields[2]),
                "qend": int(fields[3]), "strand": fields[4],
                "tname": fields[5], "tstart": int(fields[7]),
                "tend": int(fields[8]), "cigar": tags.get("cg", ""),
                "chain_id": tags.get("ch", "."),
                "bundle_type": tags.get("bt", "."),
            })
    return rows


def sample_reference_paf(row, sample: str):
    if row["qname"] == sample and row["tname"] == "reference":
        return row
    return None


def cigar_insertions(row) -> list[tuple[int, int, int]]:
    query = row["qstart"]
    target = row["tstart"]
    insertions = []
    for length_text, operation in CIGAR_RE.findall(row["cigar"]):
        length = int(length_text)
        if operation == "I":
            insertions.append((target, query, length))
        if operation in "=XM I".replace(" ", ""):
            query += length
        if operation in "=XMD":
            target += length
    return insertions


def paf_evidence(paf, sample: str, ref_interval, copies, dup_length: int):
    ref0, ref1 = ref_interval
    nearby = []
    copy_support = [0, 0]
    insertion_lengths = []
    for raw in paf:
        row = sample_reference_paf(raw, sample)
        if not row or not overlap(row["tstart"], row["tend"], ref0, ref1):
            continue
        nearby.append(row)
        ref_overlap = overlap(row["tstart"], row["tend"], ref0, ref1)
        for index, (q0, q1) in enumerate(copies):
            query_overlap = overlap(row["qstart"], row["qend"], q0, q1)
            if ref_overlap >= max(1, dup_length // 2) and query_overlap >= max(1, dup_length // 2):
                copy_support[index] += 1
        for target, query, length in cigar_insertions(row):
            if ref0 - 500 <= target <= ref1 + 500 and (
                overlap(query, query + length, copies[0][0], copies[1][1]) or
                copies[0][0] - 500 <= query <= copies[1][1] + 500
            ):
                insertion_lengths.append(length)
    size_match = any(0.7 * dup_length <= length <= 1.3 * dup_length
                     for length in insertion_lengths)
    duplicate_mapping = all(value > 0 for value in copy_support)
    return nearby, copy_support, insertion_lengths, size_match or duplicate_mapping


def allele_size(ref: str, alt: str, info: dict[str, str]) -> int:
    if "SVLEN" in info:
        try:
            return abs(int(info["SVLEN"].split(",")[0]))
        except ValueError:
            pass
    if not alt.startswith("<"):
        return abs(len(alt.split(",")[0]) - len(ref))
    try:
        return abs(int(info.get("END", "0")))
    except ValueError:
        return 0


def load_graph_calls(benchmark: Path, sample: str):
    path = benchmark / "samples" / sample / "graph.sample.vcf.gz"
    calls = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue
            genotype = fields[9].split(":", 1)[0]
            if genotype in {".", "./.", "0", "0/0", "0|0"}:
                continue
            info = parse_info(fields[7])
            calls.append({
                "position": int(fields[1]) - 1,
                "size": allele_size(fields[3], fields[4], info),
                "id": fields[2],
                "ref": fields[3],
                "alt": fields[4].split(",", 1)[0],
            })
    return calls


def read_fasta_sequence(path: Path) -> str:
    sequence = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                sequence.append(line.strip().upper())
    return "".join(sequence)


def inserted_sequence(call) -> str:
    ref = call["ref"].upper()
    alt = call["alt"].upper()
    if alt.startswith("<"):
        return ""
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < len(ref) - prefix and suffix < len(alt) - prefix and
           ref[-suffix - 1] == alt[-suffix - 1]):
        suffix += 1
    return alt[prefix:len(alt) - suffix if suffix else None]


def shifted_identity(query: str, target: str, max_shift: int = 20) -> float:
    """Ungapped identity with a small breakpoint-shift allowance."""
    if not query or not target:
        return 0.0
    denominator = max(len(query), len(target))
    best = 0
    for shift in range(-max_shift, max_shift + 1):
        query_start = max(0, shift)
        target_start = max(0, -shift)
        aligned = min(len(query) - query_start, len(target) - target_start)
        if aligned <= 0:
            continue
        matches = sum(
            query[query_start + offset] == target[target_start + offset]
            for offset in range(aligned)
        )
        best = max(best, matches)
    return best / denominator


def best_graph_dup_identity(calls, source_sequence: str):
    scored = [
        (shifted_identity(inserted_sequence(call), source_sequence), call)
        for call in calls
    ]
    return max(scored, default=(0.0, None), key=lambda item: item[0])


def parse_metric(info, name):
    try:
        return float(info[name])
    except (KeyError, ValueError):
        return None


def classify(row, correspondence_counts, chain_counts, paf_result, graph_calls,
             fn_info, direct_identity):
    nearby_paf, copy_paf, insertion_lengths, paf_dup_signal = paf_result
    sequence_similarity = parse_metric(fn_info, "PctSeqSimilarity")
    size_similarity = parse_metric(fn_info, "PctSizeSimilarity")
    if sequence_similarity is not None and sequence_similarity < 0.7:
        if direct_identity >= 0.95:
            return "benchmark_dup_resolution_mismatch"
        if size_similarity is not None and size_similarity >= 0.7:
            return "graph_allele_sequence_mismatch"
        return "graph_vcf_size_mismatch"

    ref0, ref1 = int(row["ref_start"]), int(row["ref_end"])
    length = int(row["query_end"]) - int(row["query_start"])
    dup_length = length // 2
    pad = max(500, dup_length)
    nearby_calls = [call for call in graph_calls
                    if ref0 - pad <= call["position"] <= ref1 + pad]
    compatible = [call for call in nearby_calls
                  if min(call["size"], dup_length) /
                  max(call["size"], dup_length) >= 0.7]
    if compatible:
        closest = min(compatible, key=lambda call: abs(call["position"] - ref0))
        if abs(closest["position"] - ref0) > 500:
            return "graph_vcf_breakpoint_mismatch"
        return "benchmark_matching_failure"
    if nearby_calls:
        return "graph_vcf_size_mismatch"

    if sum(correspondence_counts) == 0:
        return "assignment_missing"
    if min(correspondence_counts) == 0:
        return "one_copy_assignment_missing"
    if sum(chain_counts) == 0:
        return "chaining_missing"
    if min(chain_counts) == 0:
        return "one_copy_final_chain_missing"
    if not nearby_paf:
        return "base_alignment_missing"
    if not paf_dup_signal:
        return "paf_duplication_signal_missing"

    return "graph_vcf_call_missing"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug-dir", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--truth-coordinate-map", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    samples = ["sample01", "sample02", "sample03"]
    fn_records = read_fn_dup_records(args.benchmark_dir, samples)
    truth = [row for row in read_tsv(args.truth_coordinate_map)
             if row["svtype"] == "DUP" and row["id"] in fn_records]
    anchors = load_anchors(args.debug_dir / "grouped_anchors.tsv")
    correspondences = load_pair_anchors(
        args.debug_dir / "anchor_correspondences.tsv", anchors)
    chain_anchors = load_pair_anchors(args.debug_dir / "chain_anchors.tsv", anchors)
    paf = read_paf(args.paf)
    graph_calls = {sample: load_graph_calls(args.benchmark_dir, sample)
                   for sample in samples}
    reference_sequence = read_fasta_sequence(args.reference_fasta)

    output_rows = []
    for row in truth:
        ref0, ref1 = int(row["ref_start"]), int(row["ref_end"])
        q0, q1 = int(row["query_start"]), int(row["query_end"])
        dup_length = (q1 - q0) // 2
        copies = [(q0, q0 + dup_length), (q0 + dup_length, q1)]
        corr_counts = anchor_counts(
            correspondences, row["sample"], (ref0, ref1), copies)
        chained_counts = anchor_counts(
            chain_anchors, row["sample"], (ref0, ref1), copies)
        paf_result = paf_evidence(
            paf, row["sample"], (ref0, ref1), copies, dup_length)
        nearby_paf, copy_paf, insertion_lengths, _ = paf_result
        call_pad = max(500, dup_length // 2)
        nearby_calls = [
            call for call in graph_calls[row["sample"]]
            if ref0 - call_pad <= call["position"] <= ref1 + call_pad
        ]
        direct_identity, direct_call = best_graph_dup_identity(
            nearby_calls, reference_sequence[ref0:ref1])
        stage = classify(row, corr_counts, chained_counts, paf_result,
                         graph_calls[row["sample"]], fn_records[row["id"]],
                         direct_identity)
        output_rows.append({
            "sample": row["sample"], "sv_id": row["id"],
            "ref_start": ref0, "ref_end": ref1, "dup_length": dup_length,
            "query_start": q0, "query_end": q1,
            "copy1_assignment_count": corr_counts[0],
            "copy2_assignment_count": corr_counts[1],
            "copy1_chain_anchor_count": chained_counts[0],
            "copy2_chain_anchor_count": chained_counts[1],
            "copy1_paf_support": copy_paf[0],
            "copy2_paf_support": copy_paf[1],
            "nearby_paf_count": len(nearby_paf),
            "paf_insertion_lengths": ",".join(map(str, insertion_lengths)) or ".",
            "nearby_graph_call_count": len(nearby_calls),
            "nearby_graph_call_sizes": ",".join(
                str(call["size"]) for call in nearby_calls) or ".",
            "nearby_graph_call_positions": ",".join(
                str(call["position"]) for call in nearby_calls) or ".",
            "nearby_graph_call_ids": ",".join(
                call["id"] for call in nearby_calls) or ".",
            "best_graph_dup_identity": f"{direct_identity:.6f}",
            "best_graph_dup_call_id": direct_call["id"] if direct_call else ".",
            "truvari_pctseq": fn_records[row["id"]].get("PctSeqSimilarity", "."),
            "truvari_pctsize": fn_records[row["id"]].get("PctSizeSimilarity", "."),
            "truvari_start_distance": fn_records[row["id"]].get("StartDistance", "."),
            "truvari_end_distance": fn_records[row["id"]].get("EndDistance", "."),
            "failure_stage": stage,
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0]) if output_rows else []
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    summary = args.output.with_name(args.output.stem + ".summary.tsv")
    counts = Counter(row["failure_stage"] for row in output_rows)
    with summary.open("w") as handle:
        handle.write("failure_stage\tcount\n")
        for stage, count in sorted(counts.items()):
            handle.write(f"{stage}\t{count}\n")
    print(f"DUP FN: {len(output_rows)}")
    print(f"details: {args.output}")
    print(f"summary: {summary}")


if __name__ == "__main__":
    main()
