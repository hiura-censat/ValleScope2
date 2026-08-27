#!/usr/bin/env python3
"""Convert coarse SV candidates to VCF and compare aligners with Truvari."""

import argparse
import csv
import itertools
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


ALIGNERS = ("minimap2", "winnowmap2", "wfmash", "unialigner", "vallescope2")


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_vcf(path, pair, target_length, events):
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##contig=<ID={pair},length={target_length}>\n")
        handle.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        handle.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        handle.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        ordered = sorted(events, key=lambda x: (int(x["target_pos"]), x["svtype"], int(x["svlen"])))
        for index, event in enumerate(ordered, 1):
            pos = int(event["target_pos"]) + 1
            length = int(event["svlen"])
            svtype = event["svtype"]
            end = pos if svtype == "INS" else min(target_length, pos + length - 1)
            svlen = -length if svtype == "DEL" else length
            info = f"SVTYPE={svtype};SVLEN={svlen};END={end}"
            handle.write(
                f"{pair}\t{pos}\t{pair}_{index}\tN\t<{svtype}>\t.\tPASS\t"
                f"{info}\tGT\t0/1\n"
            )


def metric(summary, key):
    value = summary.get(key, 0)
    return float(value) if value is not None else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--truvari", type=Path, required=True)
    parser.add_argument("--bgzip", type=Path, required=True)
    parser.add_argument("--tabix", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    candidates = read_tsv(args.analysis_dir / "sv_candidates.tsv")
    metrics = read_tsv(args.analysis_dir / "alignment_metrics.tsv")
    lengths = {}
    for row in metrics:
        length = int(row["target_length"])
        if length:
            lengths[row["pair"]] = max(lengths.get(row["pair"], 0), length)
    pairs = sorted({row["pair"] for row in candidates} | set(lengths))
    grouped = defaultdict(list)
    for row in candidates:
        grouped[(row["pair"], row["aligner"])].append(row)

    vcf_root = args.output_dir / "vcfs"
    bench_root = args.output_dir / "bench"
    vcf_root.mkdir(exist_ok=True)
    bench_root.mkdir(exist_ok=True)
    vcfs = {}
    for pair in pairs:
        pair_dir = vcf_root / pair
        pair_dir.mkdir(exist_ok=True)
        for aligner in ALIGNERS:
            plain = pair_dir / f"{aligner}.vcf"
            write_vcf(plain, pair, lengths[pair], grouped[(pair, aligner)])
            gz = Path(str(plain) + ".gz")
            with gz.open("wb") as output:
                subprocess.run([args.bgzip, "-c", plain], stdout=output, check=True)
            subprocess.run([args.tabix, "-f", "-p", "vcf", gz], check=True)
            vcfs[(pair, aligner)] = gz

    rows = []
    for pair in pairs:
        for base, comp in itertools.combinations(ALIGNERS, 2):
            out = bench_root / pair / f"{base}_vs_{comp}"
            if out.exists():
                shutil.rmtree(out)
            out.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run([
                args.truvari, "bench",
                "-b", vcfs[(pair, base)],
                "-c", vcfs[(pair, comp)],
                "-o", out,
                "-r", "500",
                "-p", "0",
                "-P", "0.7",
                "-O", "0",
                "-N",
                "-m", "0",
                "-s", "50",
                "-S", "50",
                "--sizemax", "-1",
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            summary = json.loads((out / "summary.json").read_text())
            tp_base = int(metric(summary, "TP-base"))
            tp_call = int(metric(summary, "TP-comp"))
            fp = int(metric(summary, "FP"))
            fn = int(metric(summary, "FN"))
            matched = min(tp_base, tp_call)
            base_total = tp_base + fn
            comp_total = tp_call + fp
            union = base_total + comp_total - matched
            rows.append({
                "pair": pair, "base_method": base, "comparison_method": comp,
                "tp_base": tp_base, "tp_call": tp_call, "fp": fp, "fn": fn,
                "base_recall": metric(summary, "recall"),
                "comparison_precision": metric(summary, "precision"),
                "f1": metric(summary, "f1"),
                "jaccard": matched / union if union else 1.0,
                "base_total": base_total, "comparison_total": comp_total,
            })

    fields = list(rows[0])
    with (args.output_dir / "pairwise_truvari.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    aggregate = []
    for include_t2b in (True, False):
        label = "all_pairs" if include_t2b else "exclude_T2b"
        for base, comp in itertools.combinations(ALIGNERS, 2):
            subset = [
                row for row in rows
                if row["base_method"] == base and row["comparison_method"] == comp
                and (include_t2b or not row["pair"].endswith("_vs_T2b"))
            ]
            tp_base = sum(row["tp_base"] for row in subset)
            tp_call = sum(row["tp_call"] for row in subset)
            fp = sum(row["fp"] for row in subset)
            fn = sum(row["fn"] for row in subset)
            matched = min(tp_base, tp_call)
            precision = tp_call / (tp_call + fp) if tp_call + fp else 0
            recall = tp_base / (tp_base + fn) if tp_base + fn else 0
            f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0
            union = tp_base + fn + tp_call + fp - matched
            aggregate.append({
                "scope": label, "base_method": base, "comparison_method": comp,
                "tp_base": tp_base, "tp_call": tp_call, "fp": fp, "fn": fn,
                "base_recall": recall, "comparison_precision": precision, "f1": f1,
                "jaccard": matched / union if union else 1.0,
            })
    fields = list(aggregate[0])
    with (args.output_dir / "aggregate_truvari.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(aggregate)


if __name__ == "__main__":
    main()
