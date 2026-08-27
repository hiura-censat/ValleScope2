#!/usr/bin/env python3
"""Assign evidence-based confidence tiers to HG008 pairwise SV clusters."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


FAMILY = {
    "minimap2": "minimizer",
    "winnowmap2": "minimizer",
    "wfmash": "mapping",
    "unialigner": "global",
    "vallescope2": "anchor_chain",
}


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", type=Path, required=True)
    args = parser.parse_args()

    metrics = read_tsv(args.analysis_dir / "alignment_metrics.tsv")
    clusters = read_tsv(args.analysis_dir / "sv_clusters.tsv")
    evidence = read_tsv(args.analysis_dir / "breakpoint_differences.tsv")
    metric_by_pair_aligner = {
        (row["pair"], row["aligner"]): row for row in metrics
    }
    evidence_by_cluster = defaultdict(list)
    for row in evidence:
        evidence_by_cluster[row["cluster_id"]].append(row)

    results = []
    for cluster in clusters:
        rows = evidence_by_cluster[cluster["cluster_id"]]
        aligners = sorted({row["aligner"] for row in rows})
        families = sorted({FAMILY[name] for name in aligners})
        sources = sorted({row["source"] for row in rows})
        usable = []
        for aligner in aligners:
            metric = metric_by_pair_aligner[(cluster["pair"], aligner)]
            if (metric["status"] == "ok"
                    and float(metric["query_coverage"]) >= 0.50
                    and float(metric["weighted_identity"]) >= 0.85):
                usable.append(aligner)
        usable_families = sorted({FAMILY[name] for name in usable})

        length = int(cluster["median_svlen"])
        bp_range = int(cluster["breakpoint_range_bp"])
        len_range = int(cluster["svlen_range_bp"])
        bp_ok = bp_range <= max(200, round(0.10 * length))
        len_ok = len_range <= max(100, round(0.20 * length))
        cigar_support = any(row["source"] in ("cigar", "strand_block") for row in rows)
        poor_pair = cluster["pair"].endswith("_vs_T2b")

        if poor_pair:
            tier = "unresolved"
            reason = "T2b has <2% pairwise coverage or no accepted full-length alignment"
        elif len(usable_families) >= 2 and bp_ok and len_ok and cigar_support:
            tier = "high"
            reason = "concordant support from >=2 usable algorithm families including CIGAR evidence"
        elif len(families) >= 2 and bp_ok and len_ok:
            tier = "moderate"
            reason = "concordant multi-family support, but coverage/identity or CIGAR evidence is limited"
        elif len(aligners) >= 2 and bp_ok and len_ok:
            tier = "moderate"
            reason = "concordant multi-aligner support from one algorithm family"
        else:
            tier = "low"
            reason = "single-aligner or discordant evidence"

        results.append({
            **cluster,
            "algorithm_family_count": len(families),
            "algorithm_families": "&".join(families),
            "usable_aligner_count": len(usable),
            "usable_aligners": "&".join(usable),
            "evidence_sources": "&".join(sources),
            "breakpoint_concordant": int(bp_ok),
            "length_concordant": int(len_ok),
            "confidence": tier,
            "confidence_reason": reason,
        })

    fields = list(results[0])
    write_tsv(args.analysis_dir / "sv_confidence.tsv", results, fields)

    summary = []
    for pair in sorted({row["pair"] for row in results}):
        for svtype in ("INS", "DEL", "INV"):
            subset = [row for row in results if row["pair"] == pair and row["svtype"] == svtype]
            counts = Counter(row["confidence"] for row in subset)
            summary.append({
                "pair": pair, "svtype": svtype, "total": len(subset),
                "high": counts["high"], "moderate": counts["moderate"],
                "low": counts["low"], "unresolved": counts["unresolved"],
            })
    write_tsv(
        args.analysis_dir / "sv_confidence_summary.tsv",
        summary,
        ["pair", "svtype", "total", "high", "moderate", "low", "unresolved"],
    )

    tiers = ("high", "moderate", "low", "unresolved")
    colors = ("#2E7D5B", "#D69C2F", "#B95756", "#73777F")
    pairs = sorted({row["pair"] for row in results})
    fig, axis = plt.subplots(figsize=(11, 5))
    bottom = [0] * len(pairs)
    for tier, color in zip(tiers, colors):
        values = [sum(row["confidence"] == tier for row in results if row["pair"] == pair)
                  for pair in pairs]
        axis.bar(pairs, values, bottom=bottom, label=tier, color=color)
        bottom = [old + new for old, new in zip(bottom, values)]
    axis.set_ylabel("SV clusters")
    axis.tick_params(axis="x", rotation=35)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.analysis_dir / "plots" / "sv_confidence_by_pair.png", dpi=180)


if __name__ == "__main__":
    main()
