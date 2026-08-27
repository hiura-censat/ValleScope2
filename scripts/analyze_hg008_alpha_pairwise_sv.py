#!/usr/bin/env python3
"""Normalize HG008 pairwise PAFs, call coarse SVs, and compare aligners."""

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


CIGAR_RE = re.compile(r"(\d+)([=XMID])")
ALIGNERS = ("minimap2", "winnowmap2", "wfmash", "unialigner", "vallescope2")


def tags(fields):
    result = {}
    for value in fields[12:]:
        parts = value.split(":", 2)
        if len(parts) == 3:
            result[parts[0]] = parts[2]
    return result


def read_paf(path):
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                continue
            tag = tags(f)
            rows.append({
                "qname": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]),
                "strand": f[4], "tname": f[5], "tlen": int(f[6]),
                "ts": int(f[7]), "te": int(f[8]), "matches": int(f[9]),
                "block": int(f[10]), "mapq": int(f[11]), "tags": tag,
            })
    return rows


def overlap_fraction(a0, a1, b0, b1):
    overlap = max(0, min(a1, b1) - max(a0, b0))
    return overlap / max(1, min(a1 - a0, b1 - b0))


def select_records(rows, aligner):
    if aligner in ("minimap2", "winnowmap2"):
        primary = [row for row in rows if row["tags"].get("tp", "P") != "S"]
        return primary or rows
    if aligner == "unialigner":
        return [row for row in rows if row["block"] and row["matches"] / row["block"] >= 0.90]
    if aligner != "wfmash":
        return rows
    selected = []
    for row in sorted(rows, key=lambda x: (x["mapq"], x["block"]), reverse=True):
        conflict = False
        for old in selected:
            if (overlap_fraction(row["qs"], row["qe"], old["qs"], old["qe"]) > 0.5
                    or overlap_fraction(row["ts"], row["te"], old["ts"], old["te"]) > 0.5):
                conflict = True
                break
        if not conflict:
            selected.append(row)
    return selected


def union(intervals):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def gap_stats(merged, length):
    gaps = []
    cursor = 0
    for start, end in merged:
        if start > cursor:
            gaps.append(start - cursor)
        cursor = max(cursor, end)
    if cursor < length:
        gaps.append(length - cursor)
    return len(gaps), max(gaps or [0]), sum(x >= 1000 for x in gaps), sum(x >= 10000 for x in gaps)


def add_event(events, pair, aligner, kind, pos, length, source):
    if length < 50:
        return
    events.append({
        "pair": pair, "aligner": aligner, "svtype": kind,
        "target_pos": max(0, int(pos)), "svlen": int(length), "source": source,
    })


def cigar_events(row, pair, aligner, events):
    cigar = row["tags"].get("cg")
    if not cigar:
        return
    tpos = row["ts"]
    qpos = row["qs"] if row["strand"] == "+" else row["qe"]
    for length_text, operation in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if operation in "=XM":
            tpos += length
            qpos += length if row["strand"] == "+" else -length
        elif operation == "D":
            add_event(events, pair, aligner, "DEL", tpos, length, "cigar")
            tpos += length
        elif operation == "I":
            add_event(events, pair, aligner, "INS", tpos, length, "cigar")
            qpos += length if row["strand"] == "+" else -length


def call_events(rows, pair, aligner):
    events = []
    if not rows:
        return events
    strand_mass = Counter()
    for row in rows:
        strand_mass[row["strand"]] += row["block"]
        cigar_events(row, pair, aligner, events)
    main_strand = strand_mass.most_common(1)[0][0]
    for row in rows:
        if row["strand"] != main_strand:
            add_event(events, pair, aligner, "INV", row["ts"], row["te"] - row["ts"], "strand_block")

    main = sorted((row for row in rows if row["strand"] == main_strand), key=lambda x: x["ts"])
    for left, right in zip(main, main[1:]):
        tgap = right["ts"] - left["te"]
        if main_strand == "+":
            qgap = right["qs"] - left["qe"]
        else:
            qgap = left["qs"] - right["qe"]
        if tgap < 0 or qgap < 0:
            continue
        difference = qgap - tgap
        if difference >= 50:
            add_event(events, pair, aligner, "INS", left["te"], difference, "interblock")
        elif difference <= -50:
            add_event(events, pair, aligner, "DEL", left["te"], -difference, "interblock")

    dedup = {}
    for event in events:
        key = (event["svtype"], round(event["target_pos"] / 20), round(event["svlen"] / 20))
        dedup.setdefault(key, event)
    return list(dedup.values())


def compatible(a, b):
    if a["pair"] != b["pair"] or a["svtype"] != b["svtype"]:
        return False
    pos_tolerance = max(200, int(0.1 * max(a["svlen"], b["svlen"])))
    length_tolerance = max(100, int(0.2 * max(a["svlen"], b["svlen"])))
    return (abs(a["target_pos"] - b["target_pos"]) <= pos_tolerance
            and abs(a["svlen"] - b["svlen"]) <= length_tolerance)


def cluster_events(events):
    clusters = []
    for event in sorted(events, key=lambda x: (x["pair"], x["svtype"], x["target_pos"], x["svlen"])):
        match = None
        for cluster in reversed(clusters):
            if cluster["pair"] != event["pair"] or cluster["svtype"] != event["svtype"]:
                continue
            representative = cluster["events"][0]
            if compatible(representative, event):
                match = cluster
                break
            if representative["target_pos"] < event["target_pos"] - max(10000, event["svlen"]):
                break
        if match is None:
            match = {"pair": event["pair"], "svtype": event["svtype"], "events": []}
            clusters.append(match)
        match["events"].append(event)
    return clusters


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    out = args.root / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    plots = out / "plots"
    plots.mkdir(exist_ok=True)

    pairs = list(csv.DictReader((args.root / "pairs.tsv").open(), delimiter="\t"))
    metrics = []
    all_events = []
    selected_cache = {}
    for pair_row in pairs:
        pair = pair_row["pair"]
        for aligner in ALIGNERS:
            paf = args.root / "pairs" / pair / aligner / "alignment.paf"
            raw_rows = read_paf(paf)
            rows = select_records(raw_rows, aligner)
            selected_cache[(pair, aligner)] = rows
            status = "ok" if rows else (
                "failed_full_length" if aligner == "vallescope2"
                else "low_identity_global" if aligner == "unialigner" and raw_rows
                else "no_alignment")
            qlen = max((row["qlen"] for row in rows), default=0)
            tlen = max((row["tlen"] for row in rows), default=0)
            qmerged = union((row["qs"], row["qe"]) for row in rows)
            tmerged = union((row["ts"], row["te"]) for row in rows)
            qcov = sum(end - start for start, end in qmerged)
            tcov = sum(end - start for start, end in tmerged)
            qg = gap_stats(qmerged, qlen) if qlen else (0, 0, 0, 0)
            tg = gap_stats(tmerged, tlen) if tlen else (0, 0, 0, 0)
            identity_num = sum(row["matches"] for row in rows)
            identity_den = sum(row["block"] for row in rows)
            events = call_events(rows, pair, aligner)
            all_events.extend(events)
            metrics.append({
                "pair": pair, "aligner": aligner, "status": status, "blocks": len(rows),
                "query_length": qlen, "query_covered_bp": qcov,
                "query_coverage": qcov / qlen if qlen else 0,
                "query_gap_count": qg[0], "query_largest_gap_bp": qg[1],
                "query_gaps_ge_1kb": qg[2], "query_gaps_ge_10kb": qg[3],
                "target_length": tlen, "target_covered_bp": tcov,
                "target_coverage": tcov / tlen if tlen else 0,
                "target_gap_count": tg[0], "target_largest_gap_bp": tg[1],
                "target_gaps_ge_1kb": tg[2], "target_gaps_ge_10kb": tg[3],
                "weighted_identity": identity_num / identity_den if identity_den else 0,
                "sv_count": len(events),
            })

    metric_fields = list(metrics[0])
    event_fields = ["pair", "aligner", "svtype", "target_pos", "svlen", "source"]
    write_tsv(out / "alignment_metrics.tsv", metrics, metric_fields)
    write_tsv(out / "sv_candidates.tsv", all_events, event_fields)

    clusters = cluster_events(all_events)
    cluster_rows = []
    breakpoint_rows = []
    set_counts = Counter()
    for index, cluster in enumerate(clusters, 1):
        aligners = sorted(set(event["aligner"] for event in cluster["events"]))
        positions = [event["target_pos"] for event in cluster["events"]]
        lengths = [event["svlen"] for event in cluster["events"]]
        key = "&".join(aligners)
        set_counts[key] += 1
        cluster_rows.append({
            "cluster_id": index, "pair": cluster["pair"], "svtype": cluster["svtype"],
            "median_target_pos": sorted(positions)[len(positions) // 2],
            "median_svlen": sorted(lengths)[len(lengths) // 2],
            "aligner_count": len(aligners), "aligners": key,
            "breakpoint_range_bp": max(positions) - min(positions),
            "svlen_range_bp": max(lengths) - min(lengths),
        })
        for event in cluster["events"]:
            breakpoint_rows.append({
                "cluster_id": index, **event,
                "delta_from_min_breakpoint": event["target_pos"] - min(positions),
            })
    cluster_fields = [
        "cluster_id", "pair", "svtype", "median_target_pos", "median_svlen",
        "aligner_count", "aligners", "breakpoint_range_bp", "svlen_range_bp",
    ]
    write_tsv(out / "sv_clusters.tsv", cluster_rows, cluster_fields)
    write_tsv(out / "breakpoint_differences.tsv", breakpoint_rows,
              ["cluster_id"] + event_fields + ["delta_from_min_breakpoint"])

    count_rows = []
    for pair_row in pairs:
        for aligner in ALIGNERS:
            subset = [event for event in all_events
                      if event["pair"] == pair_row["pair"] and event["aligner"] == aligner]
            counts = Counter(event["svtype"] for event in subset)
            count_rows.append({
                "pair": pair_row["pair"], "aligner": aligner,
                "total": len(subset), "INS": counts["INS"],
                "DEL": counts["DEL"], "INV": counts["INV"],
            })
    write_tsv(out / "sv_counts.tsv", count_rows,
              ["pair", "aligner", "total", "INS", "DEL", "INV"])

    assignment_rows = []
    for tumor in ("T1", "T2a", "T2b"):
        for aligner in ALIGNERS:
            candidates = []
            for normal in ("N1", "N2"):
                pair = normal + "_vs_" + tumor
                metric = next((row for row in metrics
                               if row["pair"] == pair and row["aligner"] == aligner), None)
                if metric:
                    score = (2 * metric["query_coverage"] * metric["target_coverage"]
                             / max(1e-12, metric["query_coverage"] + metric["target_coverage"])
                             * metric["weighted_identity"])
                    candidates.append((score, normal, metric))
            if candidates:
                score, normal, metric = max(candidates)
                assignment_rows.append({
                    "tumor_segment": tumor, "aligner": aligner,
                    "best_normal": normal, "score": score,
                    "query_coverage": metric["query_coverage"],
                    "target_coverage": metric["target_coverage"],
                    "weighted_identity": metric["weighted_identity"],
                })
    write_tsv(out / "best_normal_assignment.tsv", assignment_rows,
              ["tumor_segment", "aligner", "best_normal", "score",
               "query_coverage", "target_coverage", "weighted_identity"])

    # SV counts by type.
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for axis, kind in zip(axes, ("INS", "DEL", "INV")):
        values = [sum(1 for event in all_events
                      if event["aligner"] == aligner and event["svtype"] == kind)
                  for aligner in ALIGNERS]
        axis.bar(ALIGNERS, values, color=["#3366A8", "#D97732", "#4C9A68", "#8B6BB1", "#C34A4A"])
        axis.set_title(kind)
        axis.tick_params(axis="x", rotation=35)
        axis.set_ylabel("SV candidates")
    fig.tight_layout()
    fig.savefig(plots / "sv_count_by_type.png", dpi=180)
    plt.close(fig)

    # SV length distribution.
    fig, axis = plt.subplots(figsize=(9, 5))
    for aligner in ALIGNERS:
        lengths = [event["svlen"] for event in all_events if event["aligner"] == aligner]
        if lengths:
            axis.hist(lengths, bins=40, histtype="step", linewidth=1.5,
                      label=aligner, log=True)
    axis.set_xscale("log")
    axis.set_xlabel("SV length (bp)")
    axis.set_ylabel("Count")
    axis.legend()
    fig.tight_layout()
    fig.savefig(plots / "sv_length_distribution.png", dpi=180)
    plt.close(fig)

    # Reference-position histograms.
    position_max = max((row["target_length"] for row in metrics), default=1)
    fig, axes = plt.subplots(len(ALIGNERS), 1, figsize=(12, 11), sharex=True)
    for axis, aligner in zip(axes, ALIGNERS):
        positions = [event["target_pos"] for event in all_events if event["aligner"] == aligner]
        axis.hist(positions, bins=60, range=(0, position_max), color="#497AA5")
        axis.set_ylabel(aligner)
    axes[-1].set_xlabel("Normal alpha-satellite coordinate (bp)")
    fig.tight_layout()
    fig.savefig(plots / "sv_position_histogram.png", dpi=180)
    plt.close(fig)

    # UpSet-like intersection bars.
    fig, axis = plt.subplots(figsize=(11, 6))
    labels, values = zip(*set_counts.most_common()) if set_counts else ([], [])
    axis.bar(range(len(values)), values, color="#3D6F8E")
    axis.set_xticks(range(len(labels)))
    axis.set_xticklabels(labels, rotation=55, ha="right")
    axis.set_ylabel("SV clusters")
    axis.set_title("Aligner support combinations")
    fig.tight_layout()
    fig.savefig(plots / "sv_aligner_upset.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
