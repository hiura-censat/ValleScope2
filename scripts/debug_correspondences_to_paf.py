#!/usr/bin/env python3
"""Temporary helper: ValleScope2 anchor_correspondences.tsv -> anchor-level PAF."""

import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def load_anchor_positions(grouped_anchors):
    anchors = {}
    for row in read_tsv(grouped_anchors):
        start = int(row["start"])
        end = int(row["end"])
        anchors[row["anchor_id"]] = {
            "sequence_id": row["sequence_id"],
            "start": start,
            "end": end,
            "span": max(1, end - start),
        }
    return anchors


def load_lengths(sequence_table):
    lengths = {}
    with open(sequence_table, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            lengths[row["sequence_id"]] = int(row["length"])
    return lengths


def optional_float(text):
    if text == ".":
        return None
    return float(text)


def choose_score_and_strand(row):
    f_score = optional_float(row["forward_score"])
    r_score = optional_float(row["reverse_score"])
    f_strand = row["forward_strand"] if row["forward_strand"] in "+-" else "+"
    r_strand = row["reverse_strand"] if row["reverse_strand"] in "+-" else "+"
    support = row["support_direction"]

    if support == "forward_only":
        return f_score, f_strand
    if support == "reverse_only":
        return r_score, r_strand
    if r_score is not None and (f_score is None or r_score > f_score):
        return r_score, r_strand
    return f_score, f_strand


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--correspondences", default="anchor_correspondences.tsv")
    parser.add_argument("--grouped-anchors", default="grouped_anchors.tsv")
    parser.add_argument("--sequences", default="vallescope2-debug/sequences.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--sample-a", default=None)
    parser.add_argument("--sample-b", default=None)
    parser.add_argument("--support", default=None,
                        help="comma-separated support filters, e.g. both,forward_only")
    args = parser.parse_args()

    anchors = load_anchor_positions(args.grouped_anchors)
    lengths = load_lengths(args.sequences)
    support_filter = set(args.support.split(",")) if args.support else None

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as output:
        for row in read_tsv(args.correspondences):
            if args.sample_a and row["sample_a"] != args.sample_a:
                continue
            if args.sample_b and row["sample_b"] != args.sample_b:
                continue
            if support_filter and row["support_direction"] not in support_filter:
                continue

            target = anchors[row["anchor_a"]]
            query = anchors[row["anchor_b"]]
            score, strand = choose_score_and_strand(row)
            if score is None:
                continue

            qname = query["sequence_id"]
            tname = target["sequence_id"]
            qlen = lengths.get(qname, query["end"])
            tlen = lengths.get(tname, target["end"])
            qstart, qend = query["start"], query["end"]
            tstart, tend = target["start"], target["end"]
            alnlen = max(query["span"], target["span"])
            nmatch = min(query["span"], target["span"])
            mapq = 60 if row["support_direction"] == "both" else 30

            tags = [
                f"as:f:{score}",
                f"sd:Z:{row['support_direction']}",
                f"ab:Z:{row['anchor_a']},{row['anchor_b']}",
            ]
            output.write(
                "\t".join(map(str, [
                    qname, qlen, qstart, qend, strand,
                    tname, tlen, tstart, tend, nmatch, alnlen, mapq,
                    *tags,
                ])) + "\n"
            )


if __name__ == "__main__":
    main()
