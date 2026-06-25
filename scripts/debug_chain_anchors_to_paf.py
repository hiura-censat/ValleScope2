#!/usr/bin/env python3
"""Temporary helper: ValleScope2 chain_anchors.tsv -> anchor-level PAF."""

import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def load_anchor_positions(grouped_anchors):
    anchors = {}
    for row in read_tsv(grouped_anchors):
        anchors[row["anchor_id"]] = {
            "sequence_id": row["sequence_id"],
            "start": int(row["start"]),
            "end": int(row["end"]),
            "span": max(1, int(row["end"]) - int(row["start"])),
        }
    return anchors


def load_lengths(sequence_table):
    lengths = {}
    for row in read_tsv(sequence_table):
        lengths[row["sequence_id"]] = int(row["length"])
    return lengths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain-anchors", default="chain_anchors.tsv")
    parser.add_argument("--grouped-anchors", default="grouped_anchors.tsv")
    parser.add_argument("--sequences", default="vallescope2-debug/sequences.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--chain-id", default=None)
    parser.add_argument("--sample-a", default=None)
    parser.add_argument("--sample-b", default=None)
    args = parser.parse_args()

    anchors = load_anchor_positions(args.grouped_anchors)
    lengths = load_lengths(args.sequences)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out, "w") as output:
        for row in read_tsv(args.chain_anchors):
            if args.chain_id is not None and row["chain_id"] != str(args.chain_id):
                continue
            if args.sample_a and row["sample_a"] != args.sample_a:
                continue
            if args.sample_b and row["sample_b"] != args.sample_b:
                continue

            target = anchors[row["anchor_a"]]
            query = anchors[row["anchor_b"]]
            qname = query["sequence_id"]
            tname = target["sequence_id"]
            qlen = lengths.get(qname, query["end"])
            tlen = lengths.get(tname, target["end"])
            strand = row["assign_strand"] if row["assign_strand"] in "+-" else "+"
            alnlen = max(query["span"], target["span"])
            nmatch = min(query["span"], target["span"])
            score = float(row["candidate_score"])
            mapq = 60 if row["support_direction"] == "both" else 30
            tags = [
                f"as:f:{score}",
                f"sd:Z:{row['support_direction']}",
                f"ch:i:{row['chain_id']}",
                f"rk:i:{row['rank']}",
                f"ab:Z:{row['anchor_a']},{row['anchor_b']}",
            ]
            output.write(
                "\t".join(map(str, [
                    qname, qlen, query["start"], query["end"], strand,
                    tname, tlen, target["start"], target["end"],
                    nmatch, alnlen, mapq, *tags,
                ])) + "\n"
            )


if __name__ == "__main__":
    main()
