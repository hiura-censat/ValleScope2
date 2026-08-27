#!/usr/bin/env python3
"""Project genomic CDR intervals through a minimap2 PAF CIGAR."""

import argparse
import re
from pathlib import Path


CIGAR = re.compile(r"(\d+)([=XMID])")


def paf_tags(fields):
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--cdr-bed", type=Path, required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    intervals = []
    with args.cdr_bed.open() as handle:
        for line in handle:
            chrom, start, end, label = line.rstrip().split("\t")[:4]
            if chrom == args.chromosome:
                intervals.append((int(start), int(end), label))

    projected = []
    with args.paf.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12 or fields[5] != args.chromosome:
                continue
            tags = paf_tags(fields)
            cigar = tags.get("cg")
            if not cigar:
                continue
            strand = fields[4]
            qpos = int(fields[2]) if strand == "+" else int(fields[3])
            tpos = int(fields[7])
            for length_text, operation in CIGAR.findall(cigar):
                length = int(length_text)
                if operation in "=XM":
                    for start, end, label in intervals:
                        ov_start = max(start, tpos)
                        ov_end = min(end, tpos + length)
                        if ov_start >= ov_end:
                            continue
                        if strand == "+":
                            q_start = qpos + ov_start - tpos
                            q_end = qpos + ov_end - tpos
                        else:
                            q_start = qpos - (ov_end - tpos)
                            q_end = qpos - (ov_start - tpos)
                        projected.append(
                            (args.chromosome, q_start, q_end, label, start, end)
                        )
                    tpos += length
                    qpos += length if strand == "+" else -length
                elif operation == "I":
                    qpos += length if strand == "+" else -length
                elif operation == "D":
                    tpos += length

    projected.sort(key=lambda row: (row[3], row[1], row[2]))
    merged = []
    for row in projected:
        if (merged and row[3] == merged[-1][3]
                and row[1] <= merged[-1][2] + 1):
            merged[-1] = (
                row[0], merged[-1][1], max(merged[-1][2], row[2]),
                row[3], row[4], row[5],
            )
        else:
            merged.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for row in merged:
            handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
