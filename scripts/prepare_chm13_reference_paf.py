#!/usr/bin/env python3
"""Orient all-to-all PAF records as sample query versus one reference target."""

import argparse
import re
from pathlib import Path


def short_name(name):
    return name.split("|", 1)[0]


def swap_cigar(cigar):
    return re.sub(r"[ID]", lambda match: "D" if match.group() == "I" else "I", cigar)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-token", default="chm13")
    parser.add_argument("--reference-name", default="chm13v2.0")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with args.input.open() as source, args.output.open("w") as dest:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            query_is_ref = args.reference_token.lower() in fields[0].lower()
            target_is_ref = args.reference_token.lower() in fields[5].lower()
            if query_is_ref == target_is_ref:
                continue
            if query_is_ref:
                fields[:12] = [
                    fields[5], fields[6], fields[7], fields[8], fields[4],
                    fields[0], fields[1], fields[2], fields[3],
                    fields[9], fields[10], fields[11],
                ]
                for index in range(12, len(fields)):
                    if fields[index].startswith("cg:Z:"):
                        fields[index] = "cg:Z:" + swap_cigar(fields[index][5:])
            fields[0] = short_name(fields[0])
            fields[5] = args.reference_name
            dest.write("\t".join(fields) + "\n")
            written += 1
    print(f"records\t{written}")


if __name__ == "__main__":
    main()
