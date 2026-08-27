#!/usr/bin/env python3
"""Keep vg deconstruct records containing an allele length change >= threshold."""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-sv-bp", type=int, default=50)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open() as source, args.output.open("w") as output:
        for line in source:
            if line.startswith("#"):
                output.write(line)
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 5:
                continue
            ref_length = len(fields[3])
            if any(abs(len(alt) - ref_length) >= args.min_sv_bp
                   for alt in fields[4].split(",") if not alt.startswith("<")):
                output.write(line)


if __name__ == "__main__":
    main()
