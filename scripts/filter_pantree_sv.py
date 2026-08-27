#!/usr/bin/env python3
"""Filter a Pantree VCF to supported structural-variant records."""

import argparse
import gzip
from collections import Counter
from pathlib import Path
from typing import Dict


def open_text(path: Path, mode: str):
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def allele_length(allele: str) -> int:
    return 0 if allele == "." else len(allele)


def parse_info(value: str) -> Dict[str, str]:
    return {
        key: item_value
        for item in value.split(";")
        if "=" in item
        for key, item_value in [item.split("=", 1)]
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Keep supported Pantree INS/DEL records above a size threshold."
    )
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--summary-tsv", required=True, type=Path)
    parser.add_argument("--min-sv-bp", type=int, default=50)
    parser.add_argument(
        "--allow-ac-zero",
        action="store_true",
        help="Retain records whose INFO/AC is zero.",
    )
    args = parser.parse_args()

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()

    with open_text(args.input_vcf, "rt") as source, open_text(
        args.output_vcf, "wt"
    ) as output:
        for line in source:
            if line.startswith("#"):
                output.write(line)
                continue

            counts["input_records"] += 1
            fields = line.rstrip("\n").split("\t")
            info = parse_info(fields[7])
            variant_type = info.get("VT", "")
            if variant_type not in {"INS", "DEL"}:
                counts["non_indel"] += 1
                continue

            alt_lengths = [allele_length(alt) for alt in fields[4].split(",")]
            size = max(
                abs(alt_length - allele_length(fields[3]))
                for alt_length in alt_lengths
            )
            if size < args.min_sv_bp:
                counts["below_min_size"] += 1
                continue

            ac_values = [
                int(value) for value in info.get("AC", "0").split(",") if value != "."
            ]
            if not args.allow_ac_zero and not any(value > 0 for value in ac_values):
                counts["ac_zero"] += 1
                continue

            output.write(line)
            counts["output_records"] += 1
            counts[variant_type] += 1

    with args.summary_tsv.open("w") as summary:
        summary.write(
            "input_vcf\tmin_sv_bp\tinput_records\toutput_records\tINS\tDEL\t"
            "non_indel\tbelow_min_size\tac_zero\n"
        )
        summary.write(
            f"{args.input_vcf}\t{args.min_sv_bp}\t{counts['input_records']}\t"
            f"{counts['output_records']}\t{counts['INS']}\t{counts['DEL']}\t"
            f"{counts['non_indel']}\t{counts['below_min_size']}\t"
            f"{counts['ac_zero']}\n"
        )


if __name__ == "__main__":
    main()
