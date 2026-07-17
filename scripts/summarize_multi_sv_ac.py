#!/usr/bin/env python3
"""Summarize SV counts and AC/AN/AF across multiple normalized VCFs."""

import argparse
import gzip
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_vcf(path):
    summary = Counter()
    keys = set()
    sizes, ac_values, an_values, af_values = [], [], [], []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            delta = len(fields[4]) - len(fields[3])
            summary["records"] += 1
            summary["INS" if delta > 0 else "DEL"] += 1
            keys.add((fields[0], fields[1], fields[3].upper(), fields[4].upper()))
            sizes.append(abs(delta))
            if "AC" in info:
                ac_values.append(int(info["AC"].split(",")[0]))
            if "AN" in info:
                an_values.append(int(info["AN"]))
            if "AF" in info:
                af_values.append(float(info["AF"].split(",")[0]))
    sizes.sort()
    return {
        "records": summary["records"], "unique": len(keys),
        "duplicates": summary["records"] - len(keys), "INS": summary["INS"],
        "DEL": summary["DEL"], "median_size": sizes[len(sizes) // 2],
        "max_size": max(sizes), "mean_ac": sum(ac_values) / len(ac_values),
        "mean_an": sum(an_values) / len(an_values),
        "mean_af": sum(af_values) / len(af_values), "ac": Counter(ac_values),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", action="append", type=Path, required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if len(args.vcf) != len(args.label):
        parser.error("--vcf and --label counts must match")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = [read_vcf(path) for path in args.vcf]

    with (args.output_dir / "sv_ac_summary.tsv").open("w") as handle:
        handle.write("method\trecords\tunique_alleles\tduplicate_records\tINS\tDEL\tmedian_sv_bp\tmax_sv_bp\tmean_AC\tmean_AN\tmean_AF\n")
        for label, item in zip(args.label, data):
            handle.write(
                f"{label}\t{item['records']}\t{item['unique']}\t{item['duplicates']}\t"
                f"{item['INS']}\t{item['DEL']}\t{item['median_size']}\t{item['max_size']}\t"
                f"{item['mean_ac']:.6f}\t{item['mean_an']:.6f}\t{item['mean_af']:.6f}\n"
            )
    with (args.output_dir / "ac_distribution.tsv").open("w") as handle:
        handle.write("AC\t" + "\t".join(args.label) + "\n")
        for ac in range(1, 10):
            handle.write(str(ac) + "\t" + "\t".join(str(item["ac"][ac]) for item in data) + "\n")

    x = np.arange(1, 10)
    width = 0.24
    colors = ["#2F6B9A", "#D97732", "#3A8D68"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=180)
    for index, (label, item) in enumerate(zip(args.label, data)):
        axes[0].bar(x + (index - 1) * width, [item["ac"][ac] for ac in x], width,
                    label=label, color=colors[index])
        axes[1].bar(x + (index - 1) * width,
                    [item["ac"][ac] / item["records"] for ac in x], width,
                    label=label, color=colors[index])
    axes[0].set_title("AC distribution: SV record counts")
    axes[0].set_ylabel("SV records")
    axes[1].set_title("AC distribution: within-method fraction")
    axes[1].set_ylabel("Fraction of SV records")
    for axis in axes:
        axis.set_xlabel("Alternative allele count (AC)")
        axis.set_xticks(x)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output_dir / "ac_distribution.png", bbox_inches="tight")


if __name__ == "__main__":
    main()
