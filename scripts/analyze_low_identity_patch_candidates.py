#!/usr/bin/env python3
"""Characterize low-identity ValleScope2 patch CIGARs and threshold sweeps."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([=XIDM])")
GAP_BINS = [0, 500, 1000, 2000, 5000, 10000, 20000, 40000, 70001]
GAP_LABELS = ["<0.5k", "0.5-1k", "1-2k", "2-5k", "5-10k", "10-20k", "20-40k", "40-70k"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gaps", required=True, type=Path)
    parser.add_argument("--patch-intervals", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cigar_segments(cigar: str) -> tuple[list[dict[str, float]], list[tuple[str, int]]]:
    segments: list[dict[str, float]] = []
    long_indels: list[tuple[str, int]] = []
    aligned_bp = matches = short_indel_bp = 0
    for length_text, operation in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if operation in "ID" and length >= 500:
            segments.append(
                {
                    "aligned_bp": aligned_bp,
                    "matches": matches,
                    "short_indel_bp": short_indel_bp,
                }
            )
            long_indels.append((operation, length))
            aligned_bp = matches = short_indel_bp = 0
        elif operation in "=XM":
            aligned_bp += length
            if operation == "=":
                matches += length
        elif operation in "ID":
            aligned_bp += length
            short_indel_bp += length
    segments.append(
        {
            "aligned_bp": aligned_bp,
            "matches": matches,
            "short_indel_bp": short_indel_bp,
        }
    )
    return segments, long_indels


def segment_identity(segment: dict[str, float]) -> float:
    return segment["matches"] / segment["aligned_bp"] if segment["aligned_bp"] else 0.0


def single_pass(row: dict[str, str], flank: int, identity: float, extra: int, density: float) -> bool:
    return (
        int(row["rescue_long_indel_count"]) == 1
        and min(int(row["rescue_left_bp"]), int(row["rescue_right_bp"])) >= flank
        and min(float(row["rescue_left_identity"]), float(row["rescue_right_identity"])) >= identity
        and int(row["rescue_extra_indel_bp"]) <= extra
        and min(float(row["left_endpoint_identity"]), float(row["right_endpoint_identity"])) >= 0.85
        and float(row["max_short_error_density"]) <= density
    )


def multi_pass(
    row: dict[str, str], min_segment_bp: int, segment_identity_min: float,
    short_indel_cap: int, density: float, operation_mode: str,
) -> bool:
    segments, long_indels = cigar_segments(row["cigar"])
    if len(long_indels) <= 1:
        return False
    operations = {operation for operation, _ in long_indels}
    if operation_mode == "insertions" and operations != {"I"}:
        return False
    if operation_mode == "same_operation" and len(operations) != 1:
        return False
    return (
        all(
            segment["aligned_bp"] >= min_segment_bp
            and segment_identity(segment) >= segment_identity_min
            for segment in segments
        )
        and sum(int(segment["short_indel_bp"]) for segment in segments) <= short_indel_cap
        and min(float(row["left_endpoint_identity"]), float(row["right_endpoint_identity"])) >= 0.85
        and float(row["max_short_error_density"]) <= density
    )


def gap_bin(value: int) -> str:
    for lower, upper, label in zip(GAP_BINS, GAP_BINS[1:], GAP_LABELS):
        if lower <= value < upper:
            return label
    return ">=70k"


def quantile(values: list[int], fraction: float) -> int:
    ordered = sorted(values)
    return ordered[round((len(ordered) - 1) * fraction)]


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_distributions(
    output: Path, features: list[dict[str, object]], profiles: list[dict[str, object]]
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    class_order = [
        "clean_single_long_indel",
        "clean_multi_long_indel",
        "near_global_identity_threshold",
        "noisy_or_weak_flank_long_indel",
        "diffuse_or_no_long_indel",
    ]
    colors = {
        "clean_single_long_indel": "#2a9d8f",
        "clean_multi_long_indel": "#3a86ff",
        "near_global_identity_threshold": "#e9c46a",
        "noisy_or_weak_flank_long_indel": "#e76f51",
        "diffuse_or_no_long_indel": "#6c757d",
    }
    labels = {
        "clean_single_long_indel": "clean single long I/D",
        "clean_multi_long_indel": "clean multi long I/D",
        "near_global_identity_threshold": "near global threshold",
        "noisy_or_weak_flank_long_indel": "weak-flank/noisy long I/D",
        "diffuse_or_no_long_indel": "diffuse/no long I/D",
    }

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))
    x = np.arange(len(GAP_LABELS))
    bottom = np.zeros(len(GAP_LABELS))
    for candidate_class in class_order:
        values = np.array([
            sum(
                row["gap_bin"] == gap_label and row["candidate_class"] == candidate_class
                for row in features
            )
            for gap_label in GAP_LABELS
        ])
        axes[0].bar(
            x, values, bottom=bottom, color=colors[candidate_class],
            label=labels[candidate_class], width=0.78,
        )
        bottom += values
    axes[0].set_xticks(x, GAP_LABELS, rotation=40, ha="right")
    axes[0].set_ylabel("Low-identity patch count")
    axes[0].set_xlabel("max(ref gap, query gap)")
    axes[0].set_title("A. Gap length distribution")
    axes[0].legend(frameon=False, fontsize=8)

    identity_bins = np.linspace(0, 0.95, 20)
    for long_group, color, label in [
        (0, "#6c757d", "no long I/D"),
        (1, "#2a9d8f", "single long I/D"),
        (2, "#3a86ff", "multiple long I/D"),
    ]:
        values = [
            float(row["global_identity"])
            for row in features
            if (int(row["long_indel_count"]) if long_group < 2 else min(2, int(row["long_indel_count"]))) == long_group
        ]
        axes[1].hist(values, bins=identity_bins, alpha=0.72, color=color, label=label)
    axes[1].axvline(0.95, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Global WFA identity")
    axes[1].set_ylabel("Count")
    axes[1].set_title("B. Identity by CIGAR structure")
    axes[1].legend(frameon=False, fontsize=8)

    selected_names = [
        "single_current", "single_flank_100", "single_relaxed",
        "multi_current_insertions", "multi_identity_0.99",
        "multi_identity_0.98", "multi_identity_0.95",
        "global_identity_0.945", "global_identity_0.940",
    ]
    selected = [row for name in selected_names for row in profiles if row["profile"] == name]
    y = np.arange(len(selected))
    bar_colors = [
        "#2a9d8f" if row["scope"] == "single_long_indel"
        else "#3a86ff" if row["scope"] == "multi_long_indel"
        else "#e9c46a"
        for row in selected
    ]
    counts = [int(row["recoverable"]) for row in selected]
    axes[2].barh(y, counts, color=bar_colors)
    axes[2].set_yticks(y, [str(row["profile"]) for row in selected], fontsize=8)
    axes[2].invert_yaxis()
    axes[2].set_xlabel("Potentially recoverable patches")
    axes[2].set_title("C. Threshold sweep")
    for index, count in enumerate(counts):
        axes[2].text(count + 0.4, index, str(count), va="center", fontsize=8)
    axes[2].set_xlim(0, max(counts) * 1.18 + 1)

    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#dddddd", linewidth=0.5, alpha=0.7)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    patches = {
        (row["left_chain_id"], row["right_chain_id"]): row
        for row in read_tsv(args.patch_intervals)
    }
    gap_rows = [row for row in read_tsv(args.gaps) if row["cause"] == "patch_interval_low_identity"]

    features: list[dict[str, object]] = []
    source_rows: list[dict[str, str]] = []
    for gap in gap_rows:
        patch = patches[(gap["patch_left_chain_id"], gap["patch_right_chain_id"])]
        source_rows.append(patch)
        segments, long_indels = cigar_segments(patch["cigar"])
        long_count = len(long_indels)
        min_segment_bp = min(int(segment["aligned_bp"]) for segment in segments)
        min_segment_identity = min(segment_identity(segment) for segment in segments)
        segment_short_indel_bp = sum(int(segment["short_indel_bp"]) for segment in segments)
        same_operation = len({operation for operation, _ in long_indels}) <= 1
        max_gap = max(int(patch["ref_gap"]), int(patch["query_gap"]))

        single_balanced = single_pass(patch, 100, 0.95, 100, 0.15)
        multi_balanced = multi_pass(patch, 200, 0.98, 100, 0.15, "same_operation")
        near_global = (
            long_count == 0
            and float(patch["identity"]) >= 0.94
            and float(patch["max_short_error_density"]) <= 0.15
            and min(float(patch["left_endpoint_identity"]), float(patch["right_endpoint_identity"])) >= 0.85
        )
        if single_balanced:
            candidate_class = "clean_single_long_indel"
        elif multi_balanced:
            candidate_class = "clean_multi_long_indel"
        elif near_global:
            candidate_class = "near_global_identity_threshold"
        elif long_count == 0:
            candidate_class = "diffuse_or_no_long_indel"
        else:
            candidate_class = "noisy_or_weak_flank_long_indel"

        features.append(
            {
                "sample": gap["sample"],
                "left_chain_id": patch["left_chain_id"],
                "right_chain_id": patch["right_chain_id"],
                "ref_gap": patch["ref_gap"],
                "query_gap": patch["query_gap"],
                "max_gap": max_gap,
                "gap_difference": abs(int(patch["ref_gap"]) - int(patch["query_gap"])),
                "gap_bin": gap_bin(max_gap),
                "global_identity": patch["identity"],
                "long_indel_count": long_count,
                "long_indel_operations": ",".join(operation for operation, _ in long_indels) or ".",
                "long_indel_lengths": ",".join(str(length) for _, length in long_indels) or ".",
                "largest_long_indel": max((length for _, length in long_indels), default=0),
                "single_left_flank_bp": patch["rescue_left_bp"],
                "single_right_flank_bp": patch["rescue_right_bp"],
                "single_left_flank_identity": patch["rescue_left_identity"],
                "single_right_flank_identity": patch["rescue_right_identity"],
                "single_extra_indel_bp": patch["rescue_extra_indel_bp"],
                "multi_min_segment_bp": min_segment_bp,
                "multi_min_segment_identity": f"{min_segment_identity:.6f}",
                "multi_short_indel_bp": segment_short_indel_bp,
                "multi_same_operation": int(same_operation),
                "left_endpoint_identity": patch["left_endpoint_identity"],
                "right_endpoint_identity": patch["right_endpoint_identity"],
                "max_short_error_density": patch["max_short_error_density"],
                "low_quality_cigar": patch["low_quality_cigar"],
                "candidate_class": candidate_class,
                "cigar": patch["cigar"],
            }
        )
    write_tsv(args.output_dir / "low_identity_patch_features.tsv", features)

    bin_rows: list[dict[str, object]] = []
    classes = sorted({str(row["candidate_class"]) for row in features})
    for label in GAP_LABELS:
        selected = [row for row in features if row["gap_bin"] == label]
        bin_rows.append(
            {"gap_bin": label, "total": len(selected), **{
                candidate_class: sum(row["candidate_class"] == candidate_class for row in selected)
                for candidate_class in classes
            }}
        )
    write_tsv(args.output_dir / "low_identity_gap_length_bins.tsv", bin_rows)

    class_summary: list[dict[str, object]] = []
    for candidate_class in sorted({str(row["candidate_class"]) for row in features}):
        selected = [row for row in features if row["candidate_class"] == candidate_class]
        gap_lengths = [int(row["max_gap"]) for row in selected]
        class_summary.append({
            "candidate_class": candidate_class,
            "count": len(selected),
            "gap_min": min(gap_lengths),
            "gap_q25": quantile(gap_lengths, 0.25),
            "gap_median": quantile(gap_lengths, 0.50),
            "gap_q75": quantile(gap_lengths, 0.75),
            "gap_max": max(gap_lengths),
        })
    write_tsv(args.output_dir / "candidate_class_gap_summary.tsv", class_summary)

    failure_counts: Counter[str] = Counter()
    for row in source_rows:
        if int(row["rescue_long_indel_count"]) != 1:
            continue
        failures = []
        if int(row["rescue_left_bp"]) < 200:
            failures.append("left_flank_lt_200")
        if int(row["rescue_right_bp"]) < 200:
            failures.append("right_flank_lt_200")
        if float(row["rescue_left_identity"]) < 0.95:
            failures.append("left_flank_identity_lt_0.95")
        if float(row["rescue_right_identity"]) < 0.95:
            failures.append("right_flank_identity_lt_0.95")
        if int(row["rescue_extra_indel_bp"]) > 100:
            failures.append("extra_indel_gt_100")
        if float(row["max_short_error_density"]) > 0.15:
            failures.append("short_error_density_gt_0.15")
        failure_counts["+".join(failures) or "passes_current"] += 1
    write_tsv(
        args.output_dir / "single_long_indel_failure_reasons.tsv",
        [{"failure_combination": key, "count": value} for key, value in failure_counts.most_common()],
    )

    profiles: list[dict[str, object]] = []
    single_profiles = [
        ("single_current", 200, 0.95, 100, 0.15),
        ("single_flank_100", 100, 0.95, 100, 0.15),
        ("single_flank_50", 50, 0.95, 100, 0.15),
        ("single_relaxed", 50, 0.90, 2000, 0.25),
    ]
    for name, flank, identity, extra, density in single_profiles:
        profiles.append({
            "profile": name, "scope": "single_long_indel",
            "flank_or_segment_bp": flank, "identity": identity,
            "extra_or_short_indel_cap": extra, "max_short_error_density": density,
            "operation_mode": "single", "recoverable": sum(
                single_pass(row, flank, identity, extra, density) for row in source_rows
            ),
        })
    multi_profiles = [
        ("multi_current_insertions", 200, 1.00, 0, 0.00, "insertions"),
        ("multi_exact_same_operation", 200, 1.00, 0, 0.00, "same_operation"),
        ("multi_identity_0.99", 200, 0.99, 100, 0.15, "same_operation"),
        ("multi_identity_0.98", 200, 0.98, 100, 0.15, "same_operation"),
        ("multi_identity_0.95", 200, 0.95, 100, 0.15, "same_operation"),
        ("multi_segment_100", 100, 0.95, 250, 0.15, "same_operation"),
    ]
    for name, segment_bp, identity, short_cap, density, mode in multi_profiles:
        profiles.append({
            "profile": name, "scope": "multi_long_indel",
            "flank_or_segment_bp": segment_bp, "identity": identity,
            "extra_or_short_indel_cap": short_cap, "max_short_error_density": density,
            "operation_mode": mode, "recoverable": sum(
                multi_pass(row, segment_bp, identity, short_cap, density, mode)
                for row in source_rows
            ),
        })
    for threshold in (0.949, 0.945, 0.94, 0.93, 0.92, 0.90):
        profiles.append({
            "profile": f"global_identity_{threshold:.3f}", "scope": "no_long_indel",
            "flank_or_segment_bp": ".", "identity": threshold,
            "extra_or_short_indel_cap": ".", "max_short_error_density": 0.15,
            "operation_mode": "none", "recoverable": sum(
                int(row["rescue_long_indel_count"]) == 0
                and float(row["identity"]) >= threshold
                and float(row["max_short_error_density"]) <= 0.15
                and min(float(row["left_endpoint_identity"]), float(row["right_endpoint_identity"])) >= 0.85
                for row in source_rows
            ),
        })
    write_tsv(args.output_dir / "rescue_threshold_sweep.tsv", profiles)
    plot_distributions(args.output_dir / "low_identity_patch_distributions.png", features, profiles)

    print(f"low_identity={len(features)}")
    print("long_indel_count=" + str(Counter(int(row["long_indel_count"]) for row in features)))
    print("candidate_class=" + str(Counter(str(row["candidate_class"]) for row in features)))


if __name__ == "__main__":
    main()
