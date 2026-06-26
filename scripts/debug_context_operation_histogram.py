#!/usr/bin/env python3

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_alignment(path):
    by_rank = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        has_rank = "rank" in reader.fieldnames
        for row in reader:
            rank = int(row["rank"]) if has_rank else 0
            by_rank[rank].append(row)
    for rows in by_rank.values():
        rows.sort(key=lambda row: int(row.get("alignment_index", row.get("index", 0))))
    return by_rank


def operation_histogram(by_rank):
    counts = defaultdict(Counter)
    query_tokens = {}
    insertion_before_first = Counter()

    for rank, rows in by_rank.items():
        query_index = -1
        for row in rows:
            operation = row["operation"]
            operation = {
                "sub": "substitution",
                "del": "deletion",
                "ins": "insertion",
            }.get(operation, operation)
            query_token = row.get("query_token", row.get("context1_token", "."))

            if operation == "insertion":
                if query_index < 0:
                    insertion_before_first[rank] += 1
                else:
                    counts[query_index]["insertion_after"] += 1
                continue

            query_index += 1
            query_tokens.setdefault(query_index, query_token)
            if operation == "match":
                counts[query_index]["match"] += 1
            elif operation == "substitution":
                counts[query_index]["substitution"] += 1
            elif operation == "deletion":
                counts[query_index]["deletion"] += 1
            else:
                counts[query_index][operation] += 1

    return counts, query_tokens, sum(insertion_before_first.values())


def write_tsv(path, counts, query_tokens):
    fieldnames = [
        "query_index",
        "query_token",
        "match",
        "substitution",
        "deletion",
        "insertion_after",
        "edited_total",
        "aligned_total",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index in sorted(query_tokens):
            row_counts = counts[index]
            edited_total = (
                row_counts["substitution"]
                + row_counts["deletion"]
                + row_counts["insertion_after"]
            )
            aligned_total = (
                row_counts["match"]
                + row_counts["substitution"]
                + row_counts["deletion"]
            )
            writer.writerow(
                {
                    "query_index": index,
                    "query_token": query_tokens[index],
                    "match": row_counts["match"],
                    "substitution": row_counts["substitution"],
                    "deletion": row_counts["deletion"],
                    "insertion_after": row_counts["insertion_after"],
                    "edited_total": edited_total,
                    "aligned_total": aligned_total,
                }
            )


def svg_escape(text):
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_svg(path, counts, query_tokens):
    width = 1300
    height = 420
    margin_left = 55
    margin_right = 20
    margin_top = 35
    margin_bottom = 80
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    n = max(query_tokens) + 1 if query_tokens else 0
    max_count = max(
        (
            counts[i]["substitution"]
            + counts[i]["deletion"]
            + counts[i]["insertion_after"]
            for i in query_tokens
        ),
        default=1,
    ) or 1
    bar_width = plot_width / max(n, 1)
    colors = {
        "substitution": "#d95f02",
        "deletion": "#7570b3",
        "insertion_after": "#1b9e77",
    }
    labels = {
        "substitution": "substitution",
        "deletion": "deletion",
        "insertion_after": "insertion after",
    }

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin_left}" y="22" font-family="sans-serif" font-size="16">Edit operation histogram by query token index</text>',
        f'<line x1="{margin_left}" y1="{margin_top + plot_height}" x2="{width - margin_right}" y2="{margin_top + plot_height}" stroke="black"/>',
        f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_height}" stroke="black"/>',
    ]
    for tick in range(0, max_count + 1, max(1, max_count // 5 or 1)):
        y = margin_top + plot_height - tick / max_count * plot_height
        lines.append(f'<line x1="{margin_left - 4}" y1="{y:.1f}" x2="{margin_left}" y2="{y:.1f}" stroke="black"/>')
        lines.append(f'<text x="{margin_left - 8}" y="{y + 4:.1f}" text-anchor="end" font-family="sans-serif" font-size="10">{tick}</text>')

    for index in sorted(query_tokens):
        x = margin_left + index * bar_width
        y_bottom = margin_top + plot_height
        stacked = 0
        for key in ["substitution", "deletion", "insertion_after"]:
            value = counts[index][key]
            if value == 0:
                continue
            h = value / max_count * plot_height
            y = y_bottom - stacked - h
            lines.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(bar_width - 1, 1):.2f}" height="{h:.2f}" fill="{colors[key]}">'
                f'<title>index {index} {query_tokens[index]} {labels[key]}={value}</title></rect>'
            )
            stacked += h
        if index % 10 == 0:
            tx = x + bar_width / 2
            lines.append(f'<text x="{tx:.1f}" y="{height - 55}" text-anchor="middle" font-family="sans-serif" font-size="10">{index}</text>')

    legend_x = margin_left
    legend_y = height - 35
    for offset, key in enumerate(["substitution", "deletion", "insertion_after"]):
        x = legend_x + offset * 170
        lines.append(f'<rect x="{x}" y="{legend_y - 10}" width="14" height="14" fill="{colors[key]}"/>')
        lines.append(f'<text x="{x + 20}" y="{legend_y + 2}" font-family="sans-serif" font-size="12">{svg_escape(labels[key])}</text>')
    lines.append(f'<text x="{width / 2:.1f}" y="{height - 15}" text-anchor="middle" font-family="sans-serif" font-size="12">query token index</text>')
    lines.append("</svg>")
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Build a per-query-token operation histogram from context alignment TSV."
    )
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-svg")
    args = parser.parse_args()

    by_rank = read_alignment(args.alignment)
    counts, query_tokens, insertion_before_first = operation_histogram(by_rank)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_tsv, counts, query_tokens)
    if args.output_svg:
        Path(args.output_svg).parent.mkdir(parents=True, exist_ok=True)
        write_svg(args.output_svg, counts, query_tokens)
    print("ranks", len(by_rank))
    print("query_tokens", len(query_tokens))
    print("insertion_before_first", insertion_before_first)
    print("output_tsv", args.output_tsv)
    if args.output_svg:
        print("output_svg", args.output_svg)


if __name__ == "__main__":
    main()
