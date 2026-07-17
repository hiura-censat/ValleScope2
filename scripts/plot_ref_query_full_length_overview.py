#!/usr/bin/env python3
"""Draw sample-to-reference PAF ribbons and coverage summaries."""

import argparse
import re
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


CIGAR_RE = re.compile(r"(\d+)([=XMID])")


def read_paf(path, long_indel):
    records = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            tags = {}
            for value in fields[12:]:
                parts = value.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            record = {
                "qname": fields[0], "qlen": int(fields[1]), "qstart": int(fields[2]),
                "qend": int(fields[3]), "strand": fields[4], "tname": fields[5],
                "tlen": int(fields[6]), "tstart": int(fields[7]), "tend": int(fields[8]),
                "matches": int(fields[9]), "block": int(fields[10]), "tags": tags,
            }
            records.extend(split_record(record, long_indel))
    return records


def split_record(record, threshold):
    cigar = record["tags"].get("cg")
    if not cigar or not any(op in "ID" and int(length) >= threshold for length, op in CIGAR_RE.findall(cigar)):
        return [record]
    qpos = record["qstart"] if record["strand"] == "+" else record["qend"]
    tpos = record["tstart"]
    q0, t0 = qpos, tpos
    chunks = []

    def emit(q1, t1):
        if q1 == q0 or t1 == t0:
            return
        item = dict(record)
        item["qstart"], item["qend"] = sorted((q0, q1))
        item["tstart"], item["tend"] = sorted((t0, t1))
        item["block"] = max(item["qend"] - item["qstart"], item["tend"] - item["tstart"])
        chunks.append(item)

    for length_text, operation in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if operation in "ID" and length >= threshold:
            emit(qpos, tpos)
            if operation == "I":
                qpos += length if record["strand"] == "+" else -length
            else:
                tpos += length
            q0, t0 = qpos, tpos
            continue
        if operation in "=XM":
            qpos += length if record["strand"] == "+" else -length
            tpos += length
        elif operation == "I":
            qpos += length if record["strand"] == "+" else -length
        elif operation == "D":
            tpos += length
    emit(qpos, tpos)
    return chunks or [record]


def union_length(intervals):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged), merged


def largest_gap(merged, length):
    cursor = 0
    maximum = 0
    for start, end in merged:
        maximum = max(maximum, start - cursor)
        cursor = max(cursor, end)
    return max(maximum, length - cursor)


def gap_summary(merged, length):
    gaps = []
    cursor = 0
    for start, end in merged:
        if start > cursor:
            gaps.append(start - cursor)
        cursor = max(cursor, end)
    if cursor < length:
        gaps.append(length - cursor)
    return len(gaps), sum(gap >= 1_000 for gap in gaps), sum(
        gap >= 10_000 for gap in gaps)


def record_identity(record, method):
    if record["matches"] > 0 and record["block"] > 0:
        return min(1.0, record["matches"] / record["block"])
    if method == "centrolign" and "id" in record["tags"]:
        return min(1.0, float(record["tags"]["id"]) / 100.0)
    return 1.0


def draw_panel(records, sample, label, method, target_label, width, height):
    image = Image.new("RGBA", (width, height), "white")
    draw = ImageDraw.Draw(image, "RGBA")
    font = ImageFont.load_default()
    margin, right = 85, width - 25
    yq, yt, track = 76, height - 64, 8
    qlen = max(record["qlen"] for record in records)
    tlen = max(record["tlen"] for record in records)
    scale = max(qlen, tlen)
    x = lambda pos: margin + (right - margin) * pos / scale
    qcov, qmerged = union_length((r["qstart"], r["qend"]) for r in records)
    tcov, tmerged = union_length((r["tstart"], r["tend"]) for r in records)

    draw.rectangle((margin, yq - track / 2, x(qlen), yq + track / 2), fill=(205, 70, 70, 255))
    draw.rectangle((margin, yt - track / 2, x(tlen), yt + track / 2), fill=(205, 70, 70, 255))
    for start, end in qmerged:
        draw.rectangle((x(start), yq - track / 2, x(end), yq + track / 2), fill=(40, 40, 40, 255))
    for start, end in tmerged:
        draw.rectangle((x(start), yt - track / 2, x(end), yt + track / 2), fill=(40, 40, 40, 255))
    for record in records:
        q0, q1 = x(record["qstart"]), x(record["qend"])
        t0, t1 = x(record["tstart"]), x(record["tend"])
        identity = record_identity(record, method)
        alpha = int(35 + 105 * max(0.0, min(1.0, identity)))
        color = (38, 105, 170, alpha) if record["strand"] == "+" else (190, 55, 65, alpha)
        lower = (t1, t0) if record["strand"] == "+" else (t0, t1)
        draw.polygon(((q0, yq + track), (q1, yq + track), (lower[0], yt - track), (lower[1], yt - track)), fill=color)
    draw.text((10, 8), label, fill=(0, 0, 0, 255), font=font)
    draw.text((10, 27), f"blocks={len(records)}  query={qcov/qlen:.1%}  ref={tcov/tlen:.1%}", fill=(30, 30, 30, 255), font=font)
    draw.text((10, yq - 7), sample, fill=(0, 0, 0, 255), font=font)
    draw.text((10, yt - 7), target_label, fill=(0, 0, 0, 255), font=font)
    for tick in range(0, int(scale) + 500_000, 500_000):
        if tick > scale:
            break
        tx = x(tick)
        draw.line((tx, yt + 9, tx, yt + 14), fill=(100, 100, 100, 255))
        draw.text((tx - 10, yt + 16), f"{tick/1e6:g}", fill=(80, 80, 80, 255), font=font)
    qgap_count, qgap_1k, qgap_10k = gap_summary(qmerged, qlen)
    tgap_count, tgap_1k, tgap_10k = gap_summary(tmerged, tlen)
    return image.convert("RGB"), (
        qlen, qcov, largest_gap(qmerged, qlen), qgap_count, qgap_1k, qgap_10k,
        tlen, tcov, largest_gap(tmerged, tlen), tgap_count, tgap_1k, tgap_10k,
        len(records))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", action="append", type=Path, required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--long-indel-bp", type=int, default=10_000)
    parser.add_argument("--target-label", default="CHM13")
    parser.add_argument("--panel-width", type=int, default=900)
    parser.add_argument("--panel-height", type=int, default=260)
    args = parser.parse_args()
    if len(args.paf) != len(args.label):
        parser.error("--paf and --label counts must match")
    methods = [label.lower().split()[0] for label in args.label]
    datasets = [read_paf(path, args.long_indel_bp) for path in args.paf]
    samples = sorted(set.intersection(*[{r["qname"] for r in records} for records in datasets]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    combined_rows = []
    for sample in samples:
        panels = []
        for label, method, records in zip(args.label, methods, datasets):
            selected = [record for record in records if record["qname"] == sample]
            panel, stats = draw_panel(
                selected, sample, label, method, args.target_label,
                args.panel_width, args.panel_height)
            panels.append(panel)
            (qlen, qcov, qgap, qgap_count, qgap_1k, qgap_10k,
             tlen, tcov, tgap, tgap_count, tgap_1k, tgap_10k, blocks) = stats
            rows.append((
                sample, label, blocks, qlen, qcov, qcov / qlen, qgap,
                qgap_count, qgap_1k, qgap_10k, tlen, tcov, tcov / tlen, tgap,
                tgap_count, tgap_1k, tgap_10k))
        canvas = Image.new("RGB", (args.panel_width * len(panels), args.panel_height), "white")
        for index, panel in enumerate(panels):
            canvas.paste(panel, (index * args.panel_width, 0))
        canvas.save(args.output_dir / f"{sample}.full_length.png")
        combined_rows.append(canvas)
    combined = Image.new("RGB", (args.panel_width * len(args.paf), args.panel_height * len(combined_rows)), "white")
    for index, image in enumerate(combined_rows):
        combined.paste(image, (0, index * args.panel_height))
    combined.save(args.output_dir / "all_samples.full_length.png")
    with (args.output_dir / "coverage_summary.tsv").open("w") as handle:
        handle.write("sample\tmethod\tblocks\tquery_length\tquery_covered_bp\tquery_coverage\tquery_largest_uncovered_bp\tquery_gap_count\tquery_gap_count_ge_1kb\tquery_gap_count_ge_10kb\treference_length\treference_covered_bp\treference_coverage\treference_largest_uncovered_bp\treference_gap_count\treference_gap_count_ge_1kb\treference_gap_count_ge_10kb\n")
        for row in rows:
            handle.write("\t".join(map(str, row[:5])) + f"\t{row[5]:.6f}\t" +
                         "\t".join(map(str, row[6:12])) + f"\t{row[12]:.6f}\t" +
                         "\t".join(map(str, row[13:])) + "\n")


if __name__ == "__main__":
    main()
