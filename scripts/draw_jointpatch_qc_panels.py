#!/usr/bin/env python3
"""Compose local truth-left QC panels and full-length SV overview PNGs."""

import argparse
import csv
import gzip
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def parse_info(info):
    result = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item, "true"
        result[key] = value
    return result


def infer_sample(path):
    name = Path(path).name
    for suffix in (".vcf.gz", ".vcf"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if "." in name:
        name = name.split(".")[0]
    return name


def read_vcf(path, kind):
    sample = infer_sample(path)
    records = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 9:
                    sample = fields[9]
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            info = parse_info(fields[7])
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(",")[0]
            svtype = info.get("SVTYPE")
            if not svtype:
                if len(ref) > len(alt):
                    svtype = "DEL"
                elif len(alt) > len(ref):
                    svtype = "INS"
                else:
                    svtype = "SNV"
            end = int(info.get("END", pos + max(1, len(ref)) - 1))
            try:
                svlen = int(info.get("SVLEN", len(alt) - len(ref)).split(",")[0])
            except ValueError:
                svlen = len(alt) - len(ref)
            records.append({
                "sample": sample,
                "chrom": fields[0],
                "pos": pos,
                "start0": max(0, pos - 1),
                "end0": max(pos, end),
                "id": fields[2] if fields[2] != "." else f"{kind}_{pos}_{svtype}",
                "svtype": svtype,
                "svlen": svlen,
                "kind": kind,
            })
    return records


def load_sequences(path):
    lengths = {}
    for row in read_tsv(path):
        lengths[row["sequence_id"]] = int(row["length"])
    return lengths


def read_paf(path):
    records = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for item in fields[12:]:
                parts = item.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            records.append({
                "qname": fields[0],
                "qlen": int(fields[1]),
                "qstart": int(fields[2]),
                "qend": int(fields[3]),
                "strand": fields[4],
                "tname": fields[5],
                "tlen": int(fields[6]),
                "tstart": int(fields[7]),
                "tend": int(fields[8]),
                "mode": tags.get("am", "."),
                "source": tags.get("bt", "."),
            })
    return records


def sv_color(svtype, kind="truth"):
    if kind == "fp":
        return (245, 130, 35)
    if svtype.startswith("DEL"):
        return (215, 55, 55)
    if svtype.startswith("INS"):
        return (55, 135, 210)
    if svtype.startswith("DUP"):
        return (210, 55, 170)
    if svtype.startswith("INV"):
        return (120, 70, 210)
    return (80, 150, 80)


def draw_truth_panel(row, width, height):
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()
    kind = row.get("kind", "truth")
    svtype = row["svtype"]
    color = sv_color(svtype, kind)
    title = "TRUTH SV" if kind == "truth" else "FP CALL"
    draw.text((24, 24), title, fill=(0, 0, 0), font=font)
    draw.text((24, 48), f"{row['sample']}  {row['sv_id']}", fill=(0, 0, 0), font=font)
    draw.text((24, 72), f"{svtype}  pos={int(row['pos']):,}", fill=color, font=font)
    draw.text((24, 96), f"end={int(row['end']):,}  svlen={int(row['svlen']):,}", fill=(0, 0, 0), font=font)

    y = 190
    x0 = 45
    x1 = width - 45
    draw.line([x0, y, x1, y], fill=(45, 45, 45), width=6)
    mid = (x0 + x1) // 2
    span = max(1, abs(int(row["end"]) - int(row["pos"])))
    half = max(8, min(90, int(math.log10(span + 10) * 18)))
    a = mid - half
    b = mid + half
    if kind == "fp":
        draw.rectangle([a, y - 38, b, y + 38], outline=color, width=4)
        draw.text((a, y - 62), "called variant", fill=color, font=font)
    elif svtype.startswith("INS"):
        draw.polygon([(mid, y - 55), (mid - 28, y + 25), (mid + 28, y + 25)], fill=color)
        draw.line([mid, y - 68, mid, y + 68], fill=color, width=2)
    elif svtype.startswith("DEL"):
        draw.rectangle([a, y - 22, b, y + 22], fill=(255, 225, 225), outline=color, width=3)
        draw.line([a, y, b, y], fill="white", width=8)
    elif svtype.startswith("DUP"):
        draw.rectangle([a, y - 28, b, y + 28], fill=(250, 220, 245), outline=color, width=3)
        draw.rectangle([a + 18, y - 46, b + 18, y + 10], outline=color, width=3)
    elif svtype.startswith("INV"):
        draw.rectangle([a, y - 28, b, y + 28], fill=(235, 225, 255), outline=color, width=3)
        draw.line([a + 10, y, b - 10, y], fill=color, width=3)
        draw.polygon([(a + 10, y), (a + 28, y - 10), (a + 28, y + 10)], fill=color)
        draw.polygon([(b - 10, y), (b - 28, y - 10), (b - 28, y + 10)], fill=color)
    else:
        draw.rectangle([a, y - 28, b, y + 28], fill=(225, 245, 225), outline=color, width=3)

    draw.text((24, y + 92), "Left: expected event", fill=(60, 60, 60), font=font)
    draw.text((24, y + 116), "Right: ValleScope2 local evidence", fill=(60, 60, 60), font=font)
    return image


def compose_local(summary_path, out_dir, panel_width):
    out_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    rows_out = []
    for row in read_tsv(summary_path):
        source = Path(row["png"])
        if not source.exists():
            continue
        right = Image.open(source).convert("RGB")
        left = draw_truth_panel(row, panel_width, right.height)
        canvas = Image.new("RGB", (left.width + right.width, right.height), "white")
        canvas.paste(left, (0, 0))
        canvas.paste(right, (left.width, 0))
        rel = source.relative_to(summary_path.parent)
        dest = out_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        canvas.save(dest)
        copied = dict(row)
        copied["png"] = str(dest)
        rows_out.append(copied)
        count += 1
    with open(out_dir / "summary.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows_out[0].keys()))
        writer.writeheader()
        writer.writerows(rows_out)
    return count


def ribbon_color(record):
    if record["source"] == "patched":
        return (40, 145, 75, 85)
    if record["source"] == "clean_patch":
        return (45, 120, 210, 75)
    if record["strand"] == "-":
        return (120, 70, 210, 70)
    return (90, 90, 90, 55)


def draw_full_length(paf_records, svs, lengths, out_png, sample, width, height):
    ref_name = "chr1_124500000-125500000"
    records = [r for r in paf_records if r["tname"] == ref_name and r["qname"] == sample]
    if not records:
        return False
    qlen = lengths.get(sample, max(r["qlen"] for r in records))
    tlen = lengths.get(ref_name, max(r["tlen"] for r in records))
    scale_len = max(qlen, tlen)
    image = Image.new("RGBA", (width, height), "white")
    overlay = Image.new("RGBA", image.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)
    base = ImageDraw.Draw(image)
    font = ImageFont.load_default()
    margin = 110
    right = width - 70
    yq = 150
    yt = height - 150
    track_h = 12

    def x(pos):
        return margin + (right - margin) * pos / max(1, scale_len)

    for r in records:
        q0, q1 = x(r["qstart"]), x(r["qend"])
        t0, t1 = x(r["tstart"]), x(r["tend"])
        if r["strand"] == "+":
            poly = [(q0, yq + track_h), (q1, yq + track_h),
                    (t1, yt - track_h), (t0, yt - track_h)]
        else:
            poly = [(q0, yq + track_h), (q1, yq + track_h),
                    (t0, yt - track_h), (t1, yt - track_h)]
        draw.polygon(poly, fill=ribbon_color(r))
    image = Image.alpha_composite(image, overlay)
    base = ImageDraw.Draw(image)
    base.rectangle([margin, yq - 5, x(qlen), yq + 5], fill=(40, 40, 40))
    base.rectangle([margin, yt - 5, x(tlen), yt + 5], fill=(40, 40, 40))
    sample_svs = [sv for sv in svs if sv["sample"] == sample]
    for sv in sample_svs:
        sx = x(sv["start0"])
        ex = x(max(sv["end0"], sv["start0"] + 1))
        if ex - sx < 3:
            ex = sx + 3
        color = sv_color(sv["svtype"], sv["kind"])
        base.rectangle([sx, yt - 42, ex, yt - 18], fill=color)
        base.line([sx, yt - 56, sx, yt + 56], fill=color, width=1)
    tick_step = 100000
    for tick in range(0, scale_len + tick_step, tick_step):
        tx = x(min(tick, scale_len))
        base.line([tx, yq + 22, tx, yq + 32], fill=(155, 155, 155))
        base.line([tx, yt - 32, tx, yt - 22], fill=(155, 155, 155))
        if tick % 200000 == 0:
            base.text((tx - 18, yq + 38), f"{tick/1000000:.1f}M", fill=(70, 70, 70), font=font)
    base.text((margin, 34), f"Full-length overview: {sample} vs reference", fill=(0, 0, 0), font=font)
    base.text((margin, 56), f"PAF ribbons={len(records)}  truth SV markers={len(sample_svs)}", fill=(0, 0, 0), font=font)
    base.text((margin, yq - 35), f"query: {sample} ({qlen:,} bp)", fill=(0, 0, 0), font=font)
    base.text((margin, yt + 25), f"reference: {ref_name} ({tlen:,} bp)", fill=(0, 0, 0), font=font)
    legend = [("DEL", sv_color("DEL")), ("INS", sv_color("INS")),
              ("DUP", sv_color("DUP")), ("INV", sv_color("INV")),
              ("patched/clean PAF", (40, 145, 75))]
    lx = margin
    for label, color in legend:
        base.rectangle([lx, height - 55, lx + 20, height - 42], fill=color)
        base.text((lx + 26, height - 57), label, fill=(0, 0, 0), font=font)
        lx += 145
    out_png.parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").save(out_png)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--local-summary", required=True)
    parser.add_argument("--local-out-dir", required=True)
    parser.add_argument("--panel-width", type=int, default=430)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--truth-vcf", action="append", required=True)
    parser.add_argument("--full-out-dir", required=True)
    parser.add_argument("--width", type=int, default=2200)
    parser.add_argument("--height", type=int, default=900)
    args = parser.parse_args()

    local_count = compose_local(Path(args.local_summary), Path(args.local_out_dir), args.panel_width)
    lengths = load_sequences(args.sequences)
    paf_records = read_paf(args.paf)
    svs = []
    for vcf in args.truth_vcf:
        svs.extend(read_vcf(vcf, "truth"))
    full_out = Path(args.full_out_dir)
    full_count = 0
    full_paths = []
    for sample in sorted({sv["sample"] for sv in svs}):
        path = full_out / f"{sample}.full_length.png"
        if draw_full_length(
            paf_records, svs, lengths, path, sample, args.width, args.height):
            full_count += 1
            full_paths.append(path)
    if full_paths:
        images = [Image.open(path).convert("RGB") for path in full_paths]
        combined = Image.new(
            "RGB",
            (max(image.width for image in images), sum(image.height for image in images)),
            "white")
        y = 0
        for image in images:
            combined.paste(image, (0, y))
            y += image.height
        combined.save(full_out / "all_samples.full_length.png")
    print(f"Wrote {local_count} local truth-left PNG files")
    print(f"Wrote {full_count} full-length PNG files")


if __name__ == "__main__":
    main()
