#!/usr/bin/env python3
"""Temporary helper: draw a simple two-track PAF ribbon plot as PNG using PIL."""

import argparse
import random
from collections import Counter, defaultdict

from PIL import Image, ImageDraw, ImageFont


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
                "mapq": int(fields[11]),
                "support": tags.get("sd", "."),
            })
    return records


def infer_pair(records):
    counts = Counter((r["qname"], r["tname"]) for r in records)
    if not counts:
        raise SystemExit("No PAF records.")
    return counts.most_common(1)[0][0]


def support_color(support):
    if support == "both":
        return (45, 120, 220, 72)
    if support == "forward_only":
        return (30, 170, 80, 60)
    if support == "reverse_only":
        return (240, 150, 40, 60)
    if support == "discordant":
        return (220, 40, 70, 70)
    return (120, 120, 120, 55)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paf")
    parser.add_argument("png")
    parser.add_argument("--qname", default=None)
    parser.add_argument("--tname", default=None)
    parser.add_argument("--max-ribbons", type=int, default=20000)
    parser.add_argument("--width", type=int, default=1800)
    parser.add_argument("--height", type=int, default=650)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    records = read_paf(args.paf)
    qname, tname = (args.qname, args.tname) if args.qname and args.tname else infer_pair(records)
    records = [r for r in records if r["qname"] == qname and r["tname"] == tname]
    if not records:
        raise SystemExit(f"No records for {qname} -> {tname}")

    if len(records) > args.max_ribbons:
        random.seed(1)
        records = random.sample(records, args.max_ribbons)

    qlen = max(r["qlen"] for r in records)
    tlen = max(r["tlen"] for r in records)

    image = Image.new("RGBA", (args.width, args.height), "white")
    overlay = Image.new("RGBA", image.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)
    base = ImageDraw.Draw(image)

    margin = 110
    right = args.width - 60
    yq = 170
    yt = 470
    track_h = 12

    def xq(pos):
        return margin + (right - margin) * pos / max(1, qlen)

    def xt(pos):
        return margin + (right - margin) * pos / max(1, tlen)

    for r in records:
        q0, q1 = xq(r["qstart"]), xq(r["qend"])
        t0, t1 = xt(r["tstart"]), xt(r["tend"])
        if r["strand"] == "+":
            poly = [(q0, yq + track_h), (q1, yq + track_h),
                    (t1, yt - track_h), (t0, yt - track_h)]
        else:
            poly = [(q0, yq + track_h), (q1, yq + track_h),
                    (t0, yt - track_h), (t1, yt - track_h)]
        draw.polygon(poly, fill=support_color(r["support"]))

    image = Image.alpha_composite(image, overlay)
    base = ImageDraw.Draw(image)
    font = ImageFont.load_default()

    base.rectangle([margin, yq - track_h // 2, right, yq + track_h // 2], fill=(40, 40, 40, 255))
    base.rectangle([margin, yt - track_h // 2, right, yt + track_h // 2], fill=(40, 40, 40, 255))
    base.text((margin, yq - 35), f"query: {qname} ({qlen:,} bp)", fill=(0, 0, 0, 255), font=font)
    base.text((margin, yt + 22), f"target: {tname} ({tlen:,} bp)", fill=(0, 0, 0, 255), font=font)
    title = args.title or f"{len(records):,} PAF anchor correspondences"
    base.text((margin, 35), title, fill=(0, 0, 0, 255), font=font)

    legend = [("both", support_color("both")), ("forward_only", support_color("forward_only")),
              ("reverse_only", support_color("reverse_only")), ("discordant", support_color("discordant"))]
    lx = margin
    for label, color in legend:
        base.rectangle([lx, args.height - 55, lx + 22, args.height - 40], fill=color)
        base.text((lx + 28, args.height - 57), label, fill=(0, 0, 0, 255), font=font)
        lx += 160

    image.convert("RGB").save(args.png)


if __name__ == "__main__":
    main()
