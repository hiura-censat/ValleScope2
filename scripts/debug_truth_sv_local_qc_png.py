#!/usr/bin/env python3
"""Draw local ValleScope2 evidence around each truth SV.

This is a temporary QC helper.  It overlays intermediate ValleScope2 results
around truth VCF records:

  - assignment correspondences
  - first-pass chain anchors/bundles
  - refined chain anchors/bundles
  - final bundle alignment PAF

The plot uses a local reference window and an inferred local query window for
the corresponding sample.
"""

import argparse
import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REFERENCE_SAMPLE = "reference"


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


def read_truth_vcf(path):
    records = []
    sample = None
    with open(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                sample = fields[9] if len(fields) > 9 else Path(path).stem.split(".")[0]
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            info = parse_info(fields[7])
            pos = int(fields[1])
            svtype = info.get("SVTYPE", "NA")
            end = int(info.get("END", pos))
            svlen = int(info.get("SVLEN", "0").split(",")[0])
            if end < pos:
                end = pos
            records.append({
                "sample": sample,
                "chrom": fields[0],
                "pos": pos,
                "start0": max(0, pos - 1),
                "end0": max(pos, end),
                "id": fields[2],
                "svtype": svtype,
                "svlen": svlen,
                "info": info,
            })
    return records


def load_sequences(path):
    lengths = {}
    samples = {}
    for row in read_tsv(path):
        lengths[row["sequence_id"]] = int(row["length"])
        samples[row["sample"]] = row["sequence_id"]
    return lengths, samples


def load_anchors(path):
    anchors = {}
    for row in read_tsv(path):
        start = int(row["start"])
        end = int(row["end"])
        anchors[row["anchor_id"]] = {
            "seq": row["sequence_id"],
            "start": start,
            "end": end,
            "center": (start + end) / 2.0,
        }
    return anchors


def load_correspondences(path, anchors):
    rows = []
    for row in read_tsv(path):
        a = anchors.get(row["anchor_a"])
        b = anchors.get(row["anchor_b"])
        if not a or not b:
            continue
        rows.append({
            "sample_a": row["sample_a"],
            "sample_b": row["sample_b"],
            "anchor_a": row["anchor_a"],
            "anchor_b": row["anchor_b"],
            "a": a,
            "b": b,
            "support": row["support_direction"],
            "score": max(parse_float(row.get("forward_score", ".")),
                         parse_float(row.get("reverse_score", "."))),
        })
    return rows


def parse_float(value):
    try:
        if value == ".":
            return float("-inf")
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def load_chain_anchors(path, anchors):
    rows = []
    for row in read_tsv(path):
        a = anchors.get(row["anchor_a"])
        b = anchors.get(row["anchor_b"])
        if not a or not b:
            continue
        rows.append({
            "chain_id": row["chain_id"],
            "sample_a": row["sample_a"],
            "sample_b": row["sample_b"],
            "a": a,
            "b": b,
            "strand": row.get("assign_strand", "+"),
            "score": parse_float(row.get("candidate_score", ".")),
            "support": row.get("support_direction", "."),
        })
    return rows


def load_bundles(path, refined=False):
    rows = []
    for row in read_tsv(path):
        rows.append({
            "chain_id": row["chain_id"],
            "sample_a": row["sample_a"],
            "sample_b": row["sample_b"],
            "sequence_a": row["sequence_a"],
            "sequence_b": row["sequence_b"],
            "strand": row.get("assign_strand", "+"),
            "ref_start": int(row["ref_start"]),
            "ref_end": int(row["ref_end"]),
            "query_start": int(row["query_start"]),
            "query_end": int(row["query_end"]),
            "score": parse_float(row.get("chain_score", ".")),
            "kind": row.get("refinement_type", "first_pass" if not refined else "refined"),
        })
    return rows


def read_paf(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for tag in fields[12:]:
                parts = tag.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            rows.append({
                "qname": fields[0],
                "qlen": int(fields[1]),
                "qstart": int(fields[2]),
                "qend": int(fields[3]),
                "strand": fields[4],
                "tname": fields[5],
                "tlen": int(fields[6]),
                "tstart": int(fields[7]),
                "tend": int(fields[8]),
                "matches": int(fields[9]),
                "block": int(fields[10]),
                "mapq": int(fields[11]),
                "chain_id": tags.get("ch", "."),
                "align_mode": tags.get("am", "."),
                "truth_variant": tags.get("tv", "."),
            })
    return rows


def overlap(a0, a1, b0, b1):
    return max(a0, b0) < min(a1, b1)


def clip_interval_to_ref(item, ref0, ref1):
    t0 = item["ref_start"] if "ref_start" in item else item["tstart"]
    t1 = item["ref_end"] if "ref_end" in item else item["tend"]
    q0 = item["query_start"] if "query_start" in item else item["qstart"]
    q1 = item["query_end"] if "query_end" in item else item["qend"]
    strand = item.get("strand", "+")
    if not overlap(t0, t1, ref0, ref1):
        return None
    ct0 = max(t0, ref0)
    ct1 = min(t1, ref1)
    span = max(1, t1 - t0)
    f0 = (ct0 - t0) / span
    f1 = (ct1 - t0) / span
    if strand == "-":
        cq0 = q1 - (q1 - q0) * f1
        cq1 = q1 - (q1 - q0) * f0
    else:
        cq0 = q0 + (q1 - q0) * f0
        cq1 = q0 + (q1 - q0) * f1
    if cq0 > cq1:
        cq0, cq1 = cq1, cq0
    return ct0, ct1, cq0, cq1


def pair_anchor(row, sample):
    a, b = row["a"], row["b"]
    if a["seq"] == "chr1_124500000-125500000" and b["seq"] == sample:
        return a, b
    if b["seq"] == "chr1_124500000-125500000" and a["seq"] == sample:
        return b, a
    return None


def choose_query_window(values, qlen, pad):
    if not values:
        return 0, min(qlen, max(1, pad * 2))
    q0 = max(0, int(min(values) - pad))
    q1 = min(qlen, int(max(values) + pad))
    if q1 <= q0:
        q1 = min(qlen, q0 + pad * 2)
    return q0, q1


def safe_name(text):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in text)


def draw_sv_plot(sv, out_png, lengths, evidence, window, query_pad, width, height):
    ref_seq = sv["chrom"]
    sample = sv["sample"]
    ref_len = lengths.get(ref_seq, sv["end0"] + window)
    q_len = lengths.get(sample, 1)
    ref0 = max(0, sv["start0"] - window)
    ref1 = min(ref_len, sv["end0"] + window)

    assignments = []
    for row in evidence["correspondences"]:
        pair = pair_anchor(row, sample)
        if not pair:
            continue
        ref_a, q_a = pair
        if ref0 <= ref_a["center"] <= ref1:
            assignments.append((ref_a, q_a, row))

    chain_anchors = []
    for row in evidence["chain_anchors"]:
        pair = pair_anchor(row, sample)
        if not pair:
            continue
        ref_a, q_a = pair
        if ref0 <= ref_a["center"] <= ref1:
            chain_anchors.append((ref_a, q_a, row))

    refined_anchors = []
    for row in evidence["refined_chain_anchors"]:
        pair = pair_anchor(row, sample)
        if not pair:
            continue
        ref_a, q_a = pair
        if ref0 <= ref_a["center"] <= ref1:
            refined_anchors.append((ref_a, q_a, row))

    bundles = [
        b for b in evidence["chains"]
        if b["sequence_a"] == ref_seq and b["sequence_b"] == sample
        and overlap(b["ref_start"], b["ref_end"], ref0, ref1)
    ]
    refined = [
        b for b in evidence["refined_chains"]
        if b["sequence_a"] == ref_seq and b["sequence_b"] == sample
        and overlap(b["ref_start"], b["ref_end"], ref0, ref1)
    ]
    pafs = [
        p for p in evidence["paf"]
        if p["tname"] == ref_seq and p["qname"] == sample
        and overlap(p["tstart"], p["tend"], ref0, ref1)
    ]

    query_values_assignment = []
    query_values_chain = []
    query_values_bundle = []
    query_values_paf = []
    for _, q_a, _ in assignments:
        query_values_assignment.extend([q_a["start"], q_a["end"]])
    for _, q_a, _ in chain_anchors + refined_anchors:
        query_values_chain.extend([q_a["start"], q_a["end"]])
    for collection in (bundles, refined):
        for item in collection:
            clipped = clip_interval_to_ref(item, ref0, ref1)
            if clipped:
                query_values_bundle.extend([clipped[2], clipped[3]])
    for item in pafs:
        clipped = clip_interval_to_ref(item, ref0, ref1)
        if clipped:
            query_values_paf.extend([clipped[2], clipped[3]])
    # Prefer final PAF/bundle coordinates for a readable local view.  Raw
    # assignment anchors in repeats can map to many remote copies and otherwise
    # stretch the query axis to nearly the entire chromosome.
    query_values = (
        query_values_paf or query_values_bundle or query_values_chain
        or query_values_assignment
    )
    q0, q1 = choose_query_window(query_values, q_len, query_pad)

    image = Image.new("RGBA", (width, height), "white")
    overlay = Image.new("RGBA", image.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)
    base = ImageDraw.Draw(image)
    font = ImageFont.load_default()

    left = 115
    right = width - 55
    y_ref = 95
    y_q = height - 95
    track_h = 10
    ref_span = max(1, ref1 - ref0)
    query_span = max(1, q1 - q0)
    common_span = max(ref_span, query_span)
    bp_scale = (right - left) / common_span
    ref_right = left + ref_span * bp_scale
    query_right = left + query_span * bp_scale

    def xr(pos):
        return left + (pos - ref0) * bp_scale

    def xq(pos):
        return left + (pos - q0) * bp_scale

    def in_q(pos):
        return q0 <= pos <= q1

    def draw_ribbon(item, fill, outline=None, width_line=1):
        clipped = clip_interval_to_ref(item, ref0, ref1)
        if not clipped:
            return
        rt0, rt1, qq0, qq1 = clipped
        if qq1 < q0 or qq0 > q1:
            return
        qq0 = max(q0, min(q1, qq0))
        qq1 = max(q0, min(q1, qq1))
        poly = [(xq(qq0), y_q - track_h),
                (xq(qq1), y_q - track_h),
                (xr(rt1), y_ref + track_h),
                (xr(rt0), y_ref + track_h)]
        draw.polygon(poly, fill=fill)
        if outline:
            draw.line(poly + [poly[0]], fill=outline, width=width_line)

    for item in bundles:
        draw_ribbon(item, (255, 150, 20, 28), (205, 100, 0, 170), 1)
    for item in refined:
        draw_ribbon(item, (150, 70, 230, 34), (110, 40, 190, 190), 2)
    for item in pafs:
        truth_variant = item.get("truth_variant", ".")
        if truth_variant == "DUP_original":
            draw_ribbon(item, (0, 140, 220, 18), (0, 95, 180, 230), 3)
        elif truth_variant == "DUP_extra_copy":
            draw_ribbon(item, (230, 30, 180, 18), (190, 0, 145, 230), 3)
        elif truth_variant == "INV":
            draw_ribbon(item, (125, 55, 210, 28), (90, 35, 170, 220), 3)
        elif truth_variant.endswith("_source"):
            draw_ribbon(item, (245, 145, 20, 25), (200, 95, 0, 220), 3)
        elif truth_variant == "flank":
            draw_ribbon(item, (90, 90, 90, 36), (120, 120, 120, 120), 1)
        else:
            color = (
                (90, 90, 90, 85)
                if item["align_mode"] != "anchor_guided"
                else (40, 120, 60, 78)
            )
            draw_ribbon(item, color)

    def draw_anchor_line(ref_a, q_a, color, radius=2):
        if not in_q(q_a["center"]):
            return
        x0 = xr(ref_a["center"])
        x1 = xq(q_a["center"])
        draw.line([x0, y_ref + 18, x1, y_q - 18], fill=color, width=1)
        draw.ellipse([x0 - radius, y_ref - radius, x0 + radius, y_ref + radius], fill=color)
        draw.ellipse([x1 - radius, y_q - radius, x1 + radius, y_q + radius], fill=color)

    max_anchor_lines = 600
    for ref_a, q_a, _ in assignments[:max_anchor_lines]:
        draw_anchor_line(ref_a, q_a, (40, 100, 220, 75), 1)
    for ref_a, q_a, _ in chain_anchors[:max_anchor_lines]:
        draw_anchor_line(ref_a, q_a, (230, 125, 10, 120), 2)
    for ref_a, q_a, _ in refined_anchors[:max_anchor_lines]:
        draw_anchor_line(ref_a, q_a, (120, 40, 210, 150), 2)

    image = Image.alpha_composite(image, overlay)
    base = ImageDraw.Draw(image)

    base.rectangle([left, y_ref - track_h // 2, ref_right, y_ref + track_h // 2],
                   fill=(45, 45, 45, 255))
    base.rectangle([left, y_q - track_h // 2, query_right, y_q + track_h // 2],
                   fill=(45, 45, 45, 255))

    sv0 = max(ref0, sv["start0"])
    sv1 = min(ref1, max(sv["end0"], sv["start0"] + 1))
    sx0, sx1 = xr(sv0), xr(sv1)
    if sx1 - sx0 < 3:
        sx1 = sx0 + 3
    base.rectangle([sx0, y_ref - 28, sx1, y_ref + 28],
                   fill=(230, 40, 50, 90), outline=(180, 0, 0, 255))
    base.line([sx0, y_ref - 38, sx0, y_q + 38], fill=(230, 40, 50, 120), width=1)

    for y, start, end, mapper, label in [
        (y_ref, ref0, ref1, xr, f"ref: {ref_seq}"),
        (y_q, q0, q1, xq, f"query: {sample}"),
    ]:
        span = end - start
        raw_step = max(1, common_span // 5)
        mag = 10 ** int(math.floor(math.log10(raw_step))) if raw_step else 1
        step = max(mag, round(raw_step / mag) * mag)
        tick = (start // step) * step
        while tick <= end + step:
            if start <= tick <= end:
                x = mapper(tick)
                base.line([x, y - 18, x, y + 18], fill=(150, 150, 150, 255), width=1)
                base.text((x - 22, y + 25 if y == y_ref else y - 42),
                          f"{tick:,}", fill=(70, 70, 70, 255), font=font)
            tick += step
        base.text((15, y - 8), label, fill=(0, 0, 0, 255), font=font)

    title = (
        f"{sv['id']}  {sv['svtype']}  pos={sv['pos']:,} "
        f"end={sv['end0']:,} svlen={sv['svlen']:,}"
    )
    base.text((left, 20), title, fill=(0, 0, 0, 255), font=font)
    counts = (
        f"assignment={len(assignments)}  chain_anchor={len(chain_anchors)}  "
        f"refined_anchor={len(refined_anchors)}  chain={len(bundles)}  "
        f"refined={len(refined)}  final_paf={len(pafs)}"
    )
    base.text((left, 40), counts, fill=(0, 0, 0, 255), font=font)
    scale_note = (
        f"common bp scale; ref window={ref_span:,} bp, "
        f"query window={query_span:,} bp"
    )
    base.text((left, 60), scale_note, fill=(80, 80, 80, 255), font=font)

    legend = [
        ("truth SV", (230, 40, 50, 130)),
        ("assignment", (40, 100, 220, 130)),
        ("chain", (230, 125, 10, 160)),
        ("refined", (120, 40, 210, 170)),
        ("final PAF", (40, 120, 60, 150)),
        ("DUP copy", (230, 30, 180, 130)),
    ]
    lx = left
    for label, color in legend:
        base.rectangle([lx, height - 35, lx + 20, height - 22], fill=color)
        base.text((lx + 25, height - 36), label, fill=(0, 0, 0, 255), font=font)
        lx += 145

    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").save(out_png)

    all_sample_paf = [
        p for p in evidence["paf"]
        if p["tname"] == ref_seq and p["qname"] == sample
    ]
    sv_start = sv["start0"]
    sv_end = max(sv["end0"], sv["start0"] + 1)
    spanning_paf = [
        p for p in all_sample_paf
        if p["tstart"] <= sv_start and p["tend"] >= sv_end
    ]
    left_paf = [
        p for p in all_sample_paf
        if p["tend"] <= sv_start and sv_start - p["tend"] <= window
    ]
    right_paf = [
        p for p in all_sample_paf
        if p["tstart"] >= sv_end and p["tstart"] - sv_end <= window
    ]
    if not assignments:
        status = "assignment_missing"
    elif not chain_anchors and not refined_anchors:
        status = "chaining_missing"
    elif not pafs:
        status = "base_alignment_missing"
    elif spanning_paf:
        status = "paf_spans_truth_interval"
    elif left_paf and right_paf:
        status = "paf_flanks_truth_interval"
    else:
        status = "paf_nearby_only"

    return {
        "sample": sample,
        "sv_id": sv["id"],
        "svtype": sv["svtype"],
        "pos": sv["pos"],
        "end": sv["end0"],
        "svlen": sv["svlen"],
        "ref_window_start": ref0,
        "ref_window_end": ref1,
        "query_window_start": q0,
        "query_window_end": q1,
        "assignment_count": len(assignments),
        "chain_anchor_count": len(chain_anchors),
        "refined_anchor_count": len(refined_anchors),
        "chain_bundle_count": len(bundles),
        "refined_bundle_count": len(refined),
        "final_paf_count": len(pafs),
        "paf_spanning_count": len(spanning_paf),
        "paf_left_flank_count": len(left_paf),
        "paf_right_flank_count": len(right_paf),
        "status": status,
        "png": str(out_png),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug-dir", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--truth-vcf", action="append", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--window-bp", type=int, default=10000)
    parser.add_argument("--query-pad-bp", type=int, default=10000)
    parser.add_argument("--width", type=int, default=1500)
    parser.add_argument("--height", type=int, default=760)
    parser.add_argument("--max-svs", type=int, default=None)
    parser.add_argument("--paf-only", action="store_true",
                        help="draw only the supplied PAF, without intermediate ValleScope2 evidence")
    args = parser.parse_args()

    debug_dir = Path(args.debug_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    lengths, _ = load_sequences(debug_dir / "sequences.tsv")
    if args.paf_only:
        evidence = {
            "correspondences": [],
            "chain_anchors": [],
            "refined_chain_anchors": [],
            "chains": [],
            "refined_chains": [],
            "paf": read_paf(args.paf),
        }
    else:
        anchors = load_anchors(debug_dir / "grouped_anchors.tsv")
        evidence = {
            "correspondences": load_correspondences(debug_dir / "anchor_correspondences.tsv", anchors),
            "chain_anchors": load_chain_anchors(debug_dir / "chain_anchors.tsv", anchors),
            "refined_chain_anchors": load_chain_anchors(debug_dir / "refined_chain_anchors.tsv", anchors),
            "chains": load_bundles(debug_dir / "chains.tsv"),
            "refined_chains": load_bundles(debug_dir / "refined_chains.tsv", refined=True),
            "paf": read_paf(args.paf),
        }

    svs = []
    for vcf in args.truth_vcf:
        svs.extend(read_truth_vcf(vcf))
    svs.sort(key=lambda r: (r["sample"], r["pos"], r["id"]))
    if args.max_svs is not None:
        svs = svs[:args.max_svs]

    summaries = []
    for sv in svs:
        sample_dir = out_dir / sv["sample"]
        name = f"{sv['pos']:07d}_{safe_name(sv['id'])}_{sv['svtype']}.png"
        png = sample_dir / name
        summaries.append(draw_sv_plot(
            sv, png, lengths, evidence, args.window_bp, args.query_pad_bp,
            args.width, args.height))

    summary_path = out_dir / "summary.tsv"
    with open(summary_path, "w", newline="") as handle:
        fieldnames = [
            "sample", "sv_id", "svtype", "pos", "end", "svlen",
            "ref_window_start", "ref_window_end",
            "query_window_start", "query_window_end",
            "assignment_count", "chain_anchor_count", "refined_anchor_count",
            "chain_bundle_count", "refined_bundle_count", "final_paf_count",
            "paf_spanning_count", "paf_left_flank_count", "paf_right_flank_count",
            "status", "png",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)

    print(f"Wrote {len(summaries)} PNG files")
    print(summary_path)


if __name__ == "__main__":
    main()
