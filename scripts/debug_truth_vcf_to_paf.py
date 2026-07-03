#!/usr/bin/env python3
"""Temporary helper: convert simulated truth VCFs into PAF-like truth blocks.

This is not intended as a general VCF-to-PAF converter.  It creates a compact
visual truth model for the current simulation:

  - unchanged reference segments are emitted as colinear '+' PAF records
  - DEL advances reference only
  - INS advances query only
  - DUP emits two query copies against the same reference interval
  - INV emits one '-' record for the inverted interval
  - records with INFO/SOURCE emit the inserted query segment against SOURCE

The output is useful for drawing ribbon plots next to ValleScope2 output.
"""

import argparse
import re
from pathlib import Path


def read_fasta_lengths(path):
    lengths = {}
    name = None
    length = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
    return lengths


def parse_info(text):
    info = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item, "true"
        info[key] = value
    return info


def parse_source(value):
    # Example: chr1_124500000-125500000:18135-19353
    match = re.match(r"(.+):(\d+)-(\d+)$", value)
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def parse_target(value):
    # Example: chr1_124500000-125500000:322210
    match = re.match(r"(.+):(\d+)$", value)
    if not match:
        return None
    return match.group(1), int(match.group(2))


def read_vcf(path):
    sample = None
    records = []
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
            pos1 = int(fields[1])
            svtype = info.get("SVTYPE", "NA")
            svlen = int(info.get("SVLEN", "0").split(",")[0])
            end1 = int(info.get("END", pos1))
            if svtype == "INS":
                ref_start = pos1 - 1
                ref_end = ref_start
                query_delta = abs(svlen) if svlen else max(0, len(fields[4]) - len(fields[3]))
            else:
                ref_start = pos1 - 1
                ref_end = max(ref_start, end1)
                query_delta = abs(svlen)
            if info.get("SOURCE") and info.get("TARGET"):
                target = parse_target(info["TARGET"])
                if target:
                    fields[0] = target[0]
                    ref_start = target[1]
                    ref_end = target[1]
            records.append({
                "chrom": fields[0],
                "pos1": pos1,
                "id": fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "info": info,
                "svtype": svtype,
                "svlen": svlen,
                "ref_start": ref_start,
                "ref_end": ref_end,
                "query_delta": query_delta,
            })
    records.sort(key=lambda r: (r["ref_start"], r["ref_end"], r["id"]))
    return sample, records


def emit_paf(out, qname, qlen, qstart, qend, strand,
             tname, tlen, tstart, tend, tags):
    if qend <= qstart or tend <= tstart:
        return
    block = max(qend - qstart, tend - tstart)
    match = min(qend - qstart, tend - tstart)
    fields = [
        qname, qlen, qstart, qend, strand,
        tname, tlen, tstart, tend,
        match, block, 60,
    ]
    tag_fields = [f"{k}:{typ}:{v}" for k, typ, v in tags]
    out.write("\t".join(map(str, fields + tag_fields)) + "\n")


def build_truth_paf(vcf, ref_lengths, query_lengths, out):
    sample, records = read_vcf(vcf)
    if not records:
        return
    chrom = records[0]["chrom"]
    tlen = ref_lengths[chrom]
    qlen = query_lengths[sample]
    ref_cursor = 0
    query_cursor = 0

    for rec in records:
        if rec["ref_start"] < ref_cursor:
            # Overlapping truth records are rare in this simulation.  Keep the
            # query model monotonic by ignoring the overlapping reference span.
            rec_start = ref_cursor
        else:
            rec_start = rec["ref_start"]

        flank_len = max(0, rec_start - ref_cursor)
        emit_paf(
            out, sample, qlen, query_cursor, query_cursor + flank_len, "+",
            chrom, tlen, ref_cursor, ref_cursor + flank_len,
            [("tv", "Z", "flank"), ("sv", "Z", rec["id"])])
        query_cursor += flank_len
        ref_cursor += flank_len

        svtype = rec["svtype"]
        ref_len = max(0, rec["ref_end"] - rec["ref_start"])
        ins_len = max(0, rec["query_delta"])

        if svtype == "DEL":
            ref_cursor = max(ref_cursor, rec["ref_end"])
        elif svtype == "INV":
            length = ref_len or ins_len
            emit_paf(
                out, sample, qlen, query_cursor, query_cursor + length, "-",
                chrom, tlen, rec["ref_start"], rec["ref_start"] + length,
                [("tv", "Z", "INV"), ("sv", "Z", rec["id"])])
            query_cursor += length
            ref_cursor = max(ref_cursor, rec["ref_start"] + length)
        elif svtype == "DUP":
            length = ref_len or ins_len
            emit_paf(
                out, sample, qlen, query_cursor, query_cursor + length, "+",
                chrom, tlen, rec["ref_start"], rec["ref_start"] + length,
                [("tv", "Z", "DUP_original"), ("sv", "Z", rec["id"])])
            emit_paf(
                out, sample, qlen, query_cursor + length, query_cursor + 2 * length, "+",
                chrom, tlen, rec["ref_start"], rec["ref_start"] + length,
                [("tv", "Z", "DUP_extra_copy"), ("sv", "Z", rec["id"])])
            query_cursor += 2 * length
            ref_cursor = max(ref_cursor, rec["ref_start"] + length)
        else:
            source = parse_source(rec["info"].get("SOURCE", ""))
            if source and ins_len:
                source_name, source_start, source_end = source
                source_len = max(0, source_end - source_start)
                length = min(ins_len, source_len) if source_len else ins_len
                emit_paf(
                    out, sample, qlen, query_cursor, query_cursor + length, "+",
                    source_name, ref_lengths.get(source_name, tlen),
                    source_start, source_start + length,
                    [("tv", "Z", svtype + "_source"), ("sv", "Z", rec["id"])])
            query_cursor += ins_len
            ref_cursor = max(ref_cursor, rec["ref_end"] if svtype == "BND" and ref_len else ref_cursor)

    if ref_cursor < tlen:
        flank_len = tlen - ref_cursor
        emit_paf(
            out, sample, qlen, query_cursor, min(qlen, query_cursor + flank_len), "+",
            chrom, tlen, ref_cursor, tlen,
            [("tv", "Z", "flank"), ("sv", "Z", "terminal")])
        query_cursor += flank_len
    if query_cursor != qlen:
        print(
            f"warning: {sample} modeled query length {query_cursor} != FASTA length {qlen}",
            file=__import__("sys").stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--sample-fasta", action="append", required=True)
    parser.add_argument("--truth-vcf", action="append", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ref_lengths = read_fasta_lengths(args.reference_fasta)
    query_lengths = {}
    for fasta in args.sample_fasta:
        query_lengths.update(read_fasta_lengths(fasta))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as out:
        for vcf in args.truth_vcf:
            build_truth_paf(vcf, ref_lengths, query_lengths, out)


if __name__ == "__main__":
    main()
