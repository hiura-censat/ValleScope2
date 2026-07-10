#!/usr/bin/env python3
"""Build a truth-like GFA/PAF directly from simulation truth VCFs.

This helper is intentionally scoped to the local vsg-benchmark simulations:
hom-alt records for DEL/INS/DUP/INV on one reference contig.  It writes a GFA
with reference segments split at truth breakpoints and sample paths threaded
through the true alleles.
"""

import argparse
import csv
from pathlib import Path


RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    name = None
    seq = []
    records = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records[name] = "".join(seq).upper()
    return records


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


def read_vcf(path):
    sample = Path(path).name.split(".")[0]
    records = []
    with open(path) as handle:
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
            svtype = info.get("SVTYPE")
            if svtype not in {"DEL", "INS", "DUP", "INV"}:
                raise ValueError(f"unsupported SVTYPE={svtype} in {path}: {line[:120]}")
            pos1 = int(fields[1])
            end1 = int(info.get("END", pos1))
            svlen = int(info.get("SVLEN", "0").split(",")[0])
            rec = {
                "sample": sample,
                "chrom": fields[0],
                "pos1": pos1,
                "end1": end1,
                "id": fields[2],
                "ref": fields[3].upper(),
                "alt": fields[4].upper(),
                "info": info,
                "svtype": svtype,
                "svlen": svlen,
            }
            if svtype == "DEL":
                # The simulation VCF's REF allele spans POS through END.  Model
                # the graph allele on that full zero-based interval and let the
                # graph deconstructor choose its own normalized anchor.
                rec["start"] = pos1 - 1
                rec["end"] = end1
            elif svtype == "INS":
                rec["start"] = pos1
                rec["end"] = pos1
            else:
                # DUP/INV intervals include POS through END in 1-based VCF.
                rec["start"] = pos1 - 1
                rec["end"] = end1
            records.append(rec)
    records.sort(key=lambda r: (r["start"], r["end"], r["id"]))
    return sample, records


def wrap(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def revcomp(seq):
    return seq.translate(RC)[::-1].upper()


class GraphBuilder:
    def __init__(self, ref_name, ref_seq, sample_sequences, sample_records):
        self.ref_name = ref_name
        self.ref_seq = ref_seq
        self.sample_sequences = sample_sequences
        self.sample_records = sample_records
        self.breakpoints = {0, len(ref_seq)}
        self.segments = {}
        self.ref_interval_to_id = {}
        self.paths = {}
        self.paf_rows = []
        self.event_rows = []
        self.next_id = 1
        self.dup_extra_as_insertion = False

    def add_breakpoints(self):
        for records in self.sample_records.values():
            for rec in records:
                self.breakpoints.add(rec["start"])
                self.breakpoints.add(rec["end"])

    def add_segment(self, seq):
        if not seq:
            return None
        sid = str(self.next_id)
        self.next_id += 1
        self.segments[sid] = seq
        return sid

    def materialize_reference_segments(self):
        bps = sorted(self.breakpoints)
        for start, end in zip(bps, bps[1:]):
            if end <= start:
                continue
            sid = self.add_segment(self.ref_seq[start:end])
            self.ref_interval_to_id[(start, end)] = sid
        self.ref_intervals = [
            (start, end, self.ref_interval_to_id[(start, end)])
            for start, end in zip(bps, bps[1:])
            if end > start
        ]

    def ref_steps(self, start, end, orient="+"):
        steps = []
        for istart, iend, sid in self.ref_intervals:
            if iend <= start:
                continue
            if istart >= end:
                break
            if istart < start or iend > end:
                raise ValueError(
                    f"interval [{start}, {end}) is not aligned to reference segments")
            steps.append((sid, orient, istart, iend))
        if orient == "-":
            steps = [(sid, "-", istart, iend) for sid, _, istart, iend in reversed(steps)]
        return steps

    def add_ref_path(self):
        self.paths[self.ref_name] = [(sid, "+") for _, _, sid in self.ref_intervals]

    def add_paf(self, sample, qlen, qstart, qend, strand, tstart, tend, svtype, svid):
        if qend <= qstart or tend <= tstart:
            return
        block = max(qend - qstart, tend - tstart)
        match = min(qend - qstart, tend - tstart)
        self.paf_rows.append([
            sample, qlen, qstart, qend, strand,
            self.ref_name, len(self.ref_seq), tstart, tend,
            match, block, 60, f"tv:Z:{svtype}", f"sv:Z:{svid}",
        ])

    def add_sample_path(self, sample, records):
        sample_seq = self.sample_sequences[sample]
        qlen = len(sample_seq)
        steps = []
        ref_cursor = 0
        query_cursor = 0
        for rec in records:
            start = rec["start"]
            end = rec["end"]
            if start < ref_cursor:
                raise ValueError(
                    f"overlapping truth records in {sample}: {rec['id']} starts "
                    f"{start} before cursor {ref_cursor}")

            for sid, orient, istart, iend in self.ref_steps(ref_cursor, start):
                steps.append((sid, orient))
            flank_len = start - ref_cursor
            self.add_paf(sample, qlen, query_cursor, query_cursor + flank_len,
                         "+", ref_cursor, start, "flank", rec["id"])
            query_cursor += flank_len
            ref_cursor = start

            svtype = rec["svtype"]
            event_qstart = query_cursor
            if svtype == "DEL":
                ref_cursor = end
            elif svtype == "INS":
                inserted = rec["alt"][1:]
                sid = self.add_segment(inserted)
                if sid:
                    steps.append((sid, "+"))
                query_cursor += len(inserted)
            elif svtype == "DUP":
                dup_steps = self.ref_steps(start, end, "+")
                for sid, orient, _, _ in dup_steps:
                    steps.append((sid, orient))
                if self.dup_extra_as_insertion:
                    sid = self.add_segment(self.ref_seq[start:end])
                    if sid:
                        steps.append((sid, "+"))
                else:
                    for sid, orient, _, _ in dup_steps:
                        steps.append((sid, orient))
                query_cursor += 2 * (end - start)
                ref_cursor = end
            elif svtype == "INV":
                inv_steps = self.ref_steps(start, end, "-")
                for sid, orient, _, _ in inv_steps:
                    steps.append((sid, orient))
                query_cursor += end - start
                ref_cursor = end
            event_qend = query_cursor
            self.add_paf(sample, qlen, event_qstart, event_qend,
                         "-" if svtype == "INV" else "+",
                         start, end, svtype, rec["id"])
            self.event_rows.append({
                "sample": sample,
                "id": rec["id"],
                "svtype": svtype,
                "ref_start": start,
                "ref_end": end,
                "query_start": event_qstart,
                "query_end": event_qend,
            })

        for sid, orient, _, _ in self.ref_steps(ref_cursor, len(self.ref_seq)):
            steps.append((sid, orient))
        flank_len = len(self.ref_seq) - ref_cursor
        self.add_paf(sample, qlen, query_cursor, min(qlen, query_cursor + flank_len),
                     "+", ref_cursor, len(self.ref_seq), "flank", "terminal")
        query_cursor += flank_len
        self.paths[sample] = steps

    def spelled_path(self, steps):
        chunks = []
        for sid, orient in steps:
            seq = self.segments[sid]
            chunks.append(seq if orient == "+" else revcomp(seq))
        return "".join(chunks)

    def validate_paths(self):
        failures = []
        ref_spelled = self.spelled_path(self.paths[self.ref_name])
        if ref_spelled != self.ref_seq:
            failures.append((self.ref_name, len(ref_spelled), len(self.ref_seq)))
        for sample, expected in self.sample_sequences.items():
            got = self.spelled_path(self.paths[sample])
            if got != expected:
                first = next((i for i, (a, b) in enumerate(zip(got, expected)) if a != b),
                             min(len(got), len(expected)))
                failures.append((sample, len(got), len(expected), first))
        return failures

    def build(self):
        self.add_breakpoints()
        self.materialize_reference_segments()
        self.add_ref_path()
        for sample, records in self.sample_records.items():
            self.add_sample_path(sample, records)

    def write_gfa(self, path):
        with open(path, "w") as out:
            out.write("H\tVN:Z:1.0\n")
            for sid in sorted(self.segments, key=lambda x: int(x)):
                out.write(f"S\t{sid}\t{self.segments[sid]}\n")
            links = set()
            for steps in self.paths.values():
                for left, right in zip(steps, steps[1:]):
                    links.add((left[0], left[1], right[0], right[1]))
            for left, lo, right, ro in sorted(links, key=lambda x: tuple(int(v) if v.isdigit() else v for v in (x[0], x[2]))):
                out.write(f"L\t{left}\t{lo}\t{right}\t{ro}\t0M\n")
            for name, steps in self.paths.items():
                walk = ",".join(f"{sid}{orient}" for sid, orient in steps)
                overlaps = ",".join(["0M"] * max(0, len(steps) - 1))
                out.write(f"P\t{name}\t{walk}\t{overlaps}\n")

    def write_paf(self, path):
        with open(path, "w") as out:
            for row in self.paf_rows:
                out.write("\t".join(map(str, row)) + "\n")

    def write_event_map(self, path):
        fields = ["sample", "id", "svtype", "ref_start", "ref_end", "query_start", "query_end"]
        with open(path, "w", newline="") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            writer.writerows(self.event_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--sample-fasta", action="append", required=True)
    parser.add_argument("--truth-vcf", action="append", required=True)
    parser.add_argument("--out-gfa", required=True)
    parser.add_argument("--out-paf", required=True)
    parser.add_argument("--out-map", required=True)
    parser.add_argument("--validate", action="store_true")
    parser.add_argument(
        "--dup-extra-as-insertion", action="store_true",
        help="Represent the extra DUP copy as a private inserted segment.")
    args = parser.parse_args()

    ref_records = read_fasta(args.reference_fasta)
    if len(ref_records) != 1:
        raise ValueError("expected exactly one reference FASTA record")
    ref_name, ref_seq = next(iter(ref_records.items()))

    sample_sequences = {}
    for fasta in args.sample_fasta:
        sample_sequences.update(read_fasta(fasta))

    sample_records = {}
    for vcf in args.truth_vcf:
        sample, records = read_vcf(vcf)
        sample_records[sample] = records

    missing = set(sample_records) - set(sample_sequences)
    if missing:
        raise ValueError(f"missing sample FASTA records: {sorted(missing)}")

    builder = GraphBuilder(ref_name, ref_seq, sample_sequences, sample_records)
    builder.dup_extra_as_insertion = args.dup_extra_as_insertion
    builder.build()
    failures = builder.validate_paths()
    if failures:
        message = "truth-like GFA path validation failed: " + repr(failures[:5])
        if args.validate:
            raise RuntimeError(message)
        print("warning:", message)

    Path(args.out_gfa).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_paf).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_map).parent.mkdir(parents=True, exist_ok=True)
    builder.write_gfa(args.out_gfa)
    builder.write_paf(args.out_paf)
    builder.write_event_map(args.out_map)
    print(f"wrote {args.out_gfa}")
    print(f"wrote {args.out_paf}")
    print(f"wrote {args.out_map}")


if __name__ == "__main__":
    main()
