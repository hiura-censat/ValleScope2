#!/usr/bin/env python3
"""Project CHM13 CDRs to HG008 contigs and emit bounded extraction windows."""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--cdr-bed", required=True, type=Path)
    parser.add_argument("--window-bp", type=int, default=15_000_000)
    args = parser.parse_args()

    cdr = {}
    with args.cdr_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            cdr.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3]))

    selected = list(csv.DictReader(
        (args.root / "selected_contigs.tsv").open(), delimiter="\t"
    ))
    rows = []
    for row in selected:
        sample, chrom, contig = row["sample"], row["chrom"], row["contig"]
        cstart, cend, label = max(cdr[chrom], key=lambda x: x[1] - x[0])
        midpoint = (cstart + cend) // 2
        best = None
        with (args.root / "mappings" / f"{sample}.chm13.paf").open() as handle:
            for line in handle:
                f = line.rstrip().split("\t")
                if f[0] != chrom or f[5] != contig:
                    continue
                qstart, qend = int(f[2]), int(f[3])
                tstart, tend = int(f[7]), int(f[8])
                distance = 0 if qstart <= midpoint < qend else min(
                    abs(midpoint - qstart), abs(midpoint - qend)
                )
                primary = not any(tag == "tp:A:S" for tag in f[12:])
                key = (primary, -distance, qend - qstart, int(f[11]))
                if best is None or key > best[0]:
                    offset = max(0, min(midpoint, qend) - qstart)
                    projected = tstart + offset if f[4] == "+" else tend - offset
                    best = (key, projected, distance, f[4], int(f[11]))
        if best is None:
            raise RuntimeError(f"no CDR projection for {sample} {chrom}")
        contig_length = int(row["contig_length"])
        window_start = max(0, best[1] - args.window_bp)
        window_end = min(contig_length, best[1] + args.window_bp)
        rows.append({
            **row,
            "cdr_label": label,
            "cdr_midpoint": midpoint,
            "projected_cdr_midpoint": best[1],
            "cdr_mapping_distance_bp": best[2],
            "cdr_mapping_strand": best[3],
            "cdr_mapping_mapq": best[4],
            "window_start0": window_start,
            "window_end0": window_end,
            "window_length_bp": window_end - window_start,
        })

    with (args.root / "projected_cdr_windows.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
