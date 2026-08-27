#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_DIR="${INPUT_DIR:-$ROOT/data/HG008_GIAB_Verkkov2_20260727}"
OUT_DIR="${OUT_DIR:-$ROOT/results/HG008_verkko_input_qc_20260727}"
REF="${REF:-$ROOT/../VallePrep_v0.0.0/src/valleprep/resources/reference/chm13v2.0.fa}"
THREADS="${THREADS:-16}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
MINIMAP2="${MINIMAP2:-$ENV_BIN/minimap2}"
SAMTOOLS="${SAMTOOLS:-$ENV_BIN/samtools}"
BGZIP="${BGZIP:-$ENV_BIN/bgzip}"

mkdir -p "$OUT_DIR/fastas" "$OUT_DIR/mappings" "$OUT_DIR/logs"

declare -A INPUTS=(
  [HG008_normal_haplotype1]="$INPUT_DIR/HG008_normal_haplotype1.fasta.gz"
  [HG008_normal_haplotype2]="$INPUT_DIR/HG008_normal_haplotype2.fasta.gz"
  [HG008_tumor_path1]="$INPUT_DIR/HG008_tumor_path1.fasta"
  [HG008_tumor_path2]="$INPUT_DIR/HG008_tumor_path2.fasta"
)

for sample in "${!INPUTS[@]}"; do
  input="${INPUTS[$sample]}"
  output="$OUT_DIR/fastas/$sample.fa.gz"
  if [[ ! -s "$output" ]]; then
    if [[ "$input" == *.gz ]]; then
      gzip -dc "$input" | "$BGZIP" -@ "$THREADS" -c > "$output"
    else
      "$BGZIP" -@ "$THREADS" -c "$input" > "$output"
    fi
  fi
  "$SAMTOOLS" faidx "$output"
done

python3 - "$OUT_DIR" <<'PY'
import csv
import gzip
import sys
from pathlib import Path

out_dir = Path(sys.argv[1])
rows = []
for fasta in sorted((out_dir / "fastas").glob("*.fa.gz")):
    lengths = []
    n_bases = 0
    total = 0
    with gzip.open(fasta, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq = line.strip()
            total += len(seq)
            n_bases += seq.upper().count("N")
    with open(str(fasta) + ".fai") as handle:
        lengths = [int(line.split("\t")[1]) for line in handle]
    threshold = total / 2
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            n50 = length
            break
    rows.append({
        "assembly": fasta.name[:-6],
        "contigs": len(lengths),
        "assembly_bp": total,
        "n_bases": n_bases,
        "n_fraction": n_bases / total if total else 0.0,
        "n50_bp": n50,
    })

path = out_dir / "assembly_stats.tsv"
with path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY

CHR8="$OUT_DIR/chm13v2.0.chr8.fa"
if [[ ! -s "$CHR8" ]]; then
  "$SAMTOOLS" faidx "$REF" chr8 > "$CHR8"
fi

for fasta in "$OUT_DIR"/fastas/*.fa.gz; do
  sample="$(basename "$fasta" .fa.gz)"
  paf="$OUT_DIR/mappings/$sample.chm13_chr8.paf"
  if [[ ! -s "$paf" ]]; then
    "$MINIMAP2" -t "$THREADS" -x asm20 -N 20 --secondary=yes "$fasta" "$CHR8" \
      > "$paf" 2> "$OUT_DIR/logs/$sample.minimap2.log"
  fi
done

python3 - "$OUT_DIR" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

out_dir = Path(sys.argv[1])
summary = []
for paf in sorted((out_dir / "mappings").glob("*.chm13_chr8.paf")):
    sample = paf.name[:-len(".chm13_chr8.paf")]
    data = defaultdict(lambda: {
        "contig_length": 0, "chr8_aligned_bp": 0, "primary_chr8_aligned_bp": 0,
        "best_block_bp": 0, "best_mapq": 0,
    })
    with paf.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            qstart, qend = int(fields[2]), int(fields[3])
            target, tlen = fields[5], int(fields[6])
            block, mapq = int(fields[10]), int(fields[11])
            tags = fields[12:]
            is_primary = not any(tag == "tp:A:S" for tag in tags)
            row = data[target]
            row["contig_length"] = tlen
            row["chr8_aligned_bp"] += qend - qstart
            if is_primary:
                row["primary_chr8_aligned_bp"] += qend - qstart
            if is_primary and (block, mapq) > (row["best_block_bp"], row["best_mapq"]):
                row["best_block_bp"] = block
                row["best_mapq"] = mapq
    for target, row in data.items():
        if row["primary_chr8_aligned_bp"] >= 50_000_000 or row["best_block_bp"] >= 10_000_000:
            classification = "chromosome_scale_chr8"
        elif row["primary_chr8_aligned_bp"] >= 1_000_000:
            classification = "chr8_fragment_candidate"
        elif row["primary_chr8_aligned_bp"] >= 100_000:
            classification = "repeat_associated_candidate"
        else:
            classification = "secondary_repeat_only"
        summary.append({
            "assembly": sample,
            "contig": target,
            **row,
            "classification": classification,
        })

fields = [
    "assembly", "contig", "contig_length", "primary_chr8_aligned_bp",
    "chr8_aligned_bp", "best_block_bp", "best_mapq", "classification",
]
with (out_dir / "contig_chromosome_assignments.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(summary, key=lambda r: (r["assembly"], -r["chr8_aligned_bp"])))

with (out_dir / "chr8_candidate_contigs.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(row for row in sorted(
        summary, key=lambda r: (r["assembly"], -r["chr8_aligned_bp"])
    ) if row["classification"] != "secondary_repeat_only")
PY

echo "QC outputs: $OUT_DIR"
