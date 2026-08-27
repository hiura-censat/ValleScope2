#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QC_ROOT="${QC_ROOT:-$ROOT/results/HG008_verkko_input_qc_20260727}"
OUT_ROOT="${OUT_ROOT:-$ROOT/results/HG008_N1_vs_T1_satellites_by_chr_20260728}"
REF="${REF:-$ROOT/../VallePrep_v0.0.0/src/valleprep/resources/reference/chm13v2.0.fa}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
MINIMAP2="${MINIMAP2:-$ENV_BIN/minimap2}"
SAMTOOLS="${SAMTOOLS:-$ENV_BIN/samtools}"
DNA_BRNN="${DNA_BRNN:-$ROOT/.tools/dna-nn/dna-brnn}"
DNA_MODEL="${DNA_MODEL:-$ROOT/.tools/dna-nn/models/attcc-alpha.knm}"
CDR_BED="${CDR_BED:-$ROOT/data/chm13v1.1_cdr.bed}"
MERGE_SCRIPT="$ROOT/../VallePrep_v0.0.0/src/valleprep/resources/tools/make_satellite_regions.py"
THREADS="${THREADS:-16}"
CHROMS="${CHROMS:-chr8}"

mkdir -p "$OUT_ROOT"/{mappings,selected_contigs,contigs,dna_brnn,regions,extracted_fastas,ribbons,logs}

declare -A ASSEMBLIES=(
  [N1]="$QC_ROOT/fastas/HG008_normal_haplotype1.fa.gz"
  [T1]="$QC_ROOT/fastas/HG008_tumor_path1.fa.gz"
)

for sample in N1 T1; do
  assembly="${ASSEMBLIES[$sample]}"
  paf="$OUT_ROOT/mappings/$sample.chm13.paf"
  if [[ ! -s "$paf" ]]; then
    /usr/bin/time -v "$MINIMAP2" -t "$THREADS" -x asm20 \
      --secondary=no "$assembly" "$REF" \
      > "$paf" 2> "$OUT_ROOT/logs/$sample.chm13.minimap2.log"
  fi
done

python3 - "$OUT_ROOT" $CHROMS <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

out = Path(sys.argv[1])
chroms = sys.argv[2:]
rows = []
for sample in ("N1", "T1"):
    stats = defaultdict(lambda: defaultdict(lambda: {
        "primary_aligned_bp": 0,
        "all_aligned_bp": 0,
        "best_block_bp": 0,
        "best_mapq": 0,
        "contig_length": 0,
        "orientation_balance": 0,
    }))
    with (out / "mappings" / f"{sample}.chm13.paf").open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            chrom, qstart, qend = f[0], int(f[2]), int(f[3])
            strand, contig, clen = f[4], f[5], int(f[6])
            block, mapq = int(f[10]), int(f[11])
            if chrom not in chroms:
                continue
            primary = not any(tag == "tp:A:S" for tag in f[12:])
            row = stats[chrom][contig]
            span = qend - qstart
            row["contig_length"] = clen
            row["all_aligned_bp"] += span
            row["orientation_balance"] += span if strand == "+" else -span
            if primary:
                row["primary_aligned_bp"] += span
                if (block, mapq) > (row["best_block_bp"], row["best_mapq"]):
                    row["best_block_bp"] = block
                    row["best_mapq"] = mapq
    for chrom in chroms:
        candidates = stats.get(chrom, {})
        if not candidates:
            raise RuntimeError(f"no mapping candidate for {sample} {chrom}")
        contig, value = max(
            candidates.items(),
            key=lambda item: (
                item[1]["primary_aligned_bp"],
                item[1]["best_block_bp"],
                item[1]["all_aligned_bp"],
            ),
        )
        rows.append({
            "sample": sample,
            "chrom": chrom,
            "contig": contig,
            **value,
            "orientation": "+" if value["orientation_balance"] >= 0 else "-",
        })

fields = [
    "sample", "chrom", "contig", "contig_length", "primary_aligned_bp",
    "all_aligned_bp", "best_block_bp", "best_mapq", "orientation",
    "orientation_balance",
]
with (out / "selected_contigs.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY

tail -n +2 "$OUT_ROOT/selected_contigs.tsv" |
while IFS=$'\t' read -r sample chrom contig contig_length primary all best mapq orientation balance; do
  assembly="${ASSEMBLIES[$sample]}"
  contig_fa="$OUT_ROOT/contigs/$sample.$chrom.fa"
  alpha_bed="$OUT_ROOT/dna_brnn/$sample.$chrom.alpha.bed"
  if [[ ! -s "$contig_fa" ]]; then
    "$SAMTOOLS" faidx "$assembly" "$contig" |
      awk -v name="$contig" 'NR == 1 {$0 = ">" name} {print}' > "$contig_fa"
    "$SAMTOOLS" faidx "$contig_fa"
  fi
  if [[ ! -s "$alpha_bed" ]]; then
    "$DNA_BRNN" -Ai "$DNA_MODEL" -t "$THREADS" -L 50 "$contig_fa" \
      > "$alpha_bed" 2> "$OUT_ROOT/logs/$sample.$chrom.dna_brnn.log"
  fi
  merged_bed="$OUT_ROOT/regions/$sample.$chrom.alpha_components.bed"
  python3 "$MERGE_SCRIPT" \
    --in-bed "$alpha_bed" \
    --out-bed "$merged_bed" \
    --bin-size 100000 \
    --threshold 0.5 \
    --pad 0 \
    --fai "$contig_fa.fai" \
    --label-col 4 \
    --labels 2
done

python3 - "$OUT_ROOT" "$CDR_BED" <<'PY'
import csv
import sys
from pathlib import Path

out = Path(sys.argv[1])
cdr_path = Path(sys.argv[2])
cdr = {}
with cdr_path.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip().split("\t")
        cdr.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3]))

selected = list(csv.DictReader((out / "selected_contigs.tsv").open(), delimiter="\t"))
rows = []
for row in selected:
    sample, chrom, contig = row["sample"], row["chrom"], row["contig"]
    cdr_intervals = cdr.get(chrom)
    if not cdr_intervals:
        raise RuntimeError(f"no CDR interval for {chrom}")
    # Primary CDR is the longest interval when a chromosome has two annotations.
    cstart, cend, label = max(cdr_intervals, key=lambda x: x[1] - x[0])
    midpoint = (cstart + cend) // 2
    best = None
    with (out / "mappings" / f"{sample}.chm13.paf").open() as handle:
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
                if f[4] == "+":
                    projected = tstart + max(0, min(midpoint, qend) - qstart)
                else:
                    projected = tend - max(0, min(midpoint, qend) - qstart)
                best = (key, projected, distance, f[4])
    if best is None:
        raise RuntimeError(f"no projection alignment for {sample} {chrom}")
    projected = best[1]
    components = []
    with (out / "regions" / f"{sample}.{chrom}.alpha_components.bed").open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) >= 3:
                components.append((int(f[1]), int(f[2])))
    if not components:
        raise RuntimeError(f"no alpha component for {sample} {chrom}")
    start, end = min(
        components,
        key=lambda x: 0 if x[0] <= projected <= x[1] else
        min(abs(projected - x[0]), abs(projected - x[1])),
    )
    rows.append({
        "sample": sample,
        "chrom": chrom,
        "contig": contig,
        "cdr_label": label,
        "cdr_midpoint": midpoint,
        "projected_cdr_midpoint": projected,
        "projection_distance_bp": best[2],
        "mapping_strand": best[3],
        "start0": start,
        "end0": end,
        "length_bp": end - start,
        "distance_component_to_projection_bp": (
            0 if start <= projected <= end else
            min(abs(projected - start), abs(projected - end))
        ),
    })

fields = list(rows[0])
with (out / "extracted_intervals.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY

tail -n +2 "$OUT_ROOT/extracted_intervals.tsv" |
while IFS=$'\t' read -r sample chrom contig label cdr_mid projected projection_distance \
  strand start end length component_distance; do
  source="$OUT_ROOT/contigs/$sample.$chrom.fa"
  output="$OUT_ROOT/extracted_fastas/$sample.$chrom.alphaSat.fa"
  "$SAMTOOLS" faidx "$source" "$contig:$((start + 1))-$end" |
    awk -v name="$sample" 'NR == 1 {$0 = ">" name} {print}' > "$output"
  "$SAMTOOLS" faidx "$output"
  python3 "$ROOT/scripts/debug_paf_ribbon_png.py" \
    "$OUT_ROOT/mappings/$sample.chm13.paf" \
    "$OUT_ROOT/ribbons/$sample.$chrom.full_length.png" \
    --qname "$chrom" \
    --tname "$contig" \
    --width 2200 \
    --height 700 \
    --title "$sample: CHM13 $chrom vs selected assembly contig"
done

echo "Prepared satellite intervals: $OUT_ROOT"
