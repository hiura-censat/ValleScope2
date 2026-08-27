#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QC_DIR="${QC_DIR:-$ROOT/results/HG008_verkko_input_qc_20260727}"
OUT_DIR="${OUT_DIR:-$ROOT/results/HG008_chr8_alphaSat_extraction_20260727}"
THREADS="${THREADS:-8}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
SAMTOOLS="${SAMTOOLS:-$ENV_BIN/samtools}"
MINIMAP2="${MINIMAP2:-$ENV_BIN/minimap2}"
DNA_BRNN="${DNA_BRNN:-$ROOT/.tools/dna-nn/dna-brnn}"
DNA_MODEL="${DNA_MODEL:-$ROOT/.tools/dna-nn/models/attcc-alpha.knm}"
MERGE_SCRIPT="$ROOT/../VallePrep_v0.0.0/src/valleprep/resources/tools/make_satellite_regions.py"

mkdir -p "$OUT_DIR/contigs" "$OUT_DIR/dna_brnn" "$OUT_DIR/regions" \
  "$OUT_DIR/extracted_fastas" "$OUT_DIR/ribbons" "$OUT_DIR/logs"

cat > "$OUT_DIR/selected_contigs.tsv" <<'EOF'
assembly	contig
HG008_normal_haplotype1	haplotype1-0000012
HG008_normal_haplotype2	haplotype2-0000113
HG008_tumor_path1	haplotype1-0000005
HG008_tumor_path2	haplotype2-0000037
EOF

CHR8_FA="$QC_DIR/chm13v2.0.chr8.fa"
test -s "$CHR8_FA"
cp -f "$CHR8_FA" "$OUT_DIR/contigs/chm13v2.0_chr8.fa"
"$SAMTOOLS" faidx "$OUT_DIR/contigs/chm13v2.0_chr8.fa"

tail -n +2 "$OUT_DIR/selected_contigs.tsv" |
while IFS=$'\t' read -r assembly contig; do
  source="$QC_DIR/fastas/$assembly.fa.gz"
  output="$OUT_DIR/contigs/$assembly.chr8.fa"
  if [[ ! -s "$output" ]]; then
    "$SAMTOOLS" faidx "$source" "$contig" |
      awk -v name="$contig" 'NR==1 {$0=">" name} {print}' > "$output"
  fi
  "$SAMTOOLS" faidx "$output"
done

for fasta in "$OUT_DIR"/contigs/*.fa; do
  name="$(basename "$fasta" .fa)"
  raw="$OUT_DIR/dna_brnn/$name.alpha.bed"
  merged="$OUT_DIR/regions/$name.alpha_core_pad5m.bed"
  if [[ ! -s "$raw" ]]; then
    "$DNA_BRNN" -Ai "$DNA_MODEL" -t "$THREADS" -L 50 "$fasta" \
      > "$raw" 2> "$OUT_DIR/logs/$name.dna_brnn.log"
  fi
  python3 "$MERGE_SCRIPT" \
    --in-bed "$raw" \
    --out-bed "$merged" \
    --bin-size 100000 \
    --threshold 0.5 \
    --pad 5000000 \
    --fai "$fasta.fai" \
    --label-col 4 \
    --labels 2
done

: > "$OUT_DIR/extracted_intervals.tsv"
printf 'assembly\tcontig\tstart0\tend0\tlength_bp\n' >> "$OUT_DIR/extracted_intervals.tsv"
tail -n +2 "$OUT_DIR/selected_contigs.tsv" |
while IFS=$'\t' read -r assembly contig; do
  fasta="$OUT_DIR/contigs/$assembly.chr8.fa"
  bed="$OUT_DIR/regions/$assembly.chr8.alpha_core_pad5m.bed"
  output="$OUT_DIR/extracted_fastas/$assembly.chr8.alphaSat.fa"
  : > "$output"
  while IFS=$'\t' read -r chrom start end; do
    [[ -n "$chrom" ]] || continue
    region="${chrom}:$((start + 1))-$end"
    "$SAMTOOLS" faidx "$fasta" "$region" |
      awk -v name="$assembly|chr8|alphaSat|source=$chrom:$start-$end" \
        'NR==1 {$0=">" name} {print}' >> "$output"
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$assembly" "$chrom" "$start" "$end" "$((end - start))" \
      >> "$OUT_DIR/extracted_intervals.tsv"
  done < "$bed"
  "$SAMTOOLS" faidx "$output"
done

# Reference alpha-satellite sequence, retained for later pairwise alignment.
ref_bed="$OUT_DIR/regions/chm13v2.0_chr8.alpha_core_pad5m.bed"
ref_out="$OUT_DIR/extracted_fastas/chm13v2.0.chr8.alphaSat.fa"
: > "$ref_out"
while IFS=$'\t' read -r chrom start end; do
  [[ -n "$chrom" ]] || continue
  "$SAMTOOLS" faidx "$OUT_DIR/contigs/chm13v2.0_chr8.fa" \
    "${chrom}:$((start + 1))-$end" |
    awk -v name="chm13v2.0|chr8|alphaSat|source=$chrom:$start-$end" \
      'NR==1 {$0=">" name} {print}' >> "$ref_out"
done < "$ref_bed"
"$SAMTOOLS" faidx "$ref_out"

# Whole-chromosome ribbon plots used to validate contig selection.
tail -n +2 "$OUT_DIR/selected_contigs.tsv" |
while IFS=$'\t' read -r assembly contig; do
  paf="$QC_DIR/mappings/$assembly.chm13_chr8.paf"
  python3 "$ROOT/scripts/debug_paf_ribbon_png.py" \
    "$paf" "$OUT_DIR/ribbons/$assembly.chr8.full_length.png" \
    --qname chr8 \
    --tname "$contig" \
    --width 2200 \
    --height 700 \
    --tick-step-bp 10000000 \
    --title "$assembly: CHM13 chr8 vs selected assembly contig"
done

echo "Alpha-satellite extraction outputs: $OUT_DIR"
