#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE="${SOURCE:-$ROOT/results/HG008_chr8_alphaSat_extraction_20260727}"
OUT="${OUT:-$ROOT/results/HG008_chr8_alpha_core_only_20260727}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
SAMTOOLS="${SAMTOOLS:-$ENV_BIN/samtools}"
MIN_CORE_BP="${MIN_CORE_BP:-500000}"
MERGE_SCRIPT="$ROOT/../VallePrep_v0.0.0/src/valleprep/resources/tools/make_satellite_regions.py"

mkdir -p "$OUT/regions" "$OUT/extracted_fastas"
printf 'assembly\tcontig\tstart0\tend0\tlength_bp\n' > "$OUT/extracted_intervals.tsv"

for raw in "$SOURCE"/dna_brnn/*.alpha.bed; do
  name="$(basename "$raw" .alpha.bed)"
  fasta="$SOURCE/contigs/$name.fa"
  bed="$OUT/regions/$name.alpha_core.bed"
  python3 "$MERGE_SCRIPT" \
    --in-bed "$raw" \
    --out-bed "$bed" \
    --bin-size 100000 \
    --threshold 0.5 \
    --pad 0 \
    --fai "$fasta.fai" \
    --label-col 4 \
    --labels 2

  if [[ "$name" == chm13v2.0_chr8 ]]; then
    assembly="chm13v2.0"
    output="$OUT/extracted_fastas/chm13v2.0.chr8.alphaSat.fa"
  else
    assembly="${name%.chr8}"
    output="$OUT/extracted_fastas/$assembly.chr8.alphaSat.fa"
  fi
  : > "$output"
  while IFS=$'\t' read -r contig start end; do
    [[ -n "$contig" ]] || continue
    (( end - start >= MIN_CORE_BP )) || continue
    "$SAMTOOLS" faidx "$fasta" "${contig}:$((start + 1))-$end" |
      awk -v name="$assembly|chr8|alphaCore|source=$contig:$start-$end" \
        'NR==1 {$0=">" name} {print}' >> "$output"
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$assembly" "$contig" "$start" "$end" "$((end - start))" \
      >> "$OUT/extracted_intervals.tsv"
  done < "$bed"
  "$SAMTOOLS" faidx "$output"
done

echo "Alpha core-only outputs: $OUT"
