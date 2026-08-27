#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_ROOT="${OUT_ROOT:-$ROOT/results/graphs_by_chr_vg_deconstruct_20260727}"
ALPHA_ROOT="${ALPHA_ROOT:-$ROOT/data/whole_3samples_alphaSat}"
GENOME="${GENOME:-$ROOT/data/whole_3samples/chm13v2.0.t2t.filtered.fasta.gz}"
CDR_BED="${CDR_BED:-$ROOT/data/chm13v1.1_cdr.bed}"
SAMTOOLS="${SAMTOOLS:-/home/senescence/miniconda3/envs/vallescope_dev/bin/samtools}"
MINIMAP2="${MINIMAP2:-/home/senescence/miniconda3/envs/vallescope_dev/bin/minimap2}"
PYTHON="${PYTHON:-/home/senescence/miniconda3/envs/vallescope_dev/bin/python}"
THREADS="${THREADS:-8}"
BIN_BP="${BIN_BP:-50000}"

for chr_out in "$OUT_ROOT"/chr*; do
  [[ -d "$chr_out" ]] || continue
  chr="$(basename "$chr_out")"
  alpha="$ALPHA_ROOT/${chr}_alpha_sat/chm13v2.0.${chr}.alphaSat.all_merged.fa"
  [[ -s "$alpha" ]] || continue

  mapping_dir="$chr_out/cdr_projection"
  mkdir -p "$mapping_dir"
  "$SAMTOOLS" faidx "$GENOME" "$chr" > "$mapping_dir/${chr}.fa"
  "$MINIMAP2" -x asm5 --secondary=no -c -t "$THREADS" \
    "$mapping_dir/${chr}.fa" "$alpha" \
    > "$mapping_dir/alpha_to_chm13v2.paf" \
    2> "$mapping_dir/minimap2.log"
  "$PYTHON" "$ROOT/scripts/project_cdr_to_alpha_coordinates.py" \
    --paf "$mapping_dir/alpha_to_chm13v2.paf" \
    --cdr-bed "$CDR_BED" \
    --chromosome "$chr" \
    --output "$chr_out/cdr.projected.bed"

  args=()
  for method in centrolign unialigner vallescope2; do
    args+=(--vcf "$chr_out/$method/deconstruct.sv50.vcf" --label "$method")
  done
  cp "$chr_out/sv_position_histogram.png" \
    "$chr_out/sv_position_histogram.concatenated.png"
  cp "$chr_out/sv_position_histogram.normalized.png" \
    "$chr_out/sv_position_histogram.concatenated.normalized.png"
  cp "$chr_out/sv_position_histogram.tsv" \
    "$chr_out/sv_position_histogram.concatenated.tsv"
  MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-graphs-by-chr-cdr}" \
    "$PYTHON" "$ROOT/scripts/plot_multi_sv_genomic_histogram.py" \
      "${args[@]}" \
      --coordinate-map-paf "$mapping_dir/alpha_to_chm13v2.paf" \
      --cdr-bed "$CDR_BED" \
      --chromosome "$chr" \
      --bin-bp "$BIN_BP" \
      --output-png "$chr_out/sv_position_histogram.png" \
      --output-normalized-png "$chr_out/sv_position_histogram.normalized.png" \
      --output-tsv "$chr_out/sv_position_histogram.tsv"
done

echo "CDR-highlighted histograms: $OUT_ROOT"
