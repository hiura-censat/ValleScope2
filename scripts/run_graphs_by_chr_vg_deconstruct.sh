#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GRAPH_ROOT="${GRAPH_ROOT:-$ROOT/results/graphs_by_chr}"
OUT_ROOT="${OUT_ROOT:-$ROOT/results/graphs_by_chr_vg_deconstruct_20260727}"
VG="${VG:-/home/senescence/miniconda3/envs/vallescope_dev/bin/vg}"
PYTHON="${PYTHON:-/home/senescence/miniconda3/envs/vallescope_dev/bin/python}"
THREADS="${THREADS:-8}"
BIN_BP="${BIN_BP:-50000}"

mkdir -p "$OUT_ROOT"
printf 'chromosome\tmethod\treference_path\treference_length\tall_records\tsv50_records\n' \
  > "$OUT_ROOT/summary.tsv"

for chr_dir in "$GRAPH_ROOT"/chr*; do
  [[ -d "$chr_dir" ]] || continue
  chr="$(basename "$chr_dir")"
  chr_out="$OUT_ROOT/$chr"
  mkdir -p "$chr_out"
  reference_length=0
  plot_args=()

  for method in centrolign unialigner vallescope2; do
    gfa="$chr_dir/$method.gfa"
    [[ -s "$gfa" ]] || continue
    method_out="$chr_out/$method"
    mkdir -p "$method_out"
    reference_path="$(
      awk '$1 == "P" && $2 ~ /^chm13v2\.0(\||$)/ {print $2; exit}' "$gfa"
    )"
    if [[ -z "$reference_path" ]]; then
      echo "No CHM13 reference path in $gfa" >&2
      exit 1
    fi

    "$VG" deconstruct -a -t "$THREADS" -p "$reference_path" "$gfa" \
      > "$method_out/deconstruct.all.vcf" \
      2> "$method_out/deconstruct.log"
    "$PYTHON" "$ROOT/scripts/filter_vg_deconstruct_sv.py" \
      --input "$method_out/deconstruct.all.vcf" \
      --output "$method_out/deconstruct.sv50.vcf"

    length="$(
      sed -n 's/^##contig=.*length=\([0-9][0-9]*\).*/\1/p' \
        "$method_out/deconstruct.all.vcf" | head -n 1
    )"
    reference_length="$length"
    all_count="$(grep -vc '^#' "$method_out/deconstruct.all.vcf" || true)"
    sv_count="$(grep -vc '^#' "$method_out/deconstruct.sv50.vcf" || true)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$chr" "$method" "$reference_path" "$length" "$all_count" "$sv_count" \
      >> "$OUT_ROOT/summary.tsv"
    plot_args+=(--vcf "$method_out/deconstruct.sv50.vcf" --label "$method")
  done

  MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-graphs-by-chr}" \
    "$PYTHON" "$ROOT/scripts/plot_multi_sv_reference_histogram.py" \
      "${plot_args[@]}" \
      --reference-length "$reference_length" \
      --bin-bp "$BIN_BP" \
      --region-label "CHM13 $chr alpha-satellite" \
      --output-png "$chr_out/sv_position_histogram.png" \
      --output-normalized-png "$chr_out/sv_position_histogram.normalized.png" \
      --output-tsv "$chr_out/sv_position_histogram.tsv"
done

echo "Results: $OUT_ROOT"
