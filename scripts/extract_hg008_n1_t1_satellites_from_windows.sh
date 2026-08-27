#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_ROOT="${OUT_ROOT:-$ROOT/results/HG008_N1_vs_T1_satellites_by_chr_20260728}"
QC_ROOT="${QC_ROOT:-$ROOT/results/HG008_verkko_input_qc_20260727}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
SAMTOOLS="${SAMTOOLS:-$ENV_BIN/samtools}"
DNA_BRNN="${DNA_BRNN:-$ROOT/.tools/dna-nn/dna-brnn}"
DNA_MODEL="${DNA_MODEL:-$ROOT/.tools/dna-nn/models/attcc-alpha.knm}"
MERGE_SCRIPT="$ROOT/../VallePrep_v0.0.0/src/valleprep/resources/tools/make_satellite_regions.py"
THREADS="${THREADS:-16}"

declare -A ASSEMBLIES=(
  [N1]="$QC_ROOT/fastas/HG008_normal_haplotype1.fa.gz"
  [T1]="$QC_ROOT/fastas/HG008_tumor_path1.fa.gz"
)

python3 "$ROOT/scripts/extract_hg008_n1_t1_satellite_windows.py" \
  --root "$OUT_ROOT" \
  --cdr-bed "$ROOT/data/chm13v1.1_cdr.bed"

mkdir -p "$OUT_ROOT"/{satellite_windows,dna_brnn_windows,regions_windows,extracted_fastas,logs}
printf 'sample\tchrom\tcontig\twindow_start0\twindow_end0\tcomponent_start0\tcomponent_end0\tassembly_start0\tassembly_end0\tlength_bp\tcdr_mapping_distance_bp\tcdr_mapping_mapq\tcomponent_distance_bp\talpha_bin_threshold\n' \
  > "$OUT_ROOT/extracted_intervals.tsv"
printf 'sample\tchrom\tcontig\treason\tcdr_mapping_distance_bp\tcdr_mapping_mapq\n' \
  > "$OUT_ROOT/skipped_intervals.tsv"

tail -n +2 "$OUT_ROOT/projected_cdr_windows.tsv" |
while IFS=$'\t' read -r sample chrom contig contig_length primary all best best_mapq \
  orientation balance label cdr_mid projected map_distance map_strand mapq \
  window_start window_end window_length; do
  window_fa="$OUT_ROOT/satellite_windows/$sample.$chrom.fa"
  raw_bed="$OUT_ROOT/dna_brnn_windows/$sample.$chrom.alpha.bed"
  merged_bed="$OUT_ROOT/regions_windows/$sample.$chrom.alpha_components.bed"
  full_marker="$OUT_ROOT/satellite_windows/$sample.$chrom.full_contig"
  source="${ASSEMBLIES[$sample]}"
  if [[ -e "$full_marker" ]]; then
    window_start=0
    window_end="$contig_length"
  fi
  if [[ ! -s "$window_fa" ]]; then
    "$SAMTOOLS" faidx "$source" "$contig:$((window_start + 1))-$window_end" |
      awk -v name="$sample.$chrom.window" 'NR == 1 {$0 = ">" name} {print}' \
      > "$window_fa"
    "$SAMTOOLS" faidx "$window_fa"
  fi
  if [[ ! -s "$raw_bed" ]]; then
    "$DNA_BRNN" -Ai "$DNA_MODEL" -t "$THREADS" -L 50 "$window_fa" \
      > "$raw_bed" 2> "$OUT_ROOT/logs/$sample.$chrom.window.dna_brnn.log"
  fi
  alpha_threshold=0.5
  for threshold in 0.5 0.25 0.1; do
    python3 "$MERGE_SCRIPT" \
      --in-bed "$raw_bed" --out-bed "$merged_bed" \
      --bin-size 100000 --threshold "$threshold" --pad 0 \
      --fai "$window_fa.fai" --label-col 4 --labels 2
    if [[ -s "$merged_bed" ]]; then
      alpha_threshold="$threshold"
      break
    fi
  done
  if [[ ! -s "$merged_bed" ]]; then
    window_start=0
    window_end="$contig_length"
    "$SAMTOOLS" faidx "$source" "$contig" |
      awk -v name="$sample.$chrom.window" 'NR == 1 {$0 = ">" name} {print}' \
      > "$window_fa"
    "$SAMTOOLS" faidx "$window_fa"
    : > "$full_marker"
    "$DNA_BRNN" -Ai "$DNA_MODEL" -t "$THREADS" -L 50 "$window_fa" \
      > "$raw_bed" 2> "$OUT_ROOT/logs/$sample.$chrom.full_contig.dna_brnn.log"
    for threshold in 0.5 0.25 0.1; do
      python3 "$MERGE_SCRIPT" \
        --in-bed "$raw_bed" --out-bed "$merged_bed" \
        --bin-size 100000 --threshold "$threshold" --pad 0 \
        --fai "$window_fa.fai" --label-col 4 --labels 2
      if [[ -s "$merged_bed" ]]; then
        alpha_threshold="$threshold"
        break
      fi
    done
  fi
  if [[ ! -s "$merged_bed" ]]; then
    printf '%s\t%s\t%s\tno_dense_alpha_component\t%s\t%s\n' \
      "$sample" "$chrom" "$contig" "$map_distance" "$mapq" \
      >> "$OUT_ROOT/skipped_intervals.tsv"
    continue
  fi

  read -r component_start component_end component_distance < <(
    python3 - "$merged_bed" "$((projected - window_start))" <<'PY'
import sys
point = int(sys.argv[2])
parts = []
for line in open(sys.argv[1]):
    f = line.rstrip().split("\t")
    if len(f) >= 3:
        start, end = int(f[1]), int(f[2])
        distance = 0 if start <= point <= end else min(abs(point-start), abs(point-end))
        parts.append((distance, -(end-start), start, end))
if not parts:
    raise SystemExit("no alpha component")
distance, _, start, end = min(parts)
print(start, end, distance)
PY
  )
  output="$OUT_ROOT/extracted_fastas/$sample.$chrom.alphaSat.fa"
  record="$(cut -f1 "$window_fa.fai")"
  "$SAMTOOLS" faidx "$window_fa" "$record:$((component_start + 1))-$component_end" |
    awk -v name="$sample" 'NR == 1 {$0 = ">" name} {print}' > "$output"
  "$SAMTOOLS" faidx "$output"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$chrom" "$contig" "$window_start" "$window_end" \
    "$component_start" "$component_end" \
    "$((window_start + component_start))" "$((window_start + component_end))" \
    "$((component_end - component_start))" "$map_distance" "$mapq" \
    "$component_distance" "$alpha_threshold" >> "$OUT_ROOT/extracted_intervals.tsv"
done

echo "Extracted N1/T1 satellite components: $OUT_ROOT"
