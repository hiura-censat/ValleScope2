#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE="/home/senescence/Documents/storage2/hiura/docker/bulk_env/data/00_analysis/Repthoscope.v1.0.0/sat-annotation/batch-results/chr1/work"
RUN_NAME="${RUN_NAME:-chr1_censat_10samples_20260827}"
THREADS="${THREADS:-16}"
RUN_DIR="$ROOT/runs/$RUN_NAME"
RESULT_DIR="$ROOT/results/graphs"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
CENTROLIGN="$ROOT/.tools/centrolign/bin/centrolign"
VS2_GFA="$RESULT_DIR/$RUN_NAME.vallescope2.gfa"
CENTROLIGN_GFA="$RESULT_DIR/$RUN_NAME.centrolign.gfa"

SOURCE_FASTAS=(
  "$BASE/HG00096_hap2.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00097_hap1.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00097_hap2.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00126_hap1.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00126_hap2.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00133_hap1.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00133_hap2.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00140_hap1.t2t.filtered/regions/censat-targets.fa"
  "$BASE/HG00140_hap2.t2t.filtered/regions/censat-targets.fa"
  "$BASE/chm13v2.0.t2t.filtered/regions/censat-targets.fa"
)

NORMALIZED_DIR="$RUN_DIR/input_fastas"
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/centrolign" "$NORMALIZED_DIR" "$RESULT_DIR"

FASTAS=()
for source in "${SOURCE_FASTAS[@]}"; do
  sample="$(basename "$(dirname "$(dirname "$source")")" .t2t.filtered)"
  normalized="$NORMALIZED_DIR/$sample.fa"
  if [[ ! -s "$normalized" ]]; then
    awk -v sample="$sample" 'BEGIN{seen=0} /^>/{if (++seen > 1) exit 2; print ">" sample; next} {print}' \
      "$source" > "$normalized"
  fi
  FASTAS+=("$normalized")
done

{
  printf 'sample\tpath\tsequences\tbp\n'
  for index in "${!FASTAS[@]}"; do
    fasta="${FASTAS[$index]}"
    source="${SOURCE_FASTAS[$index]}"
    test -s "$fasta"
    sample="$(basename "$fasta" .fa)"
    read -r sequences bp < <(awk '/^>/{n++} !/^>/{gsub(/[[:space:]]/,""); b+=length($0)} END{print n+0,b+0}' "$fasta")
    [[ "$sequences" -eq 1 ]] || { echo "Expected one sequence in $fasta" >&2; exit 1; }
    printf '%s\t%s -> %s\t%s\t%s\n' "$sample" "$source" "$fasta" "$sequences" "$bp"
  done
} > "$RUN_DIR/input_manifest.tsv"

echo "[$(date -Is)] Starting $RUN_NAME"
echo "[$(date -Is)] Inputs recorded in $RUN_DIR/input_manifest.tsv"

if [[ ! -s "$VS2_GFA" ]]; then
  echo "[$(date -Is)] Stage 1/2: ValleScope2 + PGGB"
  construct_args=()
  for fasta in "${FASTAS[@]}"; do
    construct_args+=( -f "$fasta" )
  done
  /usr/bin/time -v env PATH="$ROOT/build:$ENV_BIN:$PATH" \
    "$ENV_BIN/vsg-graph" construct \
      "${construct_args[@]}" \
      --name "$RUN_NAME" \
      --output "$RESULT_DIR" \
      --run-directory "$ROOT/runs" \
      --aligner vallescope2 \
      --threads "$THREADS" \
      --aligner-extra "" \
      > "$RUN_DIR/logs/vallescope2_pggb.log" 2>&1
  test -s "$VS2_GFA"
else
  echo "[$(date -Is)] Stage 1/2 already complete: $VS2_GFA"
fi

CENTROLIGN_INPUT="$RUN_DIR/centrolign/input.fa"
if [[ ! -s "$CENTROLIGN_INPUT" ]]; then
  for fasta in "${FASTAS[@]}"; do
    cat "$fasta"
  done > "$CENTROLIGN_INPUT"
fi

if [[ ! -s "$CENTROLIGN_GFA" ]]; then
  echo "[$(date -Is)] Stage 2/2: Centrolign (no guide tree)"
  export LD_LIBRARY_PATH="$ROOT/.tools/centrolign/lib:${LD_LIBRARY_PATH:-}"
  /usr/bin/time -v "$CENTROLIGN" -v 2 "$CENTROLIGN_INPUT" \
    > "$CENTROLIGN_GFA" \
    2> "$RUN_DIR/logs/centrolign.log"
  test -s "$CENTROLIGN_GFA"
else
  echo "[$(date -Is)] Stage 2/2 already complete: $CENTROLIGN_GFA"
fi

{
  printf 'method\tgfa\n'
  printf 'vallescope2_pggb\t%s\n' "$VS2_GFA"
  printf 'centrolign\t%s\n' "$CENTROLIGN_GFA"
} > "$RUN_DIR/outputs.tsv"

echo "[$(date -Is)] Completed $RUN_NAME"
