#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE_RUN="${SOURCE_RUN:-chr1_censat_10samples_20260827}"
RUN_NAME="${RUN_NAME:-chr1_censat_10samples_tmus_streamed_20260827}"
THREADS="${THREADS:-16}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
INPUT_DIR="$ROOT/runs/$SOURCE_RUN/input_fastas"
RUN_DIR="$ROOT/runs/$RUN_NAME"
OUTPUT_DIR="$ROOT/results/graphs"
OUTPUT_GFA="$OUTPUT_DIR/$RUN_NAME.vallescope2.gfa"

mkdir -p "$RUN_DIR/logs" "$OUTPUT_DIR"
mapfile -t FASTAS < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.fa' | sort)
[[ ${#FASTAS[@]} -eq 10 ]] || { echo "Expected 10 normalized FASTAs, found ${#FASTAS[@]}" >&2; exit 1; }

{
  printf 'sample\tpath\tsequences\tbp\n'
  for fasta in "${FASTAS[@]}"; do
    sample="$(basename "$fasta" .fa)"
    read -r sequences bp < <(awk '/^>/{n++} !/^>/{gsub(/[[:space:]]/,""); b+=length($0)} END{print n+0,b+0}' "$fasta")
    printf '%s\t%s\t%s\t%s\n' "$sample" "$fasta" "$sequences" "$bp"
  done
} > "$RUN_DIR/input_manifest.tsv"

construct_args=()
for fasta in "${FASTAS[@]}"; do
  construct_args+=( -f "$fasta" )
done

echo "[$(date -Is)] Starting revised ValleScope2 + PGGB" | tee "$RUN_DIR/logs/runner.log"
/usr/bin/time -v env PATH="$ROOT/build:$ENV_BIN:$PATH" \
  "$ENV_BIN/vsg-graph" construct \
    "${construct_args[@]}" \
    --name "$RUN_NAME" \
    --output "$OUTPUT_DIR" \
    --run-directory "$ROOT/runs" \
    --aligner vallescope2 \
    --threads "$THREADS" \
    --aligner-extra "" \
    > "$RUN_DIR/logs/vallescope2_pggb.log" 2>&1

test -s "$OUTPUT_GFA"
echo "[$(date -Is)] Completed: $OUTPUT_GFA" | tee -a "$RUN_DIR/logs/runner.log"
