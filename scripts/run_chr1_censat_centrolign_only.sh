#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_NAME="${RUN_NAME:-chr1_censat_10samples_20260827}"
RUN_DIR="$ROOT/runs/$RUN_NAME"
INPUT_DIR="$RUN_DIR/input_fastas"
INPUT="$RUN_DIR/centrolign/input.fa"
OUTPUT="$ROOT/results/graphs/$RUN_NAME.centrolign.gfa"
LOG="$RUN_DIR/logs/centrolign.log"
CENTROLIGN="$ROOT/.tools/centrolign/bin/centrolign"

mkdir -p "$RUN_DIR/centrolign" "$RUN_DIR/logs" "$ROOT/results/graphs"
mapfile -t FASTAS < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.fa' | sort)
[[ ${#FASTAS[@]} -eq 10 ]] || { echo "Expected 10 normalized FASTAs, found ${#FASTAS[@]}" >&2; exit 1; }

if [[ ! -s "$INPUT" ]]; then
  for fasta in "${FASTAS[@]}"; do
    cat "$fasta"
  done > "$INPUT"
fi

echo "[$(date -Is)] Starting Centrolign without guide tree" >> "$RUN_DIR/logs/centrolign_runner.log"
export LD_LIBRARY_PATH="$ROOT/.tools/centrolign/lib:${LD_LIBRARY_PATH:-}"
/usr/bin/time -v "$CENTROLIGN" -v 2 "$INPUT" > "$OUTPUT" 2> "$LOG"
test -s "$OUTPUT"
echo "[$(date -Is)] Completed: $OUTPUT" >> "$RUN_DIR/logs/centrolign_runner.log"
