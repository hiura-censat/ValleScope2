#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_DIR="${INPUT_DIR:-$ROOT/../ValleScope-graph/reproducibility/data/samples/chr8_alpha_sat}"
RUN_NAME="${RUN_NAME:-unialigner_chr8_alphaSat_10samples_20260717}"
THREADS="${THREADS:-16}"
PAIR_JOBS="${PAIR_JOBS:-2}"
RUN_DIR="$ROOT/runs/$RUN_NAME"
PAIR_DIR="$RUN_DIR/pairs"
PGGB_DIR="$RUN_DIR/pggb/unialigner"
COMBINED="$RUN_DIR/pggb_input/$RUN_NAME.fa"
ALL_PAF="$RUN_DIR/all_to_all_paf/unialigner/all.paf"
UNIALIGNER="$ROOT/.tools/UniAligner/tandem_aligner/build/bin/tandem_aligner"
CONVERTER="$ROOT/scripts/unialigner_cigar_to_paf.py"
PGGB="${PGGB:-/home/senescence/miniconda3/envs/vallescope_dev/bin/pggb}"
export PATH="$(dirname "$PGGB"):$PATH"

mkdir -p "$PAIR_DIR" "$PGGB_DIR" "$(dirname "$COMBINED")" "$(dirname "$ALL_PAF")" \
  "$ROOT/results/graphs" "$ROOT/logs"

mapfile -t FASTAS < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.fa' | sort)
if [[ ${#FASTAS[@]} -ne 10 ]]; then
  echo "Expected 10 FASTAs in $INPUT_DIR, found ${#FASTAS[@]}" >&2
  exit 1
fi

printf '%s\n' "${FASTAS[@]}" > "$RUN_DIR/input_fastas.txt"
: > "$COMBINED"
for fasta in "${FASTAS[@]}"; do
  cat "$fasta" >> "$COMBINED"
done
samtools faidx "$COMBINED"

PAIR_MANIFEST="$RUN_DIR/pairs.tsv"
: > "$PAIR_MANIFEST"
for ((i = 0; i < ${#FASTAS[@]}; ++i)); do
  for ((j = i + 1; j < ${#FASTAS[@]}; ++j)); do
    first="${FASTAS[$i]}"
    second="${FASTAS[$j]}"
    first_id="$(basename "$first" .chr8.alphaSat.all_merged.fa)"
    second_id="$(basename "$second" .chr8.alphaSat.all_merged.fa)"
    printf '%s\t%s\t%s\n' "$first" "$second" "${first_id}__${second_id}" >> "$PAIR_MANIFEST"
  done
done

run_pair() {
  local first="$1"
  local second="$2"
  local pair="$3"
  local out="$PAIR_DIR/$pair"
  mkdir -p "$out"
  if [[ ! -s "$out/cigar.txt" ]]; then
    /usr/bin/time -v "$UNIALIGNER" --first "$first" --second "$second" -o "$out" \
      > "$out/stdout.log" 2> "$out/stderr.log"
  fi
  python3 "$CONVERTER" --target "$first" --query "$second" --cigar "$out/cigar.txt" \
    > "$out/output.paf"
}
export -f run_pair
export PAIR_DIR UNIALIGNER CONVERTER
xargs -P "$PAIR_JOBS" -n 3 bash -c 'run_pair "$0" "$1" "$2"' < "$PAIR_MANIFEST"

find "$PAIR_DIR" -mindepth 2 -maxdepth 2 -name output.paf -print0 \
  | sort -z | xargs -0 cat > "$ALL_PAF"

if [[ $(wc -l < "$ALL_PAF") -ne 45 ]]; then
  echo "Expected 45 PAF records, found $(wc -l < "$ALL_PAF")" >&2
  exit 1
fi

/usr/bin/time -v "$PGGB" \
  -i "$COMBINED" -o "$PGGB_DIR" -t "$THREADS" -s 10k -p 90 \
  -a "$ALL_PAF" -n 10 -A -S \
  > "$ROOT/logs/$RUN_NAME.pggb.log" 2>&1

FINAL_GFA="$(find "$PGGB_DIR" -maxdepth 1 -name '*.smooth.final.gfa' -print -quit)"
RAW_GFA="$(find "$PGGB_DIR" -maxdepth 1 -name '*.seqwish.gfa' -print -quit)"
test -s "$FINAL_GFA"
test -s "$RAW_GFA"
cp "$FINAL_GFA" "$ROOT/results/graphs/$RUN_NAME.unialigner.gfa"
cp "$RAW_GFA" "$ROOT/results/graphs/$RUN_NAME.unialigner.raw_seqwish.gfa"

printf 'PAF\t%s\nRAW_GFA\t%s\nFINAL_GFA\t%s\n' \
  "$ALL_PAF" \
  "$ROOT/results/graphs/$RUN_NAME.unialigner.raw_seqwish.gfa" \
  "$ROOT/results/graphs/$RUN_NAME.unialigner.gfa" \
  > "$RUN_DIR/outputs.tsv"
