#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREP="${PREP:-$ROOT/results/HG008_N1_vs_T1_satellites_by_chr_20260728}"
OUT="${OUT:-$PREP/alignments_5methods}"
THREADS="${THREADS:-16}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
MINIMAP2="${MINIMAP2:-$ENV_BIN/minimap2}"
WINNOWMAP="${WINNOWMAP:-$ENV_BIN/winnowmap}"
MERYL="${MERYL:-$ENV_BIN/meryl}"
WFMASH="${WFMASH:-$ROOT/.tools/wfmash/bin/wfmash}"
UNIALIGNER="${UNIALIGNER:-$ROOT/.tools/UniAligner/tandem_aligner/build/bin/tandem_aligner}"
VALLESCOPE2="${VALLESCOPE2:-$ROOT/build/vallescope2}"

mkdir -p "$OUT/pairs"
printf 'pair\ttarget\tquery\tchrom\n' > "$OUT/pairs.tsv"

for chrom in chr{1..22} chrX; do
  target="$PREP/extracted_fastas/N1.$chrom.alphaSat.fa"
  query="$PREP/extracted_fastas/T1.$chrom.alphaSat.fa"
  [[ -s "$target" && -s "$query" ]] || continue
  component_distance="$(awk -F'\t' -v c="$chrom" '$1=="T1" && $2==c {print $13}' "$PREP/extracted_intervals.tsv")"
  if [[ -z "$component_distance" || "$component_distance" -gt 1000000 ]]; then
    continue
  fi
  pair="${chrom}_N1_vs_T1"
  printf '%s\tN1\tT1\t%s\n' "$pair" "$chrom" >> "$OUT/pairs.tsv"
  pair_dir="$OUT/pairs/$pair"
  mkdir -p "$pair_dir"/{minimap2,winnowmap2,wfmash,unialigner,vallescope2}

  if [[ ! -s "$pair_dir/minimap2/alignment.paf" ]]; then
    /usr/bin/time -v "$MINIMAP2" -x asm5 -c --eqx --secondary=no \
      -t "$THREADS" "$target" "$query" \
      > "$pair_dir/minimap2/alignment.paf" 2> "$pair_dir/minimap2/run.log"
  fi

  repeats="$pair_dir/winnowmap2/repetitive_k19.txt"
  if [[ ! -s "$repeats" ]]; then
    "$MERYL" count k=19 output "$pair_dir/winnowmap2/target.k19.meryl" "$target" \
      > "$pair_dir/winnowmap2/meryl.count.log" 2>&1
    "$MERYL" print greater-than distinct=0.9998 \
      "$pair_dir/winnowmap2/target.k19.meryl" \
      > "$repeats" 2> "$pair_dir/winnowmap2/meryl.print.log"
  fi
  if [[ ! -e "$pair_dir/winnowmap2/timeout" &&
        ! -s "$pair_dir/winnowmap2/alignment.paf" ]]; then
    if ! /usr/bin/time -v timeout 30m \
      "$WINNOWMAP" -W "$repeats" -x asm5 -k 19 -c --eqx \
      --sv-off -f 0.01 --secondary=no -t "$THREADS" "$target" "$query" \
      > "$pair_dir/winnowmap2/alignment.paf" 2> "$pair_dir/winnowmap2/run.log"; then
      touch "$pair_dir/winnowmap2/timeout"
    fi
  fi

  if [[ ! -s "$pair_dir/wfmash/alignment.paf" ]]; then
    /usr/bin/time -v "$WFMASH" -p 90 -n 20 -l 5000 -P 100000 \
      -t "$THREADS" "$target" "$query" \
      > "$pair_dir/wfmash/alignment.paf" 2> "$pair_dir/wfmash/run.log"
  fi

  if [[ ! -e "$pair_dir/unialigner/failed" &&
        ! -s "$pair_dir/unialigner/alignment.paf" ]]; then
    if /usr/bin/time -v "$UNIALIGNER" --first "$target" --second "$query" \
      -o "$pair_dir/unialigner/work" \
      > "$pair_dir/unialigner/stdout.log" 2> "$pair_dir/unialigner/run.log"; then
      python3 "$ROOT/scripts/unialigner_cigar_to_paf.py" \
        --target "$target" --query "$query" \
        --cigar "$pair_dir/unialigner/work/cigar.txt" \
        > "$pair_dir/unialigner/alignment.paf"
    else
      touch "$pair_dir/unialigner/failed"
    fi
  fi

  if [[ ! -s "$pair_dir/vallescope2/alignment.paf" ]]; then
    PATH="$ENV_BIN:$PATH" /usr/bin/time -v "$VALLESCOPE2" --combine \
      -t "$THREADS" "$target" "$query" \
      > "$pair_dir/vallescope2/alignment.all_directions.paf" \
      2> "$pair_dir/vallescope2/run.log"
    awk '$1 == "T1" && $6 == "N1"' \
      "$pair_dir/vallescope2/alignment.all_directions.paf" \
      > "$pair_dir/vallescope2/alignment.paf"
  fi
done

echo "N1/T1 five-aligner outputs: $OUT"
