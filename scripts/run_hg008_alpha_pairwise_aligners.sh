#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_ROOT="${INPUT_ROOT:-$ROOT/results/HG008_chr8_alphaSat_extraction_20260727/extracted_fastas}"
OUT_ROOT="${OUT_ROOT:-$ROOT/results/HG008_alpha_pairwise_5aligners_20260727}"
THREADS="${THREADS:-16}"
RUN_VALLESCOPE="${RUN_VALLESCOPE:-1}"
ENV_BIN="${ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
MINIMAP2="${MINIMAP2:-$ENV_BIN/minimap2}"
WINNOWMAP="${WINNOWMAP:-$ENV_BIN/winnowmap}"
MERYL="${MERYL:-$ENV_BIN/meryl}"
WFMASH="${WFMASH:-$ROOT/.tools/wfmash/bin/wfmash}"
UNIALIGNER="${UNIALIGNER:-$ROOT/.tools/UniAligner/tandem_aligner/build/bin/tandem_aligner}"
VALLESCOPE2="${VALLESCOPE2:-$ROOT/build/vallescope2}"

mkdir -p "$OUT_ROOT/inputs" "$OUT_ROOT/pairs" "$OUT_ROOT/logs"

extract_record() {
  local source="$1"
  local record="$2"
  local name="$3"
  local output="$4"
  "$ENV_BIN/samtools" faidx "$source" "$record" |
    awk -v name="$name" 'NR==1 {$0=">" name} {print}' > "$output"
  "$ENV_BIN/samtools" faidx "$output"
}

n1_record="$(cut -f1 "$INPUT_ROOT/HG008_normal_haplotype1.chr8.alphaSat.fa.fai")"
n2_record="$(cut -f1 "$INPUT_ROOT/HG008_normal_haplotype2.chr8.alphaSat.fa.fai")"
t1_record="$(cut -f1 "$INPUT_ROOT/HG008_tumor_path1.chr8.alphaSat.fa.fai")"
t2a_record="$(sed -n '1s/\t.*//p' "$INPUT_ROOT/HG008_tumor_path2.chr8.alphaSat.fa.fai")"
t2b_record="$(sed -n '2s/\t.*//p' "$INPUT_ROOT/HG008_tumor_path2.chr8.alphaSat.fa.fai")"

extract_record "$INPUT_ROOT/HG008_normal_haplotype1.chr8.alphaSat.fa" "$n1_record" N1 "$OUT_ROOT/inputs/N1.fa"
extract_record "$INPUT_ROOT/HG008_normal_haplotype2.chr8.alphaSat.fa" "$n2_record" N2 "$OUT_ROOT/inputs/N2.fa"
extract_record "$INPUT_ROOT/HG008_tumor_path1.chr8.alphaSat.fa" "$t1_record" T1 "$OUT_ROOT/inputs/T1.fa"
extract_record "$INPUT_ROOT/HG008_tumor_path2.chr8.alphaSat.fa" "$t2a_record" T2a "$OUT_ROOT/inputs/T2a.fa"
extract_record "$INPUT_ROOT/HG008_tumor_path2.chr8.alphaSat.fa" "$t2b_record" T2b "$OUT_ROOT/inputs/T2b.fa"

cat > "$OUT_ROOT/pairs.tsv" <<'EOF'
pair	target	query	tumor_group
N1_vs_N2	N1	N2	normal
N1_vs_T1	N1	T1	T1
N1_vs_T2a	N1	T2a	T2
N1_vs_T2b	N1	T2b	T2
N2_vs_T1	N2	T1	T1
N2_vs_T2a	N2	T2a	T2
N2_vs_T2b	N2	T2b	T2
EOF

run_pair() {
  local pair="$1"
  local target_name="$2"
  local query_name="$3"
  local target="$OUT_ROOT/inputs/$target_name.fa"
  local query="$OUT_ROOT/inputs/$query_name.fa"
  local pair_dir="$OUT_ROOT/pairs/$pair"
  mkdir -p "$pair_dir"/{minimap2,winnowmap2,wfmash,unialigner,vallescope2}

  if [[ ! -s "$pair_dir/minimap2/alignment.paf" ]]; then
    /usr/bin/time -v "$MINIMAP2" -x asm5 -c --eqx --secondary=no \
      -t "$THREADS" "$target" "$query" \
      > "$pair_dir/minimap2/alignment.paf" \
      2> "$pair_dir/minimap2/run.log"
  fi

  local meryl_db="$pair_dir/winnowmap2/target.k19.meryl"
  local repeats="$pair_dir/winnowmap2/repetitive_k19.txt"
  if [[ ! -s "$repeats" ]]; then
    "$MERYL" count k=19 output "$meryl_db" "$target" \
      > "$pair_dir/winnowmap2/meryl.count.log" 2>&1
    "$MERYL" print greater-than distinct=0.9998 "$meryl_db" \
      > "$repeats" 2> "$pair_dir/winnowmap2/meryl.print.log"
  fi
  if [[ ! -s "$pair_dir/winnowmap2/alignment.paf" ]]; then
    /usr/bin/time -v "$WINNOWMAP" -W "$repeats" -x asm5 -k 19 -c --eqx \
      --sv-off -f 0.01 --secondary=no -t "$THREADS" "$target" "$query" \
      > "$pair_dir/winnowmap2/alignment.paf" \
      2> "$pair_dir/winnowmap2/run.log"
  fi

  if [[ ! -s "$pair_dir/wfmash/alignment.paf" ]]; then
    /usr/bin/time -v "$WFMASH" -p 90 -n 20 -l 5000 -P 100000 \
      -t "$THREADS" "$target" "$query" \
      > "$pair_dir/wfmash/alignment.paf" \
      2> "$pair_dir/wfmash/run.log"
  fi

  if [[ ! -s "$pair_dir/unialigner/alignment.paf" ]]; then
    /usr/bin/time -v "$UNIALIGNER" --first "$target" --second "$query" \
      -o "$pair_dir/unialigner/work" \
      > "$pair_dir/unialigner/stdout.log" \
      2> "$pair_dir/unialigner/run.log"
    python3 "$ROOT/scripts/unialigner_cigar_to_paf.py" \
      --target "$target" \
      --query "$query" \
      --cigar "$pair_dir/unialigner/work/cigar.txt" \
      > "$pair_dir/unialigner/alignment.paf"
  fi

  if [[ "$RUN_VALLESCOPE" == 1 && ! -s "$pair_dir/vallescope2/alignment.paf" ]]; then
    PATH="$ENV_BIN:$PATH" /usr/bin/time -v "$VALLESCOPE2" --combine \
      -t "$THREADS" "$target" "$query" \
      > "$pair_dir/vallescope2/alignment.all_directions.paf" \
      2> "$pair_dir/vallescope2/run.log"
    awk -v query="$query_name" -v target="$target_name" \
      '$1 == query && $6 == target' \
      "$pair_dir/vallescope2/alignment.all_directions.paf" \
      > "$pair_dir/vallescope2/alignment.paf"
  fi
}

export -f run_pair
export ROOT OUT_ROOT THREADS RUN_VALLESCOPE ENV_BIN MINIMAP2 WINNOWMAP MERYL WFMASH UNIALIGNER VALLESCOPE2

tail -n +2 "$OUT_ROOT/pairs.tsv" |
while IFS=$'\t' read -r pair target query tumor_group; do
  echo "== $pair =="
  run_pair "$pair" "$target" "$query"
done

echo "Pairwise alignments: $OUT_ROOT"
