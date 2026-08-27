#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ODGI="${ODGI:-/home/senescence/miniconda3/envs/vallescope_dev/bin/odgi}"
OUT="${1:-$ROOT/results/odgi_chr8_alphaSat_10samples_comparison_20260820}"
THREADS="${THREADS:-8}"

mkdir -p "$OUT/logs"

METHODS=(vallescope2 vallescope2_paffy unialigner centrolign_guidetree)
GFAS=(
  "$ROOT/results/graphs/vallescope2_chr8_alphaSat_10samples_20260723.vallescope2.gfa"
  "$ROOT/results/vallescope2_chr8_alphaSat_filter_baselines_20260724/paffy_chain_tile/graph.gfa"
  "$ROOT/results/graphs/unialigner_chr8_alphaSat_10samples_20260717.unialigner.gfa"
  "$ROOT/runs/centrolign_chr8_alphaSat_guidetree_20260716/centrolign_chr8_alphaSat_guidetree.gfa"
)

if [[ ! -x "$ODGI" ]]; then
  echo "odgi executable not found: $ODGI" >&2
  exit 1
fi

printf 'method\tgfa\tnucleotides\tnodes\tedges\tpaths\tsteps\n' > "$OUT/graph_stats.tsv"
printf 'row\tsample\n' > "$OUT/path_order.tsv"
for sample in chm13v2.0 HG00146_hap1 HG00358_hap2 HG00544_hap2 HG00741_hap1 HG01114_hap1 HG01255_hap2 NA18534_hap2 NA18974_hap2 NA18982_hap1; do
  printf '%s\t%s\n' "$(( $(wc -l < "$OUT/path_order.tsv") ))" "$sample" >> "$OUT/path_order.tsv"
done

for i in "${!METHODS[@]}"; do
  method="${METHODS[$i]}"
  gfa="${GFAS[$i]}"
  dir="$OUT/$method"
  mkdir -p "$dir"

  [[ -s "$gfa" ]] || { echo "GFA not found: $gfa" >&2; exit 1; }

  "$ODGI" build -g "$gfa" -o "$dir/input.og" -O -t "$THREADS" -P \
    2> "$OUT/logs/$method.build.log"
  "$ODGI" sort -i "$dir/input.og" -o "$dir/sorted.og" -Y -t "$THREADS" -P \
    2> "$OUT/logs/$method.sort.log"

  "$ODGI" paths -i "$dir/sorted.og" -L > "$dir/paths.all.txt"
  : > "$dir/paths.txt"
  for sample in chm13v2.0 HG00146_hap1 HG00358_hap2 HG00544_hap2 HG00741_hap1 HG01114_hap1 HG01255_hap2 NA18534_hap2 NA18974_hap2 NA18982_hap1; do
    awk -v sample="$sample" '$0 == sample || index($0, sample "|") == 1 { print; exit }' \
      "$dir/paths.all.txt" >> "$dir/paths.txt"
  done
  stats="$($ODGI stats -i "$dir/sorted.og" -S | tail -n 1)"
  printf '%s\t%s\t%s\n' "$method" "$gfa" "$stats" >> "$OUT/graph_stats.tsv"
  "$ODGI" stats -i "$dir/sorted.og" -W > "$dir/components.tsv"
  "$ODGI" stats -i "$dir/sorted.og" -l -s -q > "$dir/sort_quality.tsv"

  "$ODGI" viz -i "$dir/sorted.og" -o "$dir/paths.png" \
    -p "$dir/paths.txt" -x 2400 -y 720 -a 55 -b -d -l -H -t "$THREADS" -P \
    2> "$OUT/logs/$method.viz.log"
  "$ODGI" viz -i "$dir/sorted.og" -o "$dir/depth.png" \
    -x 2400 -y 260 -O -m -B Spectral:11 -t "$THREADS" -P \
    2> "$OUT/logs/$method.depth.log"
done

if command -v montage >/dev/null 2>&1; then
  montage \
    "$OUT/vallescope2/paths.png" \
    "$OUT/vallescope2_paffy/paths.png" \
    "$OUT/unialigner/paths.png" \
    "$OUT/centrolign_guidetree/paths.png" \
    -tile 1x4 -geometry +0+12 "$OUT/paths_comparison.png"
  montage \
    "$OUT/vallescope2/depth.png" \
    "$OUT/vallescope2_paffy/depth.png" \
    "$OUT/unialigner/depth.png" \
    "$OUT/centrolign_guidetree/depth.png" \
    -tile 1x4 -geometry +0+12 "$OUT/depth_comparison.png"
fi

echo "ODGI comparison written to: $OUT"
