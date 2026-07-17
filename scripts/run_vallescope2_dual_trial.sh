#!/usr/bin/env bash
# Run one ValleScope2 parameter trial on real alpha-satellite and truth simulation data.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <output_name> [vallescope2 options...]" >&2
    exit 1
fi

OUTPUT_NAME="$1"
shift
VALLESCOPE2_OPTS=("$@")

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_BIN="${VS2_TRIAL_ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
THREADS="${VS2_DUAL_TRIAL_THREADS:-8}"
ALPHA_DATA="$ROOT/data/chr8_alphaSat_chm13_3samples_gap_analysis_20260717/input_fastas"
SIM_DATA="$ROOT/data/test/chr1_alphaSat_1M_no_translocation"

if [[ "$OUTPUT_NAME" = /* ]]; then
    OUTPUT_DIR="$OUTPUT_NAME"
else
    OUTPUT_DIR="$ROOT/results/parameter_trials/$OUTPUT_NAME"
fi
if [[ -e "$OUTPUT_DIR" ]]; then
    echo "Refusing to overwrite existing output: $OUTPUT_DIR" >&2
    exit 1
fi

ALPHA_DIR="$OUTPUT_DIR/alpha_sat"
SIM_DIR="$OUTPUT_DIR/simulation"
mkdir -p "$ALPHA_DIR/logs" "$ALPHA_DIR/graphs" "$ALPHA_DIR/runs" \
         "$SIM_DIR/logs" "$SIM_DIR/graphs" "$SIM_DIR/runs"

status_file="$OUTPUT_DIR/STATUS"
printf 'running\n' > "$status_file"
on_error() {
    printf 'failed\n' > "$status_file"
    echo "FAILED: see logs under $OUTPUT_DIR" >&2
}
trap on_error ERR

printf '%q ' "$ROOT/scripts/run_vallescope2_dual_trial.sh" "$OUTPUT_NAME" \
    "${VALLESCOPE2_OPTS[@]}" > "$OUTPUT_DIR/command.sh"
printf '\n' >> "$OUTPUT_DIR/command.sh"
chmod +x "$OUTPUT_DIR/command.sh"
printf '%s\n' "${VALLESCOPE2_OPTS[@]:-}" > "$OUTPUT_DIR/options.txt"
git -C "$ROOT" rev-parse HEAD > "$OUTPUT_DIR/git_commit.txt"
git -C "$ROOT" status --short > "$OUTPUT_DIR/git_status.txt"
"$ROOT/build/vallescope2" --version > "$OUTPUT_DIR/vallescope2_version.txt"

ALPHA_FASTAS=(
    "$ALPHA_DATA/chm13v2.0.chr8.alphaSat.all_merged.fa"
    "$ALPHA_DATA/HG01114_hap1.chr8.alphaSat.all_merged.fa"
    "$ALPHA_DATA/NA18534_hap2.chr8.alphaSat.all_merged.fa"
    "$ALPHA_DATA/NA18974_hap2.chr8.alphaSat.all_merged.fa"
)
SIM_FASTAS=(
    "$SIM_DATA/input_fastas_renamed/reference.fa"
    "$SIM_DATA/input_fastas_renamed/sample01.fa"
    "$SIM_DATA/input_fastas_renamed/sample02.fa"
    "$SIM_DATA/input_fastas_renamed/sample03.fa"
)

{
    printf 'role\tdataset\tpath\tsha256\n'
    for fasta in "${ALPHA_FASTAS[@]}"; do
        printf 'input\talpha_sat\t%s\t%s\n' "$fasta" "$(sha256sum "$fasta" | cut -d' ' -f1)"
    done
    for fasta in "${SIM_FASTAS[@]}"; do
        printf 'input\tsimulation\t%s\t%s\n' "$fasta" "$(sha256sum "$fasta" | cut -d' ' -f1)"
    done
} > "$OUTPUT_DIR/input_manifest.tsv"

echo "== dual ValleScope2 trial =="
echo "output: $OUTPUT_DIR"
echo "options: ${VALLESCOPE2_OPTS[*]:-(default)}"
echo

ALPHA_NAME="alpha_sat"
ALPHA_GFA="$ALPHA_DIR/graphs/${ALPHA_NAME}.vallescope2.gfa"
ALPHA_ALL_PAF="$ALPHA_DIR/runs/${ALPHA_NAME}/all_to_all_paf/vallescope2/all.paf"
echo "== alpha-satellite: ValleScope2 + PGGB =="
alpha_start=$(date +%s)
PATH="$ROOT/build:$ENV_BIN:$PATH" "$ENV_BIN/vsg-graph" construct \
    -f "${ALPHA_FASTAS[0]}" -f "${ALPHA_FASTAS[1]}" \
    -f "${ALPHA_FASTAS[2]}" -f "${ALPHA_FASTAS[3]}" \
    --name "$ALPHA_NAME" \
    --output "$ALPHA_DIR/graphs" \
    --run-directory "$ALPHA_DIR/runs" \
    --aligner vallescope2 \
    --threads "$THREADS" \
    --aligner-extra "${VALLESCOPE2_OPTS[*]:-}" \
    > "$ALPHA_DIR/logs/construct.log" 2>&1
alpha_seconds=$(($(date +%s) - alpha_start))
test -s "$ALPHA_GFA"
test -s "$ALPHA_ALL_PAF"
echo "alpha construct: ok (${alpha_seconds}s)"

"$ENV_BIN/python" "$ROOT/scripts/prepare_chm13_reference_paf.py" \
    --input "$ALPHA_ALL_PAF" \
    --output "$ALPHA_DIR/reference.paf" \
    --reference-token chm13 \
    --reference-name chm13v2.0 \
    > "$ALPHA_DIR/logs/prepare_reference_paf.log"
"$ENV_BIN/python" "$ROOT/scripts/plot_ref_query_full_length_overview.py" \
    --paf "$ALPHA_DIR/reference.paf" \
    --label "ValleScope2" \
    --target-label "CHM13" \
    --output-dir "$ALPHA_DIR/ribbons"

echo "== simulation: ValleScope2 + PGGB + truth benchmark =="
sim_start=$(date +%s)
VS2_TRIAL_ENV_BIN="$ENV_BIN" \
VS2_TRIAL_THREADS="$THREADS" \
VS2_TRIAL_SIM_DATA="$SIM_DATA" \
VS2_TRIAL_EXACT_RUN_NAME="simulation" \
VS2_TRIAL_GRAPH_DIR="$SIM_DIR/graphs" \
VS2_TRIAL_RUN_DIR="$SIM_DIR/runs" \
VS2_TRIAL_BENCHMARK_DIR="$SIM_DIR/benchmark" \
VS2_TRIAL_LOG_DIR="$SIM_DIR/logs" \
VS2_TRIAL_SWEEP_LOG="$OUTPUT_DIR/simulation_trials.tsv" \
    "$ROOT/scripts/run_vallescope2_trial.sh" simulation \
    "${VALLESCOPE2_OPTS[@]}" | tee "$SIM_DIR/logs/runner.stdout.log"
sim_seconds=$(($(date +%s) - sim_start))

SIM_ALL_PAF="$SIM_DIR/runs/simulation/all_to_all_paf/vallescope2/all.paf"
test -s "$SIM_ALL_PAF"
test -s "$SIM_DIR/benchmark/benchmark_summary.tsv"

mapfile -t RAW_SEQWISH_GFAS < <(
    find "$SIM_DIR/runs/simulation/pggb/vallescope2" -maxdepth 1 \
        -type f -name '*.seqwish.gfa' -print)
if [[ ${#RAW_SEQWISH_GFAS[@]} -ne 1 ]]; then
    echo "Expected exactly one raw seqwish GFA, found ${#RAW_SEQWISH_GFAS[@]}" >&2
    exit 1
fi
RAW_BENCHMARK_DIR="$SIM_DIR/benchmark_raw_seqwish"
echo "== simulation raw seqwish truth benchmark =="
raw_bench_start=$(date +%s)
"$ENV_BIN/vsg-benchmark" simulate benchmark-graph-vcf \
    --gfa "${RAW_SEQWISH_GFAS[0]}" \
    --reference-path reference \
    --reference-fasta "$SIM_DATA/truth_renamed/reference/visor_reference.fa" \
    --truth-vcf "$SIM_DATA/truth_renamed/vcfs/sample01.truth.vcf" \
    --truth-vcf "$SIM_DATA/truth_renamed/vcfs/sample02.truth.vcf" \
    --truth-vcf "$SIM_DATA/truth_renamed/vcfs/sample03.truth.vcf" \
    --output "$RAW_BENCHMARK_DIR" \
    --all-snarls --refdist 500 --pctseq 0.7 --pctsize 0.7 \
    --type-ignore --dup-to-ins --threads "$THREADS" \
    --vg "$ENV_BIN/vg" --bcftools "$ENV_BIN/bcftools" \
    --bgzip "$ENV_BIN/bgzip" --tabix "$ENV_BIN/tabix" \
    --truvari "$ENV_BIN/truvari" \
    > "$SIM_DIR/logs/simulation.raw_seqwish.bench.log" 2>&1
raw_benchmark_seconds=$(($(date +%s) - raw_bench_start))
test -s "$RAW_BENCHMARK_DIR/benchmark_summary.tsv"
echo "raw seqwish benchmark: ok (${raw_benchmark_seconds}s)"

"$ENV_BIN/python" "$ROOT/scripts/prepare_chm13_reference_paf.py" \
    --input "$SIM_ALL_PAF" \
    --output "$SIM_DIR/reference.paf" \
    --reference-token reference \
    --reference-name reference \
    > "$SIM_DIR/logs/prepare_reference_paf.log"
"$ENV_BIN/python" "$ROOT/scripts/plot_ref_query_full_length_overview.py" \
    --paf "$SIM_DIR/reference.paf" \
    --label "ValleScope2" \
    --target-label "reference" \
    --output-dir "$SIM_DIR/ribbons"

"$ENV_BIN/python" "$ROOT/scripts/summarize_vallescope2_dual_trial.py" \
    --alpha-coverage "$ALPHA_DIR/ribbons/coverage_summary.tsv" \
    --simulation-coverage "$SIM_DIR/ribbons/coverage_summary.tsv" \
    --benchmark "$SIM_DIR/benchmark/benchmark_summary.tsv" \
    --raw-benchmark "$RAW_BENCHMARK_DIR/benchmark_summary.tsv" \
    --output "$OUTPUT_DIR/trial_summary.tsv"

{
    printf 'metric\tvalue\n'
    printf 'output_dir\t%s\n' "$OUTPUT_DIR"
    printf 'threads\t%s\n' "$THREADS"
    printf 'options\t%s\n' "${VALLESCOPE2_OPTS[*]:-(default)}"
    printf 'alpha_construct_seconds\t%s\n' "$alpha_seconds"
    printf 'simulation_total_seconds\t%s\n' "$sim_seconds"
    printf 'raw_seqwish_benchmark_seconds\t%s\n' "$raw_benchmark_seconds"
} > "$OUTPUT_DIR/run_metadata.tsv"

DUAL_SWEEP_LOG="${VS2_DUAL_SWEEP_LOG:-$ROOT/results/vallescope2_dual_option_trials.tsv}"
"$ENV_BIN/python" "$ROOT/scripts/append_vallescope2_dual_trial_log.py" \
    --trial-dir "$OUTPUT_DIR" \
    --log "$DUAL_SWEEP_LOG"
printf 'comparison_log\t%s\n' "$DUAL_SWEEP_LOG" >> "$OUTPUT_DIR/run_metadata.tsv"

printf 'complete\n' > "$status_file"
trap - ERR
echo
echo "== trial summary =="
column -t -s $'\t' "$OUTPUT_DIR/trial_summary.tsv"
echo
echo "complete: $OUTPUT_DIR"
echo "comparison log: $DUAL_SWEEP_LOG"
