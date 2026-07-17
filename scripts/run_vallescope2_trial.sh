#!/usr/bin/env bash
#
# Run one full vallescope2 -> pggb -> truth benchmark trial and report P/R/F1
# and wall-clock timing. Does not modify vallescope2 itself; only orchestrates
# the existing build/vallescope2 binary via `vsg-graph construct` and
# `vsg-benchmark`, exactly as done manually throughout this project's
# no_translocation experiments.
#
# Usage:
#   scripts/run_vallescope2_trial.sh <run_name> [vallescope2 options...]
#
# Example:
#   scripts/run_vallescope2_trial.sh mytrial01 --min-chain-anchors 15 --max-patch-gap-bp 100000
#   scripts/run_vallescope2_trial.sh mytrial02 --no-chain-extension
#
# <run_name> becomes the run/results directory name (suffixed with the
# dataset tag). Re-running with the same name overwrites that trial's output.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <run_name> [vallescope2 options...]" >&2
    exit 1
fi

RUN_NAME_BASE="$1"
shift
VALLESCOPE2_OPTS=("$@")

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

ENV_BIN="${VS2_TRIAL_ENV_BIN:-/home/senescence/miniconda3/envs/vallescope_dev/bin}"
DATA="${VS2_TRIAL_SIM_DATA:-data/test/chr1_alphaSat_1M_no_translocation}"
THREADS="${VS2_TRIAL_THREADS:-8}"

RUN_NAME="${VS2_TRIAL_EXACT_RUN_NAME:-${RUN_NAME_BASE}_no_translocation_$(date +%Y%m%d)}"
GRAPH_DIR="${VS2_TRIAL_GRAPH_DIR:-results/graphs}"
RUN_DIR="${VS2_TRIAL_RUN_DIR:-runs}"
RESULT_DIR="${VS2_TRIAL_BENCHMARK_DIR:-results/simulation/${RUN_NAME}_type_harmonized}"
GFA="${GRAPH_DIR}/${RUN_NAME}.vallescope2.gfa"
LOG_DIR="${VS2_TRIAL_LOG_DIR:-/tmp/claude-1001}"
mkdir -p "$LOG_DIR" "$GRAPH_DIR" "$RUN_DIR" "$(dirname "$RESULT_DIR")"

echo "== vallescope2 trial: ${RUN_NAME} =="
echo "options: ${VALLESCOPE2_OPTS[*]:-(none, defaults)}"
echo

# --- Step 1: vallescope2 + pggb via vsg-graph construct -------------------
CONSTRUCT_LOG="${LOG_DIR}/${RUN_NAME}.construct.log"
construct_start=$(date +%s)
PATH="$ROOT/build:$ENV_BIN:$PATH" "$ENV_BIN/vsg-graph" construct \
    -f "$DATA/input_fastas_renamed/reference.fa" \
    -f "$DATA/input_fastas_renamed/sample01.fa" \
    -f "$DATA/input_fastas_renamed/sample02.fa" \
    -f "$DATA/input_fastas_renamed/sample03.fa" \
    --name "$RUN_NAME" \
    --output "$GRAPH_DIR" \
    --run-directory "$RUN_DIR" \
    --aligner vallescope2 \
    --threads "$THREADS" \
    --aligner-extra "${VALLESCOPE2_OPTS[*]:-}" \
    > "$CONSTRUCT_LOG" 2>&1
construct_status=$?
construct_end=$(date +%s)
construct_seconds=$((construct_end - construct_start))

if [[ $construct_status -ne 0 ]]; then
    echo "FAILED: vsg-graph construct exited $construct_status. See $CONSTRUCT_LOG" >&2
    tail -30 "$CONSTRUCT_LOG" >&2
    exit 1
fi
if [[ ! -s "$GFA" ]]; then
    echo "FAILED: expected GFA not found at $GFA. See $CONSTRUCT_LOG" >&2
    exit 1
fi
echo "construct: ok (${construct_seconds}s)  log: $CONSTRUCT_LOG"

# --- Step 2: truth benchmark ------------------------------------------------
BENCH_LOG="${LOG_DIR}/${RUN_NAME}.bench.log"
bench_start=$(date +%s)
"$ENV_BIN/vsg-benchmark" simulate benchmark-graph-vcf \
    --gfa "$GFA" \
    --reference-path reference \
    --reference-fasta "$DATA/truth_renamed/reference/visor_reference.fa" \
    --truth-vcf "$DATA/truth_renamed/vcfs/sample01.truth.vcf" \
    --truth-vcf "$DATA/truth_renamed/vcfs/sample02.truth.vcf" \
    --truth-vcf "$DATA/truth_renamed/vcfs/sample03.truth.vcf" \
    --output "$RESULT_DIR" \
    --all-snarls --refdist 500 --pctseq 0.7 --pctsize 0.7 \
    --type-ignore --dup-to-ins --threads "$THREADS" \
    --vg "$ENV_BIN/vg" --bcftools "$ENV_BIN/bcftools" \
    --bgzip "$ENV_BIN/bgzip" --tabix "$ENV_BIN/tabix" \
    --truvari "$ENV_BIN/truvari" \
    > "$BENCH_LOG" 2>&1
bench_status=$?
bench_end=$(date +%s)
bench_seconds=$((bench_end - bench_start))

if [[ $bench_status -ne 0 ]]; then
    echo "FAILED: vsg-benchmark exited $bench_status. See $BENCH_LOG" >&2
    tail -30 "$BENCH_LOG" >&2
    exit 1
fi
echo "benchmark: ok (${bench_seconds}s)  log: $BENCH_LOG"
echo

# --- Step 3: report ----------------------------------------------------------
SUMMARY="${RESULT_DIR}/benchmark_summary.tsv"
echo "== per-sample P/R/F1 =="
column -t "$SUMMARY"
echo

aggregate=$(awk -F'\t' '
    NR>1 { tp+=$3; fp+=$5; fn+=$6 }
    END {
        p = (tp+fp>0) ? tp/(tp+fp) : 0
        r = (tp+fn>0) ? tp/(tp+fn) : 0
        f1 = (p+r>0) ? 2*p*r/(p+r) : 0
        printf "TP=%d FP=%d FN=%d  P=%.4f R=%.4f F1=%.4f", tp, fp, fn, p, r, f1
    }' "$SUMMARY")
echo "== aggregate: ${aggregate} =="
echo "== timing: construct=${construct_seconds}s benchmark=${bench_seconds}s total=$((construct_seconds + bench_seconds))s =="

# --- Step 4: append to running comparison log --------------------------------
SWEEP_LOG="${VS2_TRIAL_SWEEP_LOG:-results/vallescope2_option_trials.tsv}"
mkdir -p "$(dirname "$SWEEP_LOG")"
if [[ ! -f "$SWEEP_LOG" ]]; then
    printf 'timestamp\trun_name\toptions\tconstruct_s\tbench_s\tTP\tFP\tFN\tP\tR\tF1\n' > "$SWEEP_LOG"
fi
read -r tp fp fn p r f1 <<< "$(echo "$aggregate" | sed -E 's/TP=([0-9]+) FP=([0-9]+) FN=([0-9]+)  P=([0-9.]+) R=([0-9.]+) F1=([0-9.]+)/\1 \2 \3 \4 \5 \6/')"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -Iseconds)" "$RUN_NAME" "${VALLESCOPE2_OPTS[*]:-(default)}" \
    "$construct_seconds" "$bench_seconds" "$tp" "$fp" "$fn" "$p" "$r" "$f1" >> "$SWEEP_LOG"
echo
echo "logged to: $SWEEP_LOG"
