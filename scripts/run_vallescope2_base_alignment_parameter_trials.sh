#!/usr/bin/env bash
# One-parameter-at-a-time sweep for Base-level bundle alignment options.
#
# Usage:
#   scripts/run_vallescope2_base_alignment_parameter_trials.sh [trial_prefix]
#
# Environment:
#   VS2_SWEEP_DRY_RUN=1       Print commands without running them.
#   VS2_DUAL_TRIAL_THREADS=N  Threads passed through to each dual trial.
#
# Completed output directories are skipped, so rerunning the same prefix
# resumes the sweep. Existing incomplete directories are reported and skipped
# because run_vallescope2_dual_trial.sh intentionally never overwrites output.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PREFIX="${1:-base_alignment_sweep_$(date +%Y%m%d)}"
DRY_RUN="${VS2_SWEEP_DRY_RUN:-0}"
completed=0
skipped=0
failed=0

run_trial() {
    local suffix="$1"
    shift
    local name="${PREFIX}_${suffix}"
    local output="results/parameter_trials/${name}"
    local command=(scripts/run_vallescope2_dual_trial.sh "$name" "$@")

    if [[ -f "$output/STATUS" ]] && [[ "$(<"$output/STATUS")" == "complete" ]]; then
        echo "SKIP complete: $name"
        ((skipped += 1))
        return
    fi
    if [[ -e "$output" ]]; then
        echo "SKIP incomplete output (remove or rename it to retry): $output" >&2
        ((failed += 1))
        return
    fi

    printf 'RUN:'
    printf ' %q' "${command[@]}"
    printf '\n'
    if [[ "$DRY_RUN" == "1" ]]; then
        ((skipped += 1))
        return
    fi
    if "${command[@]}"; then
        ((completed += 1))
    else
        echo "FAILED: $name" >&2
        ((failed += 1))
    fi
}

# Chain extension (default: on, 70000 bp, 5 anchors, score 100, copy support 1)
run_trial no-chain-extension --no-chain-extension
run_trial max-chain-extension-bp20000 --max-chain-extension-bp 20000
run_trial max-chain-extension-bp150000 --max-chain-extension-bp 150000
run_trial min-chain-extension-anchors2 --min-chain-extension-anchors 2
run_trial min-chain-extension-anchors10 --min-chain-extension-anchors 10
run_trial min-chain-extension-score25 --min-chain-extension-score 25
run_trial min-chain-extension-score200 --min-chain-extension-score 200
run_trial min-copy-support-anchors0 --min-copy-support-anchors 0
run_trial min-copy-support-anchors5 --min-copy-support-anchors 5

# Bundle and patch interval bounds
run_trial max-bundle-align-bp200000 --max-bundle-align-bp 200000
run_trial max-bundle-align-bp2000000 --max-bundle-align-bp 2000000
run_trial max-patch-gap-bp20000 --max-patch-gap-bp 20000
run_trial max-patch-gap-bp150000 --max-patch-gap-bp 150000
run_trial patch-flank-bp200 --patch-flank-bp 200
run_trial patch-flank-bp2000 --patch-flank-bp 2000

# Direct patch and long-I/D rescue quality
run_trial min-patch-identity0p90 --min-patch-identity 0.90
run_trial min-patch-identity0p98 --min-patch-identity 0.98
run_trial min-patch-long-indel-bp250 --min-patch-long-indel-bp 250
run_trial min-patch-long-indel-bp1000 --min-patch-long-indel-bp 1000
run_trial min-patch-rescue-flank-bp100 --min-patch-rescue-flank-bp 100
run_trial min-patch-rescue-flank-bp500 --min-patch-rescue-flank-bp 500
run_trial min-patch-rescue-flank-identity0p90 --min-patch-rescue-flank-identity 0.90
run_trial min-patch-rescue-flank-identity0p99 --min-patch-rescue-flank-identity 0.99
run_trial max-patch-rescue-extra-indel-bp0 --max-patch-rescue-extra-indel-bp 0
run_trial max-patch-rescue-extra-indel-bp500 --max-patch-rescue-extra-indel-bp 500
run_trial patch-quality-window-bp200 --patch-quality-window-bp 200
run_trial patch-quality-window-bp1000 --patch-quality-window-bp 1000
run_trial min-patch-endpoint-identity0p75 --min-patch-endpoint-identity 0.75
run_trial min-patch-endpoint-identity0p95 --min-patch-endpoint-identity 0.95
run_trial max-patch-short-error-density0p05 --max-patch-short-error-density 0.05
run_trial max-patch-short-error-density0p25 --max-patch-short-error-density 0.25

# Clean multi-event rescue (defaults: I-only, segment identity 1.0, no errors)
run_trial patch-multi-event-allow-deletions --patch-multi-event-allow-deletions
run_trial min-patch-multi-segment-identity0p95 --min-patch-multi-segment-identity 0.95
run_trial min-patch-multi-segment-identity0p99 --min-patch-multi-segment-identity 0.99
run_trial max-patch-multi-short-indel-bp100 --max-patch-multi-short-indel-bp 100
run_trial max-patch-multi-short-indel-bp500 --max-patch-multi-short-indel-bp 500
run_trial max-patch-multi-short-error-density0p05 --max-patch-multi-short-error-density 0.05
run_trial max-patch-multi-short-error-density0p15 --max-patch-multi-short-error-density 0.15

# WFA resource ceiling (default: 64 GB)
run_trial max-wfa-memory-gb16 --max-wfa-memory-gb 16
run_trial max-wfa-memory-gb128 --max-wfa-memory-gb 128

echo
echo "Sweep finished: completed=$completed skipped=$skipped failed=$failed"
echo "Shared comparison: results/vallescope2_dual_option_trials.tsv"
((failed == 0))
