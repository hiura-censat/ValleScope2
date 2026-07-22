#!/usr/bin/env bash
# ValleScope2 parameter sweep focused on reducing unresolved reference/query gaps.
#
# Usage:
#   scripts/run_vallescope2_gap_reduction_parameter_trials.sh [trial_prefix]
#
# Environment:
#   VS2_SWEEP_DRY_RUN=1       Print commands without running them.
#   VS2_DUAL_TRIAL_THREADS=N  Threads passed through to each dual trial.
#
# Most trials change one parameter from the current defaults in a permissive
# direction. The final profiles combine related settings and are labeled
# separately. Completed outputs are skipped; incomplete outputs are not
# overwritten. One failed trial does not stop the remaining sweep.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PREFIX="${1:-gap_reduction_sweep_$(date +%Y%m%d)}"
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

# Current patch-related defaults used as the baseline:
#   chain extension: 1000000 bp, 5 anchors, score 100
#   max patch gap: 1000000 bp; patch flank: 500 bp
#   direct identity: 0.90; long I/D: 100 bp
#   single rescue: flank 200 bp / identity 0.95 / extra I/D 100 bp
#   endpoint identity: 0.85; short-error density: 0.15
#   multi rescue: segment identity 0.99 / short I/D 100 bp / density 0.15
#   multi-event deletions: on; patch Z-drop: off

# ---------------------------------------------------------------------------
# One-parameter-at-a-time trials
# ---------------------------------------------------------------------------

# Recover more partial extension paths before WFA patching.
run_trial max-chain-extension-bp1500000 --max-chain-extension-bp 1500000
run_trial max-chain-extension-bp2000000 --max-chain-extension-bp 2000000
run_trial min-chain-extension-anchors3 --min-chain-extension-anchors 3
run_trial min-chain-extension-anchors1 --min-chain-extension-anchors 1
run_trial min-chain-extension-score50 --min-chain-extension-score 50
run_trial min-chain-extension-score25 --min-chain-extension-score 25

# Attempt larger unresolved intervals and permit larger final bundle alignment.
run_trial max-patch-gap-bp1500000 --max-patch-gap-bp 1500000
run_trial max-patch-gap-bp2000000 --max-patch-gap-bp 2000000
run_trial max-bundle-align-bp1500000 --max-bundle-align-bp 1500000
run_trial max-bundle-align-bp2000000 --max-bundle-align-bp 2000000

# Accept moderately lower-identity end-to-end patches.
run_trial min-patch-identity0p85 --min-patch-identity 0.85
run_trial min-patch-identity0p80 --min-patch-identity 0.80

# Recognize smaller structural I/D operations as long-I/D rescue candidates.
run_trial min-patch-long-indel-bp50 --min-patch-long-indel-bp 50
run_trial min-patch-long-indel-bp25 --min-patch-long-indel-bp 25

# Relax single-long-I/D flank support in two steps.
run_trial min-patch-rescue-flank-bp100 --min-patch-rescue-flank-bp 100
run_trial min-patch-rescue-flank-bp50 --min-patch-rescue-flank-bp 50
run_trial min-patch-rescue-flank-identity0p90 --min-patch-rescue-flank-identity 0.90
run_trial min-patch-rescue-flank-identity0p85 --min-patch-rescue-flank-identity 0.85
run_trial max-patch-rescue-extra-indel-bp250 --max-patch-rescue-extra-indel-bp 250
run_trial max-patch-rescue-extra-indel-bp500 --max-patch-rescue-extra-indel-bp 500

# Relax endpoint and local short-error checks used by rescue classification.
run_trial min-patch-endpoint-identity0p80 --min-patch-endpoint-identity 0.80
run_trial min-patch-endpoint-identity0p75 --min-patch-endpoint-identity 0.75
run_trial max-patch-short-error-density0p20 --max-patch-short-error-density 0.20
run_trial max-patch-short-error-density0p25 --max-patch-short-error-density 0.25

# Recover clean or nearly clean multiple-long-I/D intervals.
run_trial min-patch-multi-segment-identity0p98 --min-patch-multi-segment-identity 0.98
run_trial min-patch-multi-segment-identity0p95 --min-patch-multi-segment-identity 0.95
run_trial max-patch-multi-short-indel-bp250 --max-patch-multi-short-indel-bp 250
run_trial max-patch-multi-short-indel-bp500 --max-patch-multi-short-indel-bp 500
run_trial max-patch-multi-short-error-density0p20 --max-patch-multi-short-error-density 0.20
run_trial max-patch-multi-short-error-density0p25 --max-patch-multi-short-error-density 0.25

# ---------------------------------------------------------------------------
# Combined profiles
# These are intentionally separate from OAT trials because their effects
# cannot be attributed to one parameter.
# ---------------------------------------------------------------------------

run_trial profile-balanced \
    --max-patch-gap-bp 1500000 \
    --min-patch-identity 0.85 \
    --min-patch-rescue-flank-bp 100 \
    --min-patch-rescue-flank-identity 0.90 \
    --max-patch-rescue-extra-indel-bp 250 \
    --min-patch-multi-segment-identity 0.98 \
    --max-patch-multi-short-indel-bp 250 \
    --max-patch-multi-short-error-density 0.20

run_trial profile-permissive \
    --max-chain-extension-bp 2000000 \
    --min-chain-extension-anchors 3 \
    --min-chain-extension-score 50 \
    --max-patch-gap-bp 2000000 \
    --max-bundle-align-bp 2000000 \
    --min-patch-identity 0.80 \
    --min-patch-long-indel-bp 50 \
    --min-patch-rescue-flank-bp 50 \
    --min-patch-rescue-flank-identity 0.85 \
    --max-patch-rescue-extra-indel-bp 500 \
    --min-patch-endpoint-identity 0.75 \
    --max-patch-short-error-density 0.25 \
    --min-patch-multi-segment-identity 0.95 \
    --max-patch-multi-short-indel-bp 500 \
    --max-patch-multi-short-error-density 0.25

echo
echo "Sweep finished: completed=$completed skipped=$skipped failed=$failed"
echo "Shared comparison: results/vallescope2_dual_option_trials.tsv"
((failed == 0))
