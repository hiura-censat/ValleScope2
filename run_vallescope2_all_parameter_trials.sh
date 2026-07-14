#!/usr/bin/env bash

# ValleScope2 one-parameter-at-a-time parameter sweep
# Base commit / trial prefix: 46cacc3
#
# Each command changes only one parameter from the run_vallescope2_trial.sh
# baseline. The two tested values are selected as lower/stricter and
# higher/more-permissive settings where possible.
#
# This script intentionally does not use `set -e`, so a failed trial does not
# prevent subsequent trials from running.

# ---------------------------------------------------------------------------
# GenMap
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_kmer-length5 \
  --kmer-length 5
scripts/run_vallescope2_trial.sh 46cacc3_kmer-length8 \
  --kmer-length 8

scripts/run_vallescope2_trial.sh 46cacc3_kmer-errors1 \
  --kmer-errors 1
scripts/run_vallescope2_trial.sh 46cacc3_kmer-errors2 \
  --kmer-errors 2

# --threads changes computational resources rather than the alignment model,
# but it is included because it is listed as a value-bearing README option.
scripts/run_vallescope2_trial.sh 46cacc3_threads10 \
  --threads 10
scripts/run_vallescope2_trial.sh 46cacc3_threads40 \
  --threads 40

# ---------------------------------------------------------------------------
# Anchor detection
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_anchor-length35 \
  --anchor-length 35
scripts/run_vallescope2_trial.sh 46cacc3_anchor-length75 \
  --anchor-length 75

scripts/run_vallescope2_trial.sh 46cacc3_min-center-distance25 \
  --min-center-distance 25
scripts/run_vallescope2_trial.sh 46cacc3_min-center-distance100 \
  --min-center-distance 100

scripts/run_vallescope2_trial.sh 46cacc3_target-density8000 \
  --target-density 8000
scripts/run_vallescope2_trial.sh 46cacc3_target-density30000 \
  --target-density 30000

# ---------------------------------------------------------------------------
# Structural context
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_distance-bin-size5 \
  --distance-bin-size 5
scripts/run_vallescope2_trial.sh 46cacc3_distance-bin-size25 \
  --distance-bin-size 25

scripts/run_vallescope2_trial.sh 46cacc3_context-radius-tokens25 \
  --context-radius-tokens 25
scripts/run_vallescope2_trial.sh 46cacc3_context-radius-tokens100 \
  --context-radius-tokens 100

# ---------------------------------------------------------------------------
# Correspondence assignment
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_beta-tolerance15 \
  --beta-tolerance 15
scripts/run_vallescope2_trial.sh 46cacc3_beta-tolerance60 \
  --beta-tolerance 60

scripts/run_vallescope2_trial.sh 46cacc3_min-shared-tmus1 \
  --min-shared-tmus 1
scripts/run_vallescope2_trial.sh 46cacc3_min-shared-tmus6 \
  --min-shared-tmus 6

scripts/run_vallescope2_trial.sh 46cacc3_min-candidate-score5 \
  --min-candidate-score 5
scripts/run_vallescope2_trial.sh 46cacc3_min-candidate-score20 \
  --min-candidate-score 20

scripts/run_vallescope2_trial.sh 46cacc3_primary-margin1 \
  --primary-margin 1
scripts/run_vallescope2_trial.sh 46cacc3_primary-margin15 \
  --primary-margin 15

scripts/run_vallescope2_trial.sh 46cacc3_pair-merge-union \
  --pair-merge-mode union
scripts/run_vallescope2_trial.sh 46cacc3_pair-merge-reciprocal \
  --pair-merge-mode reciprocal

# ---------------------------------------------------------------------------
# Chaining
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_chain-predecessors20 \
  --chain-predecessors 20
scripts/run_vallescope2_trial.sh 46cacc3_chain-predecessors100 \
  --chain-predecessors 100

scripts/run_vallescope2_trial.sh 46cacc3_max-chain-gap300 \
  --max-chain-gap 300
scripts/run_vallescope2_trial.sh 46cacc3_max-chain-gap1500 \
  --max-chain-gap 1500

scripts/run_vallescope2_trial.sh 46cacc3_chain-max-gap-ratio1p01 \
  --chain-max-gap-ratio 1.01
scripts/run_vallescope2_trial.sh 46cacc3_chain-max-gap-ratio1p20 \
  --chain-max-gap-ratio 1.20

scripts/run_vallescope2_trial.sh 46cacc3_gap-cost-absolute \
  --gap-cost-model absolute
scripts/run_vallescope2_trial.sh 46cacc3_gap-cost-relative \
  --gap-cost-model relative

scripts/run_vallescope2_trial.sh 46cacc3_gap-weight0p0005 \
  --gap-weight 0.0005
scripts/run_vallescope2_trial.sh 46cacc3_gap-weight0p01 \
  --gap-weight 0.01

scripts/run_vallescope2_trial.sh 46cacc3_min-chain-anchors10 \
  --min-chain-anchors 10
scripts/run_vallescope2_trial.sh 46cacc3_min-chain-anchors40 \
  --min-chain-anchors 40

scripts/run_vallescope2_trial.sh 46cacc3_min-chain-both-anchors0 \
  --min-chain-both-anchors 0
scripts/run_vallescope2_trial.sh 46cacc3_min-chain-both-anchors5 \
  --min-chain-both-anchors 5

scripts/run_vallescope2_trial.sh 46cacc3_min-chain-score10 \
  --min-chain-score 10
scripts/run_vallescope2_trial.sh 46cacc3_min-chain-score100 \
  --min-chain-score 100

scripts/run_vallescope2_trial.sh 46cacc3_chain-trim-overlap0p001 \
  --chain-trim-overlap 0.001
scripts/run_vallescope2_trial.sh 46cacc3_chain-trim-overlap0p10 \
  --chain-trim-overlap 0.10

scripts/run_vallescope2_trial.sh 46cacc3_chain-reciprocal-weight0p4 \
  --chain-reciprocal-weight 0.4
scripts/run_vallescope2_trial.sh 46cacc3_chain-reciprocal-weight1p0 \
  --chain-reciprocal-weight 1.0

scripts/run_vallescope2_trial.sh 46cacc3_chain-multimap-overlap0p25 \
  --chain-multimap-overlap 0.25
scripts/run_vallescope2_trial.sh 46cacc3_chain-multimap-overlap0p75 \
  --chain-multimap-overlap 0.75

# ---------------------------------------------------------------------------
# Refinement
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_refinement-window10000 \
  --refinement-window 10000
scripts/run_vallescope2_trial.sh 46cacc3_refinement-window100000 \
  --refinement-window 100000

scripts/run_vallescope2_trial.sh 46cacc3_refinement-min-chain-anchors3 \
  --refinement-min-chain-anchors 3
scripts/run_vallescope2_trial.sh 46cacc3_refinement-min-chain-anchors10 \
  --refinement-min-chain-anchors 10

# ---------------------------------------------------------------------------
# Base alignment switches
# These runs test pipeline behavior rather than a continuous numeric value.
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_base-align-on \
  --base-align
scripts/run_vallescope2_trial.sh 46cacc3_base-align-off \
  --no-base-align

scripts/run_vallescope2_trial.sh 46cacc3_chain-extension-on \
  --chain-extension
scripts/run_vallescope2_trial.sh 46cacc3_chain-extension-off \
  --no-chain-extension

# ---------------------------------------------------------------------------
# Chain extension
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_max-chain-extension-bp20000 \
  --max-chain-extension-bp 20000
scripts/run_vallescope2_trial.sh 46cacc3_max-chain-extension-bp150000 \
  --max-chain-extension-bp 150000

scripts/run_vallescope2_trial.sh 46cacc3_min-chain-extension-anchors1 \
  --min-chain-extension-anchors 1
scripts/run_vallescope2_trial.sh 46cacc3_min-chain-extension-anchors5 \
  --min-chain-extension-anchors 5

scripts/run_vallescope2_trial.sh 46cacc3_min-chain-extension-score10 \
  --min-chain-extension-score 10
scripts/run_vallescope2_trial.sh 46cacc3_min-chain-extension-score100 \
  --min-chain-extension-score 100

scripts/run_vallescope2_trial.sh 46cacc3_min-copy-support-anchors1 \
  --min-copy-support-anchors 1
scripts/run_vallescope2_trial.sh 46cacc3_min-copy-support-anchors5 \
  --min-copy-support-anchors 5

# ---------------------------------------------------------------------------
# Bundle / patch alignment
# ---------------------------------------------------------------------------

scripts/run_vallescope2_trial.sh 46cacc3_max-bundle-align-bp200000 \
  --max-bundle-align-bp 200000
scripts/run_vallescope2_trial.sh 46cacc3_max-bundle-align-bp2000000 \
  --max-bundle-align-bp 2000000

scripts/run_vallescope2_trial.sh 46cacc3_max-patch-gap-bp20000 \
  --max-patch-gap-bp 20000
scripts/run_vallescope2_trial.sh 46cacc3_max-patch-gap-bp150000 \
  --max-patch-gap-bp 150000

scripts/run_vallescope2_trial.sh 46cacc3_patch-flank-bp200 \
  --patch-flank-bp 200
scripts/run_vallescope2_trial.sh 46cacc3_patch-flank-bp2000 \
  --patch-flank-bp 2000

# --max-wfa-memory-gb changes a resource ceiling rather than alignment scoring.
scripts/run_vallescope2_trial.sh 46cacc3_max-wfa-memory-gb16 \
  --max-wfa-memory-gb 16
scripts/run_vallescope2_trial.sh 46cacc3_max-wfa-memory-gb128 \
  --max-wfa-memory-gb 128
