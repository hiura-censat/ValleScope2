#!/usr/bin/env bash

# ValleScope2 one-parameter-at-a-time parameter sweep
# Base commit / trial prefix: 03be8f7
#
# Each trial changes only one option from the current ValleScope2 defaults.
# Two non-default values are tested for each numeric/categorical parameter.
#
# No `set -e`: one failed trial will not stop the remaining trials.

# ===========================================================================
# GenMap
#
# Current defaults:
#   --kmer-length 6
#   --kmer-errors 0
#   --threads 20
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_kmer-length5 \
  --kmer-length 5

scripts/run_vallescope2_trial.sh 03be8f7_kmer-length8 \
  --kmer-length 8


scripts/run_vallescope2_trial.sh 03be8f7_kmer-errors1 \
  --kmer-errors 1

scripts/run_vallescope2_trial.sh 03be8f7_kmer-errors2 \
  --kmer-errors 2


# Resource parameter; mainly useful for checking determinism and runtime.
scripts/run_vallescope2_trial.sh 03be8f7_threads10 \
  --threads 10

scripts/run_vallescope2_trial.sh 03be8f7_threads40 \
  --threads 40


# ===========================================================================
# Anchor detection
#
# Current defaults:
#   --anchor-length 50
#   --min-center-distance 50
#   --target-density 15000
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_anchor-length35 \
  --anchor-length 35

scripts/run_vallescope2_trial.sh 03be8f7_anchor-length75 \
  --anchor-length 75


scripts/run_vallescope2_trial.sh 03be8f7_min-center-distance25 \
  --min-center-distance 25

scripts/run_vallescope2_trial.sh 03be8f7_min-center-distance100 \
  --min-center-distance 100


scripts/run_vallescope2_trial.sh 03be8f7_target-density8000 \
  --target-density 8000

scripts/run_vallescope2_trial.sh 03be8f7_target-density30000 \
  --target-density 30000


# ===========================================================================
# Structural context
#
# Current defaults:
#   --distance-bin-size 10
#   --context-radius-tokens 50
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_distance-bin-size5 \
  --distance-bin-size 5

scripts/run_vallescope2_trial.sh 03be8f7_distance-bin-size25 \
  --distance-bin-size 25


scripts/run_vallescope2_trial.sh 03be8f7_context-radius-tokens25 \
  --context-radius-tokens 25

scripts/run_vallescope2_trial.sh 03be8f7_context-radius-tokens100 \
  --context-radius-tokens 100


# ===========================================================================
# Correspondence assignment
#
# Current defaults:
#   --beta-tolerance 30
#   --min-shared-tmus 3
#   --min-candidate-score 10
#   --primary-margin 5
#   --pair-merge-mode union
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_beta-tolerance15 \
  --beta-tolerance 15

scripts/run_vallescope2_trial.sh 03be8f7_beta-tolerance60 \
  --beta-tolerance 60


scripts/run_vallescope2_trial.sh 03be8f7_min-shared-tmus1 \
  --min-shared-tmus 1

scripts/run_vallescope2_trial.sh 03be8f7_min-shared-tmus6 \
  --min-shared-tmus 6


scripts/run_vallescope2_trial.sh 03be8f7_min-candidate-score5 \
  --min-candidate-score 5

scripts/run_vallescope2_trial.sh 03be8f7_min-candidate-score20 \
  --min-candidate-score 20


scripts/run_vallescope2_trial.sh 03be8f7_primary-margin1 \
  --primary-margin 1

scripts/run_vallescope2_trial.sh 03be8f7_primary-margin15 \
  --primary-margin 15


# union is currently the default. The explicit union run is useful for
# confirming that the trial wrapper reproduces the default result.
scripts/run_vallescope2_trial.sh 03be8f7_pair-merge-union \
  --pair-merge-mode union

scripts/run_vallescope2_trial.sh 03be8f7_pair-merge-reciprocal \
  --pair-merge-mode reciprocal


# ===========================================================================
# Chaining
#
# Current defaults:
#   --chain-predecessors 50
#   --max-chain-gap 300
#   --chain-max-gap-ratio 1.05
#   --gap-cost-model absolute
#   --gap-weight 0.002
#   --min-chain-anchors 10
#   --min-chain-both-anchors 1
#   --min-chain-score 0
#   --chain-trim-overlap 0.01
#   --chain-reciprocal-weight 0.8
#   --chain-multimap-overlap 0.5
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_chain-predecessors20 \
  --chain-predecessors 20

scripts/run_vallescope2_trial.sh 03be8f7_chain-predecessors100 \
  --chain-predecessors 100


# Current default is 300, so test a stricter and a more permissive value.
scripts/run_vallescope2_trial.sh 03be8f7_max-chain-gap150 \
  --max-chain-gap 150

scripts/run_vallescope2_trial.sh 03be8f7_max-chain-gap750 \
  --max-chain-gap 750


scripts/run_vallescope2_trial.sh 03be8f7_chain-max-gap-ratio1p01 \
  --chain-max-gap-ratio 1.01

scripts/run_vallescope2_trial.sh 03be8f7_chain-max-gap-ratio1p20 \
  --chain-max-gap-ratio 1.20


# absolute is currently the default. This comparison has only two modes.
scripts/run_vallescope2_trial.sh 03be8f7_gap-cost-absolute \
  --gap-cost-model absolute

scripts/run_vallescope2_trial.sh 03be8f7_gap-cost-relative \
  --gap-cost-model relative


scripts/run_vallescope2_trial.sh 03be8f7_gap-weight0p0005 \
  --gap-weight 0.0005

scripts/run_vallescope2_trial.sh 03be8f7_gap-weight0p01 \
  --gap-weight 0.01


# Current default is 10.
scripts/run_vallescope2_trial.sh 03be8f7_min-chain-anchors5 \
  --min-chain-anchors 5

scripts/run_vallescope2_trial.sh 03be8f7_min-chain-anchors20 \
  --min-chain-anchors 20


scripts/run_vallescope2_trial.sh 03be8f7_min-chain-both-anchors0 \
  --min-chain-both-anchors 0

scripts/run_vallescope2_trial.sh 03be8f7_min-chain-both-anchors5 \
  --min-chain-both-anchors 5


# Current default is 0.
scripts/run_vallescope2_trial.sh 03be8f7_min-chain-score10 \
  --min-chain-score 10

scripts/run_vallescope2_trial.sh 03be8f7_min-chain-score100 \
  --min-chain-score 100


scripts/run_vallescope2_trial.sh 03be8f7_chain-trim-overlap0p001 \
  --chain-trim-overlap 0.001

scripts/run_vallescope2_trial.sh 03be8f7_chain-trim-overlap0p10 \
  --chain-trim-overlap 0.10


scripts/run_vallescope2_trial.sh 03be8f7_chain-reciprocal-weight0p4 \
  --chain-reciprocal-weight 0.4

scripts/run_vallescope2_trial.sh 03be8f7_chain-reciprocal-weight1p0 \
  --chain-reciprocal-weight 1.0


scripts/run_vallescope2_trial.sh 03be8f7_chain-multimap-overlap0p25 \
  --chain-multimap-overlap 0.25

scripts/run_vallescope2_trial.sh 03be8f7_chain-multimap-overlap0p75 \
  --chain-multimap-overlap 0.75


# ===========================================================================
# Refinement
#
# Current defaults:
#   --refinement-window 100000
#   --refinement-min-chain-anchors 20
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_refinement-window20000 \
  --refinement-window 20000

scripts/run_vallescope2_trial.sh 03be8f7_refinement-window200000 \
  --refinement-window 200000


scripts/run_vallescope2_trial.sh 03be8f7_refinement-min-chain-anchors5 \
  --refinement-min-chain-anchors 5

scripts/run_vallescope2_trial.sh 03be8f7_refinement-min-chain-anchors40 \
  --refinement-min-chain-anchors 40


# ===========================================================================
# Base-level bundle alignment switches
#
# Current defaults:
#   base alignment: on
#   refined chains: on
#   chain extension: on
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_base-align-on \
  --base-align

scripts/run_vallescope2_trial.sh 03be8f7_base-align-off \
  --no-base-align

scripts/run_vallescope2_trial.sh 03be8f7_chain-extension-on \
  --chain-extension

scripts/run_vallescope2_trial.sh 03be8f7_chain-extension-off \
  --no-chain-extension


# ===========================================================================
# Chain extension
#
# Current defaults:
#   --max-chain-extension-bp 70000
#   --min-chain-extension-anchors 5
#   --min-chain-extension-score 100
#   --min-copy-support-anchors 1
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_max-chain-extension-bp20000 \
  --max-chain-extension-bp 20000

scripts/run_vallescope2_trial.sh 03be8f7_max-chain-extension-bp150000 \
  --max-chain-extension-bp 150000


scripts/run_vallescope2_trial.sh 03be8f7_min-chain-extension-anchors2 \
  --min-chain-extension-anchors 2

scripts/run_vallescope2_trial.sh 03be8f7_min-chain-extension-anchors10 \
  --min-chain-extension-anchors 10


scripts/run_vallescope2_trial.sh 03be8f7_min-chain-extension-score25 \
  --min-chain-extension-score 25

scripts/run_vallescope2_trial.sh 03be8f7_min-chain-extension-score200 \
  --min-chain-extension-score 200


# Current default is 1.
scripts/run_vallescope2_trial.sh 03be8f7_min-copy-support-anchors0 \
  --min-copy-support-anchors 0

scripts/run_vallescope2_trial.sh 03be8f7_min-copy-support-anchors5 \
  --min-copy-support-anchors 5


# ===========================================================================
# Bundle / patch alignment
#
# Current defaults:
#   --max-bundle-align-bp 1000000
#   --max-patch-gap-bp 70000
#   --patch-flank-bp 500
#   --max-wfa-memory-gb 64
# ===========================================================================

scripts/run_vallescope2_trial.sh 03be8f7_max-bundle-align-bp200000 \
  --max-bundle-align-bp 200000

scripts/run_vallescope2_trial.sh 03be8f7_max-bundle-align-bp2000000 \
  --max-bundle-align-bp 2000000


scripts/run_vallescope2_trial.sh 03be8f7_max-patch-gap-bp20000 \
  --max-patch-gap-bp 20000

scripts/run_vallescope2_trial.sh 03be8f7_max-patch-gap-bp150000 \
  --max-patch-gap-bp 150000


scripts/run_vallescope2_trial.sh 03be8f7_patch-flank-bp200 \
  --patch-flank-bp 200

scripts/run_vallescope2_trial.sh 03be8f7_patch-flank-bp2000 \
  --patch-flank-bp 2000


# Resource ceiling rather than an alignment-scoring parameter.
scripts/run_vallescope2_trial.sh 03be8f7_max-wfa-memory-gb16 \
  --max-wfa-memory-gb 16

scripts/run_vallescope2_trial.sh 03be8f7_max-wfa-memory-gb128 \
  --max-wfa-memory-gb 128
