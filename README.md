# ValleScope2

ValleScope2 is an experimental structural-context-based genome alignment tool.
The current implementation builds low-copy anchor candidates from GenMap k-mer
profiles, constructs token-level structural contexts, assigns corresponding
anchors across samples, chains those correspondences into bundles, optionally
refines local chains, and can run WFA2-lib global alignment for first-pass
bundle intervals.

This repository is still under active development. Output formats and scoring
details may change.

## Dependencies

Build-time dependencies:

- CMake >= 3.16
- C++17 compiler
- OpenSSL development files
- HTSlib development files
- WFA2-lib submodule, included at `external/WFA2-lib`

Runtime dependency:

- GenMap executable

ValleScope2 looks for GenMap on `PATH`, or at `.tools/genmap/bin/genmap`.
You can also specify it explicitly with `--genmap PATH`.

## Clone

Use `--recursive` so that WFA2-lib is checked out.

```bash
git clone --recursive <repo-url>
cd ValleScope2_v0.0.0
```

If the repository was already cloned without submodules:

```bash
git submodule update --init --recursive
```

## Build

If HTSlib is installed in a standard location:

```bash
cmake -S . -B build
cmake --build build -j
```

If HTSlib is installed in a custom prefix:

```bash
cmake -S . -B build -DHTSLIB_ROOT=/path/to/htslib
cmake --build build -j
```

Using the local wfmash/HTSlib environment used during development:

```bash
cmake -S . -B build -DHTSLIB_ROOT="$PWD/.tools/wfmash"
cmake --build build -j
```

WFA2-lib is built automatically from the submodule by CMake.

## Basic commands

One FASTA: all-vs-all self comparison.

```bash
build/vallescope2 genome.fa
```

Two FASTAs: target-query comparison.

```bash
build/vallescope2 target.fa query.fa > aln.paf
```

Two or more FASTAs combined into one normalized multi-FASTA, followed by
all-vs-all self comparison.

```bash
build/vallescope2 --combine hap1.fa hap2.fa > aln.paf
build/vallescope2 --combine hap1.fa hap2.fa hap3.fa > aln.paf
```

Keep intermediate files:

```bash
build/vallescope2 --combine --debug hap1.fa hap2.fa hap3.fa > aln.paf
```

Base-level WFA2 alignment is enabled by default. To stop after chaining during
development:

```bash
build/vallescope2 --no-base-align --combine --debug hap1.fa hap2.fa hap3.fa
```

Increase the maximum bundle interval size allowed for WFA2 alignment:

```bash
build/vallescope2 --max-bundle-align-bp 200000 target.fa query.fa > aln.paf
```

## Important options

Input and debugging:

```text
--combine
--debug
--dump-window-scores
--genmap PATH
```

GenMap:

```text
-K, --kmer-length INT       default: 6
-E, --kmer-errors INT      default: 0
-t, --threads INT          default: 20
```

Anchor detection:

```text
-al, --anchor-length INT           default: 50
-md, --min-center-distance INT     default: 50
-td, --target-density INT          default: 15000 anchors/Mb
```

Structural context:

```text
-db, --distance-bin-size INT       default: 10
-cr, --context-radius-tokens INT   default: 50
```

Correspondence assignment:

```text
-bt, --beta-tolerance INT          default: 30
-mcs, --min-candidate-score INT    default: -10
-pm, --primary-margin INT          default: 5
--pair-merge-mode MODE             union or reciprocal, default: union
```

Chaining and refinement:

```text
--chain-predecessors INT           default: 50
--max-chain-gap INT                default: 750
--chain-max-gap-ratio FLOAT        default: 1.2
--gap-cost-model MODE              default: absolute
--gap-weight FLOAT                 default: 0.002
--min-chain-anchors INT            default: 20
--min-chain-both-anchors INT       default: 0
--min-chain-score FLOAT            default: 0
--chain-trim-overlap FLOAT         default: 0.01
--refinement-window INT            default: 50000
--refinement-min-chain-anchors INT default: 5
```

Base-level bundle alignment:

```text
--base-align
--no-base-align
--max-bundle-align-bp INT          default: 50000
--max-patch-gap-bp INT            default: 0
--patch-window-bp INT             default: 1000
--min-patch-identity FLOAT        default: 0.85
--max-wfa-memory-gb INT           default: 64
```

## Current workflow

### 1. Input preparation

ValleScope2 supports three input modes:

- one FASTA: self comparison
- two FASTAs: target-query comparison
- `--combine` with two or more FASTAs: combine inputs, then self comparison

ValleScope2 uses the first token of each FASTA header as the sequence ID,
matching the convention used by common aligners. For example:

```text
>chr1 description...
```

is used as `chr1`. Sequence IDs must be unique across combined inputs; ValleScope2
does not automatically rename duplicated FASTA headers.

With `--debug`, the prepared FASTA, `.fai`, and `sequences.tsv` are kept in
`vallescope2-debug/`.

### 2. GenMap profile

ValleScope2 runs:

```text
genmap index
genmap map -K <k> -E <e> -bg -fl -T <threads>
```

The default is `-K 6 -E 0`.

For target-query mode, target and query FASTAs are combined for GenMap so the
k-mer profile is computed over both sides together.

### 3. Anchor detection

For each sequence, ValleScope2 computes a sliding-window score:

```text
score = mean(log(1 + kmer_frequency))
```

The window width is `--anchor-length`.

Windows are processed from low score to high score. Anchors within
`--min-center-distance` of an already accepted anchor are removed. Selection
stops at approximately `--target-density` anchors/Mb.

Persistent outputs:

```text
anchors.bed
anchors.meta.json
genmap.log
```

If `--dump-window-scores` is specified:

```text
window_scores.tsv
```

### 4. Anchor grouping

For each selected anchor, ValleScope2 extracts its sequence from the GenMap
input FASTA.

The anchor group ID is based on a strand-invariant canonical sequence:

```text
anchor_seq = extracted sequence
revcomp_seq = reverse_complement(anchor_seq)
canonical_seq = min(anchor_seq, revcomp_seq)
anchor_group_id = SHA-256(canonical_seq)
```

Anchors containing `N` are retained in `grouped_anchors.tsv` but are not
assigned to a group.

Outputs:

```text
grouped_anchors.tsv
anchor_groups.tsv
```

### 5. Structural tokens and contexts

Groupable anchors are converted to an alternating token stream:

```text
A:<anchor_group_id>
D:<distance_bin>
A:<anchor_group_id>
...
```

Distance tokens use floor-based bins controlled by
`--distance-bin-size`.

For each anchor, context extends `--context-radius-tokens` tokens to the left
and right. With the default radius 50, each context has at most 101 tokens.

Sample-specific token-level minimal unique substrings are counted when fully
contained within the context. Contexts and tMUS substrings are canonicalized
against reversed token order.

Outputs:

```text
structural_tokens.tsv
structural_tokens.meta.json
anchor_contexts.tsv
context_groups.tsv
tmus.tsv
structural_contexts.meta.json
```

`anchor_contexts.tsv` includes `forward_context_tokens` for strand-invariant
edit-distance comparison.

### 6. Context-based correspondence assignment

For each ordered target-query sample pair, ValleScope2 assigns query anchors to
target anchors using:

- exact canonical context key match
- fallback search among anchors with the same anchor group ID
- bounded token-level edit distance for candidate search
- `candidate_score = q_alpha_prime + r_alpha_prime - edit_distance`

Candidates are retained when:

```text
candidate_score > --min-candidate-score
```

Directed assignments are merged per sample pair into final correspondences.

`--pair-merge-mode union` keeps anchor pairs detected in either direction.
`--pair-merge-mode reciprocal` keeps only reciprocal pairs.

Outputs:

```text
assignments.tsv
assignments.meta.json
anchor_correspondences.tsv
anchor_correspondences.meta.json
```

### 7. Bundle chaining

ValleScope2 chains merged anchor correspondences by dynamic programming.

For candidate `i`, up to `--chain-predecessors` previous candidates in target
coordinate order are considered. A predecessor is valid only if target and
query order are both consistent with the assignment strand.
In addition, both the target-side gap and query-side gap between adjacent
chained anchors must be no larger than `--max-chain-gap`.
The ratio `max(dr, dq) / min(dr, dq)` must also be no larger than
`--chain-max-gap-ratio`.

The chaining recurrence is:

```text
dp[i] = candidate_score[i] + max_j { dp[j] - gap_cost(j, i) }
```

where:

```text
x = abs(dr - dq) / gap_unit
gap_cost = gap_weight * x + log2(1 + x)
gap_unit = 10
```

With `--gap-cost-model relative`, chaining instead uses:

```text
x = abs(dr - dq) / max(1, min(dr, dq))
gap_cost = gap_weight * x + sqrt(1 + x) - 1
```

Multiple chains are extracted iteratively after removing used candidates.
`--min-chain-both-anchors` can require at least N reciprocally supported
(`support_direction=both`) anchor pairs per emitted chain; this is useful for
filtering one-direction-only repeat-copy chains during QC.

Outputs:

```text
chains.tsv
chain_anchors.tsv
chains.meta.json
```

### 8. Refinement

After first-pass chaining, unused correspondence candidates plus `primary` and
`ambiguous` rows from `assignments.tsv` are re-chained:

- around each first-pass chain within `--refinement-window`
- between adjacent first-pass chains from the same sample pair, sequence pair,
  and strand

Outputs:

```text
refined_chains.tsv
refined_chain_anchors.tsv
refined_chains.meta.json
```

### 9. Bundle-level base alignment

By default, ValleScope2 uses chain anchor pairs as fixed points for
anchor-guided base alignment. Inter-anchor segments are aligned with WFA2-lib,
while anchor sequences are first checked against the FASTA sequence and then
emitted as exact-match CIGAR operations. The final PAF is written to standard
output, so users can choose the filename with shell redirection, for example
`vallescope2 reference.fa query.fa > aln.paf`.
Both first-pass chains (`chains.tsv`) and refined chains (`refined_chains.tsv`)
are used as bundle intervals.

For `-` strand chains, the query interval is reverse-complemented before
alignment.

Before WFA2 alignment, adjacent bundles on the same sample pair, sequence pair,
and strand are tested for gap patching. If both the ref and query gaps are at
most `--max-patch-gap-bp`, and extension windows keep at least
`--min-patch-identity` by Levenshtein edit distance, the two bundles are merged
and aligned as one patched bundle.

After gap patching, the final bundle set is overlap-trimmed immediately before
base alignment. Lower-scoring bundles are trimmed only when both the ref and
query intervals overlap a higher-scoring bundle by at least
`--chain-trim-overlap`. One-sided overlaps are retained because they may
represent duplication or other structural variation.

Large intervals are skipped if either side exceeds:

```text
--max-bundle-align-bp
```

WFA2-lib is also given a per-alignment memory limit:

```text
--max-wfa-memory-gb
```

Final output:

```text
stdout PAF
```

The final PAF includes:

```text
cg:Z:<CIGAR>
ch:i:<chain_id>
bt:Z:<primary|refined|patched>
pg:i:<number_of_patched_gaps>
tr:i:<number_of_bundle_trim_events>
am:Z:anchor_guided
os:Z:<original_assigned_strand>
es:Z:<effective_sequence_inferred_strand>
af:i:<forward_anchor_match_count>
ar:i:<reverse_anchor_match_count>
ai:i:<incompatible_anchor_match_count>
as:i:<alignment_score>
```

Bundles without any validated fixed anchor pair are skipped rather than being
aligned with whole-bundle fallback.

In debug mode, strand inference and skipped-bundle reports are also written:

```text
strand_conflicts.tsv
skipped_bundles.tsv
```

## Output summary

Default persistent outputs:

```text
none
```

Standard output:

```text
PAF with cg:Z CIGAR
```

Optional persistent outputs:

```text
window_scores.tsv
```

Debug workspace:

```text
vallescope2-debug/
```

## Source layout

```text
include/vallescope2/
  alignment/       bundle-level WFA2 alignment interface
  anchor/          GenMap runner, window scoring, anchor selection, outputs
  bundle/          DP chaining and refinement
  common/          process, SHA-256, workspace helpers
  context/         anchor grouping, structural tokens, structural contexts
  correspondence/  context-based assignment and pair merging
  input/           FASTA/index/input preparation
  interface/       command-line options

src/
  alignment/
  anchor/
  bundle/
  common/
  context/
  correspondence/
  input/
  interface/

external/WFA2-lib/  WFA2-lib git submodule
scripts/            temporary/debug visualization helpers
test/               lightweight regression tests
```

## Development notes

Temporary visualization scripts under `scripts/` can convert chain TSVs to PAF
and draw simple PNG ribbon plots. They are for inspection during development
and are not part of the final alignment pipeline.

Run the lightweight regression test with:

```bash
build/anchor_grouping_test
```
