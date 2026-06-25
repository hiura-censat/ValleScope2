# ValleScope2

ValleScope2 is under development. The current implementation validates one or
more FASTA files and prepares a normalized multi-FASTA for the later pipeline.

## Build

```bash
cmake -S . -B build -DHTSLIB_ROOT=/path/to/htslib
cmake --build build
```

Using the local wfmash environment created during development:

```bash
cmake -S . -B build -DHTSLIB_ROOT="$PWD/.tools/wfmash"
cmake --build build
```

## Run

```bash
build/vallescope2 genome.fa
build/vallescope2 --debug sample1.fa sample2.fa
build/vallescope2 --combine --debug sample1.fa sample2.fa
build/vallescope2 --combine --debug sample1.fa sample2.fa sample3.fa
```

ValleScope2 runs GenMap after preparing the input. The defaults are a k-mer
length of 6 and zero errors:

```bash
build/vallescope2 -K 6 -E 0 --debug genome.fa
build/vallescope2 --kmer-length 8 --kmer-errors 1 --debug genome.fa
build/vallescope2 -K 6 -E 0 -al 50 --debug genome.fa
build/vallescope2 -md 100 -td 6000 --debug genome.fa
build/vallescope2 -t 8 --debug genome.fa
build/vallescope2 -db 10 --debug genome.fa
build/vallescope2 -cr 50 --debug genome.fa
build/vallescope2 -mcs -10 -pm 5 --debug genome.fa
build/vallescope2 --pair-merge-mode union --debug genome.fa
build/vallescope2 --chain-predecessors 50 --gap-weight 1 --debug genome.fa
build/vallescope2 --min-chain-anchors 10 --min-chain-score 0 --debug genome.fa
build/vallescope2 --refinement-window 50000 --refinement-min-chain-anchors 5 --debug genome.fa
build/vallescope2 --base-align --max-bundle-align-bp 50000 --debug genome.fa
```

`-t/--threads` controls the number of threads passed to GenMap (`-T`); the
default is 20.

The executable is found on `PATH`, in `.tools/genmap/bin/genmap`, or can be
specified explicitly with `--genmap PATH`.

Window scoring uses all k-mers fully contained in each anchor window. For an
anchor length `W` and k-mer length `K`, each score is the mean of
`log(1 + frequency)` over `W - K + 1` values. The anchor/window length defaults
to 50 and can be changed with `-al/--anchor-length`.

Anchor candidates are processed from low to high score. Candidates within the
minimum center distance of an accepted anchor are removed. The defaults are
100 bp (`-md/--min-center-distance`) and 6000 anchors/Mb
(`-td/--target-density`). The persistent outputs are `anchors.bed`,
`anchors.meta.json`, `grouped_anchors.tsv`, `anchor_groups.tsv`, and
`structural_tokens.tsv`, `structural_tokens.meta.json`, and `genmap.log`.
The structural-context phase additionally writes `anchor_contexts.tsv`,
`context_groups.tsv`, `tmus.tsv`, and `structural_contexts.meta.json`.
The context-based correspondence phase writes `assignments.tsv` and
`assignments.meta.json`. Candidate anchors are retained when
`candidate_score > min_candidate_score`, where
`candidate_score = q_alpha_prime + r_alpha_prime - edit_distance`.
`-mcs/--min-candidate-score` defaults to -10. `-pm/--primary-margin` controls
the score margin required to call a unique primary assignment and defaults to
5. Primary directed assignments are merged per sample pair into
`anchor_correspondences.tsv`; `--pair-merge-mode` can be `union` or
`reciprocal` and defaults to `union`.
The bundle-chaining phase uses `anchor_correspondences.tsv` and writes
`chains.tsv`, `chain_anchors.tsv`, and `chains.meta.json`. Each sample pair and
assignment strand emits one or more DP chains. The predecessor search limit
defaults to 50 and is controlled by `--chain-predecessors`; the gap penalty
weight defaults to 1 and is controlled by `--gap-weight`. Chains are
iteratively extracted after removing used candidates; `--min-chain-anchors`
and `--min-chain-score` control when extraction stops. The minimum first-pass
chain length defaults to 10 anchors. Unused correspondence candidates plus
primary and ambiguous rows from `assignments.tsv` are then re-chained around
first-pass chains and between adjacent first-pass chains. This refinement
writes `refined_chains.tsv`, `refined_chain_anchors.tsv`, and
`refined_chains.meta.json`; `--refinement-window` controls the around-chain
extension and defaults to 50000 bp. `--refinement-min-chain-anchors` controls
the minimum refined chain length and defaults to 5 anchors.

With `--base-align`, ValleScope2 extracts each first-pass chain interval from
the GenMap input FASTA, reverse-complements the query interval for `-` strand
chains, globally aligns the whole bundle interval, and writes
`bundle_alignments.paf` with a `cg:Z:` CIGAR tag plus
`bundle_alignments.meta.json`. The current build uses an internal global
alignment fallback unless WFA2 headers are available; `--max-bundle-align-bp`
skips very large bundle intervals and defaults to 50000 bp.

For each selected anchor, ValleScope2 extracts its sequence from the normalized
GenMap input FASTA and defines a strand-invariant canonical sequence as the
lexicographically smaller of the forward sequence and its reverse complement.
The SHA-256 of the exact canonical sequence is used as the stable anchor group
ID. Anchors containing `N` remain in `grouped_anchors.tsv` but are not assigned
to a group.

Groupable anchors are converted into an alternating structural token stream.
Anchor tokens are written as `A:<anchor_group_id>`. The center-to-center
distance between adjacent groupable anchors is divided into floor-based bins
and written as `D:<bin>`. The bin size is controlled by
`-db/--distance-bin-size` and defaults to 10 bp.

Each anchor context extends `-cr/--context-radius-tokens` tokens to either side
(default 50, for at most 101 tokens). Sample-specific token-level minimal unique
substrings are counted only when fully contained in the context. Substrings
and whole contexts are canonicalized against their reversed token order.
`anchor_contexts.tsv` also includes `forward_context_tokens`, which preserves
the observed token order for strand-invariant edit-distance assignment.

Window scores are not written by `--debug` because the file can be larger than
the genome. Use `--dump-window-scores` only when the full TSV is needed.

One FASTA selects self-comparison mode. Two FASTAs select target-query mode in
the given order. `--combine` combines two or more FASTAs for self-comparison;
without it, more than two inputs are rejected. Every sequence ID is normalized
to `sample#1#original_id` in combine mode. One or two inputs are used directly
without copying; HTSlib creates temporary indexes for them. With `--debug`,
indexes, metadata, and any combined FASTA are retained in `vallescope2-debug/`.
