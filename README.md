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
build/vallescope2 -db 50 --debug genome.fa
build/vallescope2 -cr 25 --debug genome.fa
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
`-db/--distance-bin-size` and defaults to 50 bp.

Each anchor context extends `-cr/--context-radius-tokens` tokens to either side
(default 25, for at most 51 tokens). Sample-specific token-level minimal unique
substrings are counted only when fully contained in the context. Substrings
and whole contexts are canonicalized against their reversed token order.

Window scores are not written by `--debug` because the file can be larger than
the genome. Use `--dump-window-scores` only when the full TSV is needed.

One FASTA selects self-comparison mode. Two FASTAs select target-query mode in
the given order. `--combine` combines two or more FASTAs for self-comparison;
without it, more than two inputs are rejected. Every sequence ID is normalized
to `sample#1#original_id` in combine mode. One or two inputs are used directly
without copying; HTSlib creates temporary indexes for them. With `--debug`,
indexes, metadata, and any combined FASTA are retained in `vallescope2-debug/`.
