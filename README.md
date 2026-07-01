# sam-subsampler

[![CI](https://github.com/jiangyun-fun/sam-subsampler/actions/workflows/ci.yml/badge.svg)](https://github.com/jiangyun-fun/sam-subsampler/actions/workflows/ci.yml)
[![release](https://github.com/jiangyun-fun/sam-subsampler/actions/workflows/release.yml/badge.svg)](https://github.com/jiangyun-fun/sam-subsampler/actions/workflows/release.yml)
[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![release](https://img.shields.io/github/v/release/jiangyun-fun/sam-subsampler)](https://github.com/jiangyun-fun/sam-subsampler/releases)

Subsample reads from a BAM/CRAM/SAM file by **per-reference count** and **tag
the selected reads in place** — the output is the full file with a BAM aux tag
added to a randomly chosen subset. It does **not** filter.

Two things set it apart from `samtools view -s`, `rasusa`, `picard DownSampleSam`,
and friends:

1. **Per-reference counts** — pick a different number of reads from each
   chromosome (`--config` CSV), or one global count (`--count`).
2. **Tagging, not filtering** — every record is written; selected reads are
   marked with a custom 2-character aux tag (e.g. `YS:i:1`) so downstream tools
   see the whole alignment plus a labelled subset.

## Bias fix vs. the original `bam_subsampler`

The predecessor collected read names *per record*, so a paired or
multi-alignment read (one qname on several records) had roughly N× the
selection probability of a single-record read. `sam-subsampler` samples
**unique qnames per reference**, so each read is one selection unit — and when
it is selected, **all of its records** (mate, supplementary) are tagged.

## Install

`rust-htslib` needs htslib and libclang at build time. The repo ships a
[pixi](https://pixi.sh) environment that pins compatible versions (htslib,
samtools, clang 18, gcc forced to C17).

```sh
pixi install            # sets up the toolchain
pixi run cargo build --release
# binary: target/release/sam-subsampler
```

If you already have htslib (≥1.10) and libclang on your system:

```sh
cargo install sam-subsampler        # from crates.io (no clone needed)
# or, from a clone of this repo:
cargo install --locked --path .     # or: cargo build --release
```

A prebuilt **Bioconda** package is in review — once
[bioconda/bioconda-recipes#66826](https://github.com/bioconda/bioconda-recipes/pull/66826)
merges:

```sh
conda install -c conda-forge -c bioconda sam-subsampler
```

## Usage

Subsample 1000 reads per reference into `out.bam`, tagging selected reads with
`YS`:

```sh
sam-subsampler -i in.bam -o out.bam --count 1000 --add-ssub YS --seed 42
```

Different counts per reference via a CSV config:

```sh
sam-subsampler -i in.bam -o out.bam --config refs.csv --add-ssub YS
```

```text
# refs.csv
seq_name,subsample_count
chr1,5000
chr2,2500
chrX,
```

A blank `subsample_count` (or a reference absent from the CSV) falls back to
the default of 1000.

Write CRAM (requires a reference with a `.fai` index):

```sh
samtools faidx ref.fa
sam-subsampler -i in.bam -o out.cram --reference ref.fa --count 1000 --add-ssub YS
```

Stream BAM to stdout (`-`, only for `.bam` output):

```sh
sam-subsampler -i in.bam -o - --count 100 --add-ssub YS | samtools view -b > tagged.bam
```

Add `-v` for info logging and a progress bar; repeat for more detail.

## CLI reference

| Flag            | Type     | Required | Default | Notes |
|-----------------|----------|----------|---------|-------|
| `-i, --input-bam`   | path | yes | — | BAM/CRAM/SAM file; stdin (`-`) not supported (the file is read twice) |
| `-o, --output-bam`  | path | yes | — | `-` ⇒ stdout (BAM); extension picks format (`.bam`/`.cram`/`.sam`) |
| `--config`     | path | no | — | Per-reference CSV; conflicts with `--count` |
| `--count`      | u32  | no | — | Global per-reference count |
| `--add-ssub`   | str  | yes | — | 2-char aux tag (letter then letter/digit, e.g. `YS`) |
| `--reference`  | path | no | — | Required for `.cram` output; `.fai` must exist beside it |
| `--seed`       | u64  | no | 42 | RNG seed |
| `-v`           | count| no | 0 | Verbosity (`-v`, `-vv`, `-vvv`) |

If neither `--count` nor `--config` is given, every reference uses the default
of 1000.

## Reproducibility

The selection is a pure function of **(input file, plan, seed)**: references
and qnames are sorted before sampling, and a single RNG seeded once drives all
references. The same inputs always yield the identical selected set. Note that
**adding or removing a chromosome shifts downstream selections**, because the
single RNG is drawn from sequentially across references.

## Algorithm

1. **Pass 1** — stream the file once, collecting the *unique* qname set per
   reference (unmapped reads are skipped). Dedup happens on insert, so memory
   scales with the number of unique read names, not records.
2. **Select** — for each reference, Vitter's reservoir sampling (Algorithm R)
   draws the target number of qnames from its sorted unique set, using the
   shared seeded RNG.
3. **Pass 2** — re-read the file and write every record out; records whose
   qname was selected get `Aux::I32(1)` under the chosen tag.

## Testing

```sh
pixi run cargo test            # 54 tests: unit + integration
pixi run cargo clippy --all-targets -- -D warnings
pixi run cargo fmt --all -- --check
```

The integration test builds a BAM from a SAM string with rust-htslib itself
(no `samtools` needed) and checks: record count is preserved (tagging, not
filtering); the right number of unique qnames is tagged per reference;
unmapped reads are never tagged; a paired read is one selection unit with both
mates tagged; unselected reads carry no tag; the tag value is `i32(1)`; and
the same seed reproduces the set.

## License

MIT.
