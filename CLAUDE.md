# sam-subsampler — project notes for Claude Code

Rust BAM/CRAM/SAM subsampler built on `rust-htslib`. The pixi workspace pins the
toolchain (htslib, samtools, clang 18, gcc forced to C17) so builds are
reproducible. GitHub: `jiangyun-fun/sam-subsampler`.

## Commands (always via `pixi run`)

- `pixi run cargo test --all` — unit + integration suite
- `pixi run cargo clippy --all-targets -- -D warnings`
- `pixi run cargo fmt --all -- --check`
- `pixi run cargo build --release`

`clippy::unwrap_used` is workspace-warn (`Cargo.toml`) and relaxed only under
`#[cfg(test)]`; never add `unwrap`/`expect` to production paths.

## Release process

Releases are **tag-driven** via `.github/workflows/release.yml`.

1. Bump the version in `Cargo.toml` (and `Cargo.lock`) as its own
   `chore(release): bump version to X.Y.Z` commit — mirror the existing history
   (`feat: …` first, then the version-bump commit).
2. Push `main`; wait for **CI** (`ci.yml`: fmt + clippy + test) to go green.
3. `git tag -a vX.Y.Z -m "Release X.Y.Z"` on the release commit, then
   `git push origin vX.Y.Z`.
4. `release.yml` builds the Linux x86_64 binary and publishes the GitHub Release
   with the `sam-subsampler-vX.Y.Z-x86_64-linux.tar.gz` tarball. No manual
   `gh release create` is needed.

### Registry publishing — NOT automated by the workflow

- **crates.io**: `cargo publish` is runnable and is **run by the user** (not by
  CI). A new GitHub Release/tag does **not** imply a crates.io publish — flag it
  as a separate manual step.
- **bioconda**: version bumps are handled **automatically by the bioconda bot**.
  Do **not** open a manual `bioconda/bioconda-recipes` PR.

## Conventions

- Conventional Commits; author `Yun Jiang <jiangyunacw@gmail.com>`. No
  `Co-Authored-By` trailer.
- Tag style: **annotated** (`git tag -a`), placed on the `chore(release)` commit.
- Selection logic lives in `src/selection.rs` (pure, reference-agnostic). Pass 1
  (`bam_io::read_unique_qnames_by_ref`) collects unique qnames per reference,
  pooling unmapped reads under a `*` bucket. Pass 2 (`bam_io::tag_and_write`)
  writes per `OutputMode` (default `KeepSelected`; `--keep-all` ⇒ `TagInPlace`).
- Unmapped reads are subsampled as a first-class `*` bucket in every mode.
