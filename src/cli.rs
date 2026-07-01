//! Command-line interface: clap definition, boundary validation, logger setup.

use crate::error::{AppError, Result};
use clap::Parser;
use std::path::{Path, PathBuf};

/// Default RNG seed (genuinely optional reproducibility knob).
pub const DEFAULT_SEED: u64 = 42;

/// Parsed command-line arguments.
#[derive(Parser, Debug)]
#[command(
    name = "sam-subsampler",
    version,
    about = "Subsample reads from a BAM/CRAM by per-reference count and tag selected reads in place"
)]
pub struct Cli {
    /// Input alignment file (BAM/CRAM/SAM). Cannot be stdin: the file is read twice.
    #[arg(short = 'i', long)]
    pub input_bam: PathBuf,

    /// Output file; use '-' for stdout (BAM). Extension selects the format:
    /// `.bam`, `.cram`, or `.sam`.
    #[arg(short = 'o', long)]
    pub output_bam: PathBuf,

    /// Per-reference config CSV (`seq_name,subsample_count`). Conflicts with `--count`.
    #[arg(long, conflicts_with = "count")]
    pub config: Option<PathBuf>,

    /// Global per-reference subsample count applied to every reference.
    #[arg(long)]
    pub count: Option<u32>,

    /// 2-character BAM aux tag added to selected reads (e.g. `YS`).
    #[arg(long)]
    pub add_ssub: String,

    /// Reference FASTA (with a sibling `.fai` index). Required for `.cram` output.
    #[arg(long, value_name = "FASTA")]
    pub reference: Option<PathBuf>,

    /// RNG seed for reproducible subsampling.
    #[arg(long, default_value_t = DEFAULT_SEED)]
    pub seed: u64,

    /// Increase logging verbosity (`-v`, `-vv`, `-vvv`).
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbose: u8,
}

impl Cli {
    /// Validate boundary conditions clap cannot enforce.
    ///
    /// Run this once in `main` before doing any work.
    pub fn validate(&self) -> Result<()> {
        validate_tag(&self.add_ssub)?;
        if self.input_bam.as_os_str() == "-" {
            return Err(AppError::Argument(
                "input cannot be stdin ('-'): the file must be read twice for subsampling".into(),
            ));
        }
        if self.output_is_cram() {
            let reference = self.reference.as_deref().ok_or_else(|| {
                AppError::Argument(
                    "--reference <FASTA> is required for .cram output (CRAM encodes against a reference)"
                        .into(),
                )
            })?;
            let fai = fai_path(reference);
            if !fai.exists() {
                return Err(AppError::Argument(format!(
                    "reference index '{}' not found; create it with `samtools faidx {}`",
                    fai.display(),
                    reference.display()
                )));
            }
        }
        Ok(())
    }

    /// True when the output path's extension indicates CRAM.
    pub fn output_is_cram(&self) -> bool {
        self.output_bam.extension().and_then(|e| e.to_str()) == Some("cram")
    }
}

/// Validate a BAM aux tag name: exactly two ASCII chars, first alphabetic,
/// second alphanumeric (SAM/BAM spec).
fn validate_tag(tag: &str) -> Result<()> {
    let bytes = tag.as_bytes();
    let valid =
        bytes.len() == 2 && bytes[0].is_ascii_alphabetic() && bytes[1].is_ascii_alphanumeric();
    if !valid {
        return Err(AppError::Argument(format!(
            "tag '{tag}' is invalid: provide exactly two ASCII characters (a letter then a letter/digit, e.g. YS)"
        )));
    }
    Ok(())
}

/// Path of the faidx index htslib expects beside a reference (`<fasta>.fai`).
fn fai_path(reference: &Path) -> PathBuf {
    let mut s = reference.as_os_str().to_os_string();
    s.push(".fai");
    PathBuf::from(s)
}

/// Initialize `env_logger` from a verbosity count. Safe to call once.
pub fn setup_logger(verbose: u8) {
    let level = match verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    let _ = env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level))
        .format_timestamp(None)
        .try_init();
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn cli(input: &str, output: &str, add_ssub: &str, reference: Option<&str>) -> Cli {
        Cli {
            input_bam: input.into(),
            output_bam: output.into(),
            config: None,
            count: None,
            add_ssub: add_ssub.into(),
            reference: reference.map(PathBuf::from),
            seed: DEFAULT_SEED,
            verbose: 0,
        }
    }

    // --- tag validation ---

    #[rstest]
    #[case("YS", true)]
    #[case("Z9", true)]
    #[case("ab", true)]
    #[case("Y", false)]
    #[case("YSX", false)]
    #[case("", false)]
    #[case("1S", false)] // first char must be a letter
    #[case("Y-", false)] // second char must be alphanumeric
    #[case("Y S", false)]
    fn tag_validation(#[case] tag: &str, #[case] ok: bool) {
        assert_eq!(validate_tag(tag).is_ok(), ok, "tag {tag:?}");
    }

    // --- input / output validation ---

    #[test]
    fn stdin_input_rejected() {
        let err = cli("-", "out.bam", "YS", None).validate().unwrap_err();
        assert!(err.to_string().contains("stdin"), "got: {err}");
    }

    #[test]
    fn bam_output_needs_no_reference() {
        assert!(cli("in.bam", "out.bam", "YS", None).validate().is_ok());
    }

    #[test]
    fn cram_output_without_reference_rejected() {
        let err = cli("in.bam", "out.cram", "YS", None)
            .validate()
            .unwrap_err();
        assert!(err.to_string().contains("--reference"), "got: {err}");
    }

    #[test]
    fn cram_output_with_indexed_reference_ok() {
        let dir = tempfile::tempdir().unwrap();
        let ref_path = dir.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        std::fs::File::create(fai_path(&ref_path)).unwrap();
        assert!(
            cli("in.bam", "out.cram", "YS", ref_path.to_str())
                .validate()
                .is_ok()
        );
    }

    #[test]
    fn cram_output_with_unindexed_reference_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let ref_path = dir.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let err = cli("in.bam", "out.cram", "YS", ref_path.to_str())
            .validate()
            .unwrap_err();
        assert!(err.to_string().contains("samtools faidx"), "got: {err}");
    }

    #[test]
    fn fai_path_appends_fai_suffix() {
        assert_eq!(
            fai_path(Path::new("/x/ref.fa")),
            PathBuf::from("/x/ref.fa.fai")
        );
        assert_eq!(fai_path(Path::new("/x/ref")), PathBuf::from("/x/ref.fai"));
    }

    #[test]
    fn output_is_cram_detects_extension() {
        assert!(cli("in.bam", "out.cram", "YS", None).output_is_cram());
        assert!(!cli("in.bam", "out.bam", "YS", None).output_is_cram());
        assert!(!cli("in.bam", "-", "YS", None).output_is_cram());
    }
}
