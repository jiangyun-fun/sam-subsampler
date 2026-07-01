//! sam-subsampler: two-pass BAM/CRAM subsampler that tags selected reads in place.
//!
//! Idiomatic Rust rewrite of `bam_subsampler`. See `README.md` for usage and
//! design notes (per-reference reservoir sampling, qname-dedup bias fix, CRAM
//! reference handling).

pub mod bam_io;
pub mod cli;
pub mod config;
pub mod error;
pub mod selection;

pub use error::{AppError, Result};

use clap::Parser;
use config::SubsamplePlan;
use log::{error, info};
use std::io::IsTerminal;

/// Entry point invoked by `main`. Returns the process exit code.
pub fn run() -> i32 {
    let cli = cli::Cli::parse();
    cli::setup_logger(cli.verbose);
    match try_run(&cli) {
        Ok(()) => 0,
        Err(e) => {
            error!("{e}");
            1
        }
    }
}

fn try_run(cli: &cli::Cli) -> Result<()> {
    cli.validate()?;
    let plan = build_plan(cli)?;

    info!(
        "pass 1: reading {:?} to collect unique qnames per reference",
        cli.input_bam
    );
    let (qnames_by_ref, total) = bam_io::read_unique_qnames_by_ref(&cli.input_bam)?;
    let refs = qnames_by_ref.len();
    info!("pass 1: {total} records across {refs} references");

    info!("selecting qnames (seed {})", cli.seed);
    let selected = selection::select_per_reference(qnames_by_ref, &plan, cli.seed);
    info!("selected {} unique qnames for tagging", selected.len());

    let format = bam_io::detect_format(&cli.output_bam)?;
    info!(
        "pass 2: writing {:?} (format {format:?}) with tag '{}'",
        cli.output_bam, cli.add_ssub
    );

    let show_progress = cli.verbose >= 1 && std::io::stderr().is_terminal();
    bam_io::tag_and_write(
        &cli.input_bam,
        &cli.output_bam,
        format,
        cli.reference.as_deref(),
        &selected,
        cli.add_ssub.as_bytes(),
        total,
        show_progress,
    )?;

    info!(
        "done: tagged {} unique qnames across {refs} references",
        selected.len()
    );
    Ok(())
}

/// Build the [`SubsamplePlan`] from the CLI. `--count` and `--config` conflict
/// (enforced by clap), so the `(Some, Some)` case is unreachable.
fn build_plan(cli: &cli::Cli) -> Result<SubsamplePlan> {
    match (cli.count, cli.config.as_ref()) {
        (Some(n), _) => Ok(SubsamplePlan::Global(n)),
        (None, Some(path)) => Ok(SubsamplePlan::PerRef(config::load_config_csv(path)?)),
        (None, None) => Ok(SubsamplePlan::Default),
    }
}
