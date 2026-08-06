//! BAM/CRAM/SAM input-output primitives for the two-pass subsampler.
//!
//! - [`detect_format`] maps an output path extension to an htslib format.
//! - [`read_unique_qnames_by_ref`] is pass 1: stream the file once and collect
//!   the *unique* qname set per reference (dedup happens on insert, so a read
//!   with several records — mate, supplementary — is one selection unit).
//!   Unmapped reads are pooled under [`UNMAPPED_BUCKET`] (`*`).
//! - [`tag_and_write`] is pass 2: re-read the file and write records out under
//!   the chosen [`OutputMode`], adding a BAM aux tag to selected records.

use crate::error::{AppError, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use log::{info, trace, warn};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Aux;
use std::collections::{HashMap, HashSet};
use std::path::Path;

/// Sentinel bucket key for unmapped reads.
///
/// Unmapped records carry no reference, so they are pooled under SAM's reserved
/// `*` RNAME and become a first-class selection unit in every mode (`--count`,
/// `--config`, `--total-count`, `--ratio`). A real reference can never be named
/// `*` (the SAM spec reserves it), and `tid2name` is only ever called for
/// `tid >= 0` (a real target), so this sentinel cannot collide with a
/// reference-name bucket.
pub const UNMAPPED_BUCKET: &str = "*";

/// How pass 2 emits records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    /// Write only records whose qname was selected, tagging each — the default.
    /// The output is a true subsample (smaller than the input).
    KeepSelected,
    /// Write every record, tagging the selected subset (`--keep-all`; the
    /// pre-0.3 tag-in-place behavior).
    TagInPlace,
}

/// Infer the output format from a path extension.
///
/// No extension (e.g. stdout `-`) ⇒ BAM.
pub fn detect_format(path: &Path) -> Result<bam::Format> {
    match path.extension().and_then(|e| e.to_str()) {
        None | Some("bam") => Ok(bam::Format::Bam),
        Some("cram") => Ok(bam::Format::Cram),
        Some("sam") => Ok(bam::Format::Sam),
        Some(other) => Err(AppError::Argument(format!(
            "unsupported output extension '.{other}'; use .bam, .cram, or .sam"
        ))),
    }
}

/// Pass 1: stream `path` once, returning the unique-qname set per reference
/// and the total number of records seen.
///
/// Mapped reads are keyed by reference name; unmapped reads are pooled under
/// the [`UNMAPPED_BUCKET`] sentinel (`*`) so they are subsampled like a
/// reference. Records that are neither (`tid < 0` without the unmapped flag —
/// invalid SAM) are skipped for selection but still counted in `total`.
pub fn read_unique_qnames_by_ref(path: &Path) -> Result<(crate::QnamesByRef, u64)> {
    let mut reader = bam::Reader::from_path(path)?;
    let header = reader.header().to_owned();

    let mut by_ref: HashMap<String, HashSet<Vec<u8>>> = HashMap::new();
    let mut total: u64 = 0;
    let mut skipped: u64 = 0;
    for result in reader.records() {
        let record = result?;
        total += 1;
        let bucket = if record.is_unmapped() {
            // Unmapped reads have no reference → pool under the reserved `*`.
            UNMAPPED_BUCKET.to_string()
        } else if record.tid() >= 0 {
            String::from_utf8(header.tid2name(record.tid() as u32).to_vec())?
        } else {
            // Mapped flag set but tid < 0: invalid SAM. Skip defensively (count
            // and warn once below) rather than silently inventing a bucket.
            skipped += 1;
            continue;
        };
        by_ref
            .entry(bucket)
            .or_default()
            .insert(record.qname().to_vec());
    }
    if skipped > 0 {
        warn!(
            "skipped {skipped} record(s) with invalid mapping state \
             (mapped flag set but reference id < 0); not eligible for selection"
        );
    }
    Ok((by_ref, total))
}

/// Parameters for [`tag_and_write`].
pub struct TagWrite<'a> {
    pub input: &'a Path,
    pub output: &'a Path,
    pub output_format: bam::Format,
    pub reference: Option<&'a Path>,
    pub selected: &'a HashSet<Vec<u8>>,
    pub tag: &'a [u8],
    pub total_records: u64,
    pub mode: OutputMode,
    pub show_progress: bool,
}

/// Pass 2: re-read `input` and write records to `output`, tagging those whose
/// qname is in `selected` with `Aux::I32(1)` under `tag`.
///
/// Under [`OutputMode::KeepSelected`] (the default) only selected records are
/// written — a true subsample. Under [`OutputMode::TagInPlace`] (`--keep-all`)
/// every record is written, with the selected subset tagged.
///
/// `total_records` drives the progress bar (shown only when `show_progress`).
pub fn tag_and_write(args: TagWrite<'_>) -> Result<()> {
    let mut reader = bam::Reader::from_path(args.input)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = build_writer(args.output, &header, args.output_format, args.reference)?;

    let pb = ProgressBar::new(args.total_records);
    if args.show_progress {
        pb.set_style(progress_style());
    } else {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }

    let mut written: u64 = 0;
    for result in reader.records() {
        let mut record = result?;
        let is_selected = args.selected.contains(record.qname());
        // Both modes tag selected records identically; only the write decision
        // differs. KeepSelected drops non-selected records (true subsample);
        // TagInPlace (--keep-all) writes every record.
        if is_selected {
            trace!(
                "tagging {} with {}",
                String::from_utf8_lossy(record.qname()),
                String::from_utf8_lossy(args.tag)
            );
            // Aux is not Copy; construct a fresh value per record.
            record.push_aux(args.tag, Aux::I32(1))?;
        }
        let write_record = match args.mode {
            OutputMode::KeepSelected => is_selected,
            OutputMode::TagInPlace => true,
        };
        if write_record {
            writer.write(&record)?;
            written += 1;
        }
        pb.inc(1);
    }
    pb.finish_and_clear();

    info!("wrote {written} records to {:?}", args.output);
    Ok(())
}

/// Construct the output writer, attaching a reference for CRAM when provided.
fn build_writer(
    output: &Path,
    header: &bam::Header,
    output_format: bam::Format,
    reference: Option<&Path>,
) -> Result<bam::Writer> {
    let mut writer = if output.as_os_str() == "-" {
        info!("writing BAM to stdout");
        bam::Writer::from_stdout(header, output_format)?
    } else {
        info!("writing {output:?} (format {output_format:?})");
        bam::Writer::from_path(output, header, output_format)?
    };
    if let Some(ref_path) = reference {
        writer.set_reference(ref_path)?;
    }
    Ok(writer)
}

fn progress_style() -> ProgressStyle {
    ProgressStyle::with_template("{elapsed} {wide_bar} {pos}/{len} records ({percent}%)")
        .unwrap_or_else(|_| ProgressStyle::default_bar())
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn detect_format_by_extension() {
        assert_eq!(detect_format(Path::new("a.bam")).unwrap(), bam::Format::Bam);
        assert_eq!(
            detect_format(Path::new("a.cram")).unwrap(),
            bam::Format::Cram
        );
        assert_eq!(detect_format(Path::new("a.sam")).unwrap(), bam::Format::Sam);
        // no extension (stdout) -> BAM
        assert_eq!(detect_format(Path::new("-")).unwrap(), bam::Format::Bam);
        assert!(detect_format(Path::new("a.txt")).is_err());
    }

    #[test]
    fn unmapped_bucket_sentinel_is_star() {
        // Documents the contract relied on by read_unique_qnames_by_ref: the
        // unmapped bucket key is SAM's reserved '*'.
        assert_eq!(UNMAPPED_BUCKET, "*");
    }
}
