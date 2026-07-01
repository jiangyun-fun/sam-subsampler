//! BAM/CRAM/SAM input-output primitives for the two-pass subsampler.
//!
//! - [`detect_format`] maps an output path extension to an htslib format.
//! - [`read_unique_qnames_by_ref`] is pass 1: stream the file once and collect
//!   the *unique* qname set per reference (dedup happens on insert, so a read
//!   with several records — mate, supplementary — is one selection unit).
//! - [`tag_and_write`] is pass 2: re-read the file and write every record out,
//!   adding a BAM aux tag to records whose qname was selected.

use crate::error::{AppError, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use log::{debug, info, trace};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Aux;
use std::collections::{HashMap, HashSet};
use std::path::Path;

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
/// Unmapped records and records with no reference (`tid < 0`) are skipped for
/// selection but still counted in `total`.
pub fn read_unique_qnames_by_ref(path: &Path) -> Result<(HashMap<String, HashSet<Vec<u8>>>, u64)> {
    let mut reader = bam::Reader::from_path(path)?;
    let header = reader.header().to_owned();

    let mut by_ref: HashMap<String, HashSet<Vec<u8>>> = HashMap::new();
    let mut total: u64 = 0;
    for result in reader.records() {
        let record = result?;
        total += 1;
        if !record.is_unmapped() && record.tid() >= 0 {
            let name = String::from_utf8(header.tid2name(record.tid() as u32).to_vec())?;
            by_ref
                .entry(name)
                .or_default()
                .insert(record.qname().to_vec());
        }
    }
    Ok((by_ref, total))
}

/// Pass 2: re-read `input` and write every record to `output`, tagging records
/// whose qname is in `selected` with `Aux::I32(1)` under `tag`.
///
/// `total_records` drives the progress bar (shown only when `show_progress`).
pub fn tag_and_write(
    input: &Path,
    output: &Path,
    output_format: bam::Format,
    reference: Option<&Path>,
    selected: &HashSet<Vec<u8>>,
    tag: &[u8],
    total_records: u64,
    show_progress: bool,
) -> Result<()> {
    let mut reader = bam::Reader::from_path(input)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = build_writer(output, &header, output_format, reference)?;

    let pb = ProgressBar::new(total_records);
    if show_progress {
        pb.set_style(progress_style());
    } else {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }

    let mut written: u64 = 0;
    for result in reader.records() {
        let mut record = result?;
        if selected.contains(record.qname()) {
            trace!(
                "tagging {} with {}",
                String::from_utf8_lossy(record.qname()),
                String::from_utf8_lossy(tag)
            );
            // Aux is not Copy; construct a fresh value per record.
            record.push_aux(tag, Aux::I32(1))?;
        }
        writer.write(&record)?;
        pb.inc(1);
        written += 1;
    }
    pb.finish_and_clear();

    debug!("wrote {written} records to {output:?}");
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
}
