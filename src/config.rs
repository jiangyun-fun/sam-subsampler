//! Per-reference subsample-count configuration: CSV parsing + resolution.
//!
//! Configuration CSV format (header required):
//! ```text
//! seq_name,subsample_count
//! chr1,5000
//! chr2,2500
//! chrX,
//! ```
//! A blank `subsample_count` cell falls back to [`DEFAULT_SUBSAMPLE_COUNT`].
//! References absent from the CSV also fall back to the default.

use crate::error::{AppError, Result};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

/// Fallback count when none is specified: blank CSV cell, unlisted reference,
/// or neither `--count` nor `--config` provided.
pub const DEFAULT_SUBSAMPLE_COUNT: u32 = 1000;

/// A row in the configuration CSV.
///
/// `subsample_count` is read as an `Option<String>` so that a blank cell is
/// unambiguously distinguished from a malformed value (explicit parse + clear
/// error rather than serde's implicit failure).
#[derive(Debug, Deserialize)]
struct ConfigRecord {
    seq_name: String,
    subsample_count: Option<String>,
}

/// How many reads to sample for a given reference.
#[derive(Debug, Clone)]
pub enum SubsamplePlan {
    /// One count applied to every reference (`--count N`).
    Global(u32),
    /// Per-reference counts from `--config`; unlisted refs use the default.
    PerRef(HashMap<String, u32>),
    /// No `--count` and no `--config`: default everywhere.
    Default,
}

impl SubsamplePlan {
    /// Number of reads to sample for `ref_name`.
    pub fn count_for(&self, ref_name: &str) -> usize {
        match self {
            Self::Global(n) => *n as usize,
            Self::PerRef(map) => map
                .get(ref_name)
                .copied()
                .unwrap_or(DEFAULT_SUBSAMPLE_COUNT) as usize,
            Self::Default => DEFAULT_SUBSAMPLE_COUNT as usize,
        }
    }
}

/// Parse a configuration CSV into a per-reference count map.
///
/// On a duplicate `seq_name`, the last row wins. Row numbers in errors are
/// 1-based and include the header line (so the first data row is row 2).
pub fn load_config_csv(path: &Path) -> Result<HashMap<String, u32>> {
    let mut rdr = csv::Reader::from_path(path)?;
    let mut map = HashMap::new();
    for (i, result) in rdr.deserialize().enumerate() {
        let row = i + 2;
        let record: ConfigRecord =
            result.map_err(|e| AppError::Config(format!("row {row}: could not parse: {e}")))?;
        let count = parse_count(&record.subsample_count, row)?;
        map.insert(record.seq_name, count);
    }
    Ok(map)
}

/// Parse an optional count cell: missing/empty ⇒ default, otherwise a u32.
fn parse_count(raw: &Option<String>, row: usize) -> Result<u32> {
    match raw {
        None => Ok(DEFAULT_SUBSAMPLE_COUNT),
        Some(s) => {
            let trimmed = s.trim();
            if trimmed.is_empty() {
                Ok(DEFAULT_SUBSAMPLE_COUNT)
            } else {
                trimmed.parse::<u32>().map_err(|_| {
                    AppError::Config(format!(
                        "row {row}: subsample_count '{trimmed}' is not a valid non-negative integer"
                    ))
                })
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_csv(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    // --- parse_count ---

    #[test]
    fn parse_count_none_is_default() {
        assert_eq!(parse_count(&None, 2).unwrap(), DEFAULT_SUBSAMPLE_COUNT);
    }

    #[test]
    fn parse_count_empty_is_default() {
        assert_eq!(
            parse_count(&Some(String::new()), 2).unwrap(),
            DEFAULT_SUBSAMPLE_COUNT
        );
        assert_eq!(
            parse_count(&Some("   ".into()), 2).unwrap(),
            DEFAULT_SUBSAMPLE_COUNT
        );
    }

    #[test]
    fn parse_count_valid_number() {
        assert_eq!(parse_count(&Some("5000".into()), 2).unwrap(), 5000);
        assert_eq!(parse_count(&Some(" 42 ".into()), 3).unwrap(), 42);
    }

    #[test]
    fn parse_count_rejects_garbage_and_negative() {
        assert!(parse_count(&Some("abc".into()), 2).is_err());
        assert!(parse_count(&Some("-5".into()), 2).is_err());
        assert!(parse_count(&Some("3.5".into()), 2).is_err());
    }

    // --- load_config_csv ---

    #[test]
    fn parses_csv_with_counts_and_blank_rows() {
        let f = write_csv("seq_name,subsample_count\nchr1,5000\nchr2,\nchrX,2500\n");
        let map = load_config_csv(f.path()).unwrap();
        assert_eq!(map.len(), 3);
        assert_eq!(map.get("chr1"), Some(&5000));
        assert_eq!(map.get("chr2"), Some(&DEFAULT_SUBSAMPLE_COUNT));
        assert_eq!(map.get("chrX"), Some(&2500));
    }

    #[test]
    fn duplicate_seq_name_last_wins() {
        let f = write_csv("seq_name,subsample_count\nchr1,100\nchr1,200\n");
        let map = load_config_csv(f.path()).unwrap();
        assert_eq!(map.get("chr1"), Some(&200));
    }

    #[test]
    fn malformed_count_returns_config_error() {
        let f = write_csv("seq_name,subsample_count\nchr1,notanumber\n");
        let err = load_config_csv(f.path()).unwrap_err();
        let msg = err.to_string();
        assert!(msg.starts_with("config error:"), "got: {msg}");
        assert!(msg.contains("row 2"), "got: {msg}");
        assert!(msg.contains("notanumber"), "got: {msg}");
    }

    #[test]
    fn missing_file_returns_error() {
        let err = load_config_csv(Path::new("/nonexistent/does-not-exist.csv")).unwrap_err();
        assert!(err.to_string().starts_with("CSV error:") || err.to_string().contains("I/O"));
    }

    // --- SubsamplePlan::count_for ---

    #[test]
    fn plan_global_applies_to_every_ref() {
        let plan = SubsamplePlan::Global(7);
        assert_eq!(plan.count_for("chr1"), 7);
        assert_eq!(plan.count_for("chr2"), 7);
    }

    #[test]
    fn plan_perref_listed_and_unlisted() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), 100u32);
        let plan = SubsamplePlan::PerRef(map);
        assert_eq!(plan.count_for("chr1"), 100);
        assert_eq!(plan.count_for("chr9"), DEFAULT_SUBSAMPLE_COUNT as usize);
    }

    #[test]
    fn plan_default_uses_default_everywhere() {
        let plan = SubsamplePlan::Default;
        assert_eq!(plan.count_for("anything"), DEFAULT_SUBSAMPLE_COUNT as usize);
    }
}
