//! Typed errors for sam-subsampler.
//!
//! Every error is explicit and propagated with `?`; nothing is swallowed.
//! `#[from]` wrappers carry the underlying cause; `Argument` and `Config`
//! hold human-readable messages produced at the system boundary
//! (CLI validation / config parsing).

use thiserror::Error;

/// All errors produced by sam-subsampler.
#[derive(Debug, Error)]
pub enum AppError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),

    #[error("BAM/HTSlib error: {0}")]
    Htslib(#[from] rust_htslib::errors::Error),

    #[error("UTF-8 conversion error: {0}")]
    Utf8(#[from] std::string::FromUtf8Error),

    #[error("argument error: {0}")]
    Argument(String),

    #[error("config error: {0}")]
    Config(String),
}

/// Convenience alias used throughout the crate.
pub type Result<T> = std::result::Result<T, AppError>;

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn io_error_wraps_underlying_message() {
        let err = AppError::from(std::io::Error::new(std::io::ErrorKind::NotFound, "missing"));
        let msg = err.to_string();
        assert!(msg.starts_with("I/O error:"), "got: {msg}");
        assert!(msg.contains("missing"), "got: {msg}");
    }

    #[test]
    fn argument_error_displays_message() {
        let err = AppError::Argument("tag must be 2 chars".into());
        assert_eq!(err.to_string(), "argument error: tag must be 2 chars");
    }

    #[test]
    fn config_error_displays_message() {
        let err = AppError::Config("duplicate seq_name 'chr1'".into());
        assert_eq!(err.to_string(), "config error: duplicate seq_name 'chr1'");
    }

    #[test]
    fn utf8_error_propagates_from_from_utf8() {
        let bad = vec![0xFFu8, 0xFE, 0xFD];
        let src = String::from_utf8(bad).unwrap_err();
        let err = AppError::from(src);
        let msg = err.to_string();
        assert!(msg.starts_with("UTF-8 conversion error:"), "got: {msg}");
    }

    #[test]
    fn result_alias_propagates_errors() {
        fn fallible() -> Result<()> {
            Err(AppError::Argument("boom".into()))
        }
        let err = fallible().unwrap_err();
        assert_eq!(err.to_string(), "argument error: boom");
    }
}
