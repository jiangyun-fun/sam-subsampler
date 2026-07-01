//! sam-subsampler: two-pass BAM/CRAM subsampler that tags selected reads in place.
//!
//! This crate is an idiomatic Rust rewrite of `bam_subsampler`. See `README.md`
//! for the usage and design notes (per-reference reservoir sampling, qname-dedup
//! bias fix, CRAM reference handling).

pub mod cli;
pub mod config;
pub mod error;
pub mod selection;

pub use error::{AppError, Result};
