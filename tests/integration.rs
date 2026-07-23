//! End-to-end tests for the two-pass subsample-and-tag pipeline.
//!
//! A small BAM fixture is built by converting a SAM string to BAM with
//! rust-htslib itself (no external `samtools` dependency). The fixture
//! includes a paired read (shared qname on two records), an unmapped read, and
//! known per-reference counts, so it exercises the bias fix, the unmapped
//! skip, the per-reference reservoir, and reproducibility.

#![allow(clippy::unwrap_used)]

use rust_htslib::bam::{self, Format, Read};
use sam_subsampler::{bam_io, config::SubsamplePlan, selection};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// Fixture SAM. `pair1` appears twice (a proper pair) to test the bias fix.
const SAM: &str = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:1000
r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r2\t0\tchr1\t2\t60\t10M\t*\t0\t0\tCGTACGTACG\tIIIIIIIIII
r3\t0\tchr1\t3\t60\t10M\t*\t0\t0\tGTACGTACGT\tIIIIIIIIII
pair1\t99\tchr1\t4\t60\t10M\t=\t10\t200\tTACGTACGTA\tIIIIIIIIII
pair1\t147\tchr1\t10\t60\t10M\t=\t4\t-200\tACGTACGTAC\tIIIIIIIIII
r4\t0\tchr2\t1\t60\t10M\t*\t0\t0\tTTTTGGGGCC\tIIIIIIIIII
r5\t0\tchr2\t2\t60\t10M\t*\t0\t0\tGGGGCCTTTT\tIIIIIIIIII
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tNNNNNNNNNN\tIIIIIIIIII
";

/// Build a BAM file from `SAM` via rust-htslib (SAM reader → BAM writer).
/// Returns the temp dir (keep alive) and the BAM path.
fn write_bam_from_sam() -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let sam_path = dir.path().join("in.sam");
    let bam_path = dir.path().join("in.bam");
    std::fs::write(&sam_path, SAM).unwrap();

    let mut reader = bam::Reader::from_path(&sam_path).unwrap();
    let header = bam::Header::from_template(reader.header());
    {
        let mut writer = bam::Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
        for result in reader.records() {
            writer.write(&result.unwrap()).unwrap();
        }
    }
    (dir, bam_path)
}

/// Run the public pipeline (pass1 → select → pass2) with the given plan/seed.
fn run_pipeline(input: &Path, output: &Path, plan: SubsamplePlan, seed: u64) {
    let (qnames_by_ref, total) = bam_io::read_unique_qnames_by_ref(input).unwrap();
    let selected = selection::select(qnames_by_ref, &plan, seed);
    bam_io::tag_and_write(bam_io::TagWrite {
        input,
        output,
        output_format: Format::Bam,
        reference: None,
        selected: &selected,
        tag: b"YS",
        total_records: total,
        show_progress: false,
    })
    .unwrap();
}

fn count_records(bam_path: &Path) -> usize {
    bam::Reader::from_path(bam_path).unwrap().records().count()
}

/// Unique qnames carrying the `YS` aux tag.
fn tagged_qnames(bam_path: &Path) -> HashSet<Vec<u8>> {
    let mut reader = bam::Reader::from_path(bam_path).unwrap();
    let mut out = HashSet::new();
    for result in reader.records() {
        let rec = result.unwrap();
        if has_ys(&rec) {
            out.insert(rec.qname().to_vec());
        }
    }
    out
}

/// Number of *records* carrying the `YS` tag (count of lines, not unique qnames).
fn tagged_record_count(bam_path: &Path) -> usize {
    let mut reader = bam::Reader::from_path(bam_path).unwrap();
    reader
        .records()
        .filter(|r| has_ys(r.as_ref().unwrap()))
        .count()
}

/// True when `rec` carries any value under the `YS` tag.
fn has_ys(rec: &bam::Record) -> bool {
    rec.aux(b"YS").is_ok()
}

#[test]
fn tagging_preserves_record_count() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    let n_in = count_records(&bam_path);
    run_pipeline(&bam_path, &out, SubsamplePlan::Global(2), 42);
    assert_eq!(count_records(&out), n_in, "tagging must not drop records");
}

#[test]
fn selects_exactly_count_unique_qnames_per_ref() {
    // chr1: 4 unique qnames (r1,r2,r3,pair1); chr2: 2 unique (r4,r5).
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(&bam_path, &out, SubsamplePlan::Global(2), 42);
    let tagged = tagged_qnames(&out);
    // chr1 -> 2, chr2 -> 2 (only 2 available) => 4 unique qnames tagged.
    assert_eq!(
        tagged.len(),
        4,
        "tagged qnames: {:?}",
        tagged
            .iter()
            .map(|q| String::from_utf8_lossy(q))
            .collect::<Vec<_>>()
    );
}

#[test]
fn unmapped_read_is_never_tagged() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(&bam_path, &out, SubsamplePlan::Global(1000), 42);
    let tagged = tagged_qnames(&out);
    assert!(!tagged.contains(b"unmapped".as_slice()));
}

#[test]
fn paired_read_is_one_selection_unit_and_tags_both_mates() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    // Global(1000) selects every unique qname on every ref.
    run_pipeline(&bam_path, &out, SubsamplePlan::Global(1000), 42);

    // pair1 is one selection unit, so it is selected exactly once as a qname…
    let mut reader = bam::Reader::from_path(&out).unwrap();
    let pair1_tagged_records: usize = reader
        .records()
        .filter(|r| {
            let r = r.as_ref().unwrap();
            r.qname() == b"pair1" && has_ys(r)
        })
        .count();
    // …but both of its records must carry the tag (pair-preserving bias fix).
    assert_eq!(
        pair1_tagged_records, 2,
        "both mates of pair1 must be tagged"
    );
}

#[test]
fn unselected_reads_carry_no_tag() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    let total = count_records(&bam_path);
    // chr1 has 4 unique qnames; Global(2) selects only 2, leaving reads untagged.
    run_pipeline(&bam_path, &out, SubsamplePlan::Global(2), 42);
    assert!(
        tagged_record_count(&out) < total,
        "some records must be untagged"
    );
}

#[test]
fn selected_tag_value_is_i32_one() {
    use rust_htslib::bam::record::Aux;
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(&bam_path, &out, SubsamplePlan::Global(2), 42);
    let mut reader = bam::Reader::from_path(&out).unwrap();
    let mut found = false;
    for result in reader.records() {
        if let Ok(Aux::I32(1)) = result.unwrap().aux(b"YS") {
            found = true;
            break;
        }
    }
    assert!(found, "tagged records must carry Aux::I32(1)");
}

#[test]
fn same_seed_reproduces_identical_selected_set() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out1 = bam_path.with_file_name("a.bam");
    let out2 = bam_path.with_file_name("b.bam");
    run_pipeline(&bam_path, &out1, SubsamplePlan::Global(2), 42);
    run_pipeline(&bam_path, &out2, SubsamplePlan::Global(2), 42);
    assert_eq!(tagged_qnames(&out1), tagged_qnames(&out2));
}

// --- global mode (--total-count / --ratio): reference-agnostic selection ---
//
// Fixture has 6 unique mapped qnames: chr1 -> {r1, r2, r3, pair1}, chr2 ->
// {r4, r5}; "unmapped" is skipped. 8 records total (pair1 spans two).

fn dbg_qnames(set: &HashSet<Vec<u8>>) -> Vec<String> {
    let mut v: Vec<String> = set
        .iter()
        .map(|q| String::from_utf8_lossy(q).into_owned())
        .collect();
    v.sort();
    v
}

#[test]
fn global_total_tags_exactly_target_unique_qnames() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(&bam_path, &out, SubsamplePlan::GlobalTotal(3), 42);
    let tagged = tagged_qnames(&out);
    assert_eq!(tagged.len(), 3, "tagged qnames: {:?}", dbg_qnames(&tagged));
    // unmapped is never eligible, regardless of mode.
    assert!(!tagged.contains(b"unmapped".as_slice()));
}

#[test]
fn global_total_differs_from_per_reference_count() {
    // Same N, different semantics. Per-ref Global(3) -> 3 from chr1 + 2 from
    // chr2 = 5; GlobalTotal(3) -> exactly 3 (pooled, ignores reference).
    let (_dir, bam_path) = write_bam_from_sam();
    let per_ref_out = bam_path.with_file_name("per_ref.bam");
    let global_out = bam_path.with_file_name("global.bam");
    run_pipeline(&bam_path, &per_ref_out, SubsamplePlan::Global(3), 42);
    run_pipeline(&bam_path, &global_out, SubsamplePlan::GlobalTotal(3), 42);
    assert_eq!(tagged_qnames(&per_ref_out).len(), 5);
    assert_eq!(tagged_qnames(&global_out).len(), 3);
}

#[test]
fn global_total_preserves_record_count() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    let n_in = count_records(&bam_path);
    run_pipeline(&bam_path, &out, SubsamplePlan::GlobalTotal(3), 42);
    assert_eq!(count_records(&out), n_in, "tagging must not drop records");
}

#[test]
fn global_total_zero_tags_nothing() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(&bam_path, &out, SubsamplePlan::GlobalTotal(0), 42);
    assert_eq!(tagged_qnames(&out).len(), 0);
    assert_eq!(tagged_record_count(&out), 0);
}

#[test]
fn global_ratio_tags_rounded_count() {
    // round(6 * 0.5) = 3; round(6 * 0.34) = round(2.04) = 2.
    let (_dir, bam_path) = write_bam_from_sam();
    let half = bam_path.with_file_name("half.bam");
    run_pipeline(&bam_path, &half, SubsamplePlan::GlobalRatio(0.5), 42);
    assert_eq!(tagged_qnames(&half).len(), 3);
    let thirdish = bam_path.with_file_name("thirdish.bam");
    run_pipeline(&bam_path, &thirdish, SubsamplePlan::GlobalRatio(0.34), 42);
    assert_eq!(tagged_qnames(&thirdish).len(), 2);
}

#[test]
fn global_mode_reproducible_across_runs() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out1 = bam_path.with_file_name("a.bam");
    let out2 = bam_path.with_file_name("b.bam");
    run_pipeline(&bam_path, &out1, SubsamplePlan::GlobalTotal(3), 42);
    run_pipeline(&bam_path, &out2, SubsamplePlan::GlobalTotal(3), 42);
    assert_eq!(tagged_qnames(&out1), tagged_qnames(&out2));
}

#[test]
fn global_mode_selects_all_and_tags_both_pair_mates() {
    // GlobalTotal(large) selects every unique qname; pair1 must then have both
    // records tagged (the qname-dedup bias fix holds in global mode too).
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(&bam_path, &out, SubsamplePlan::GlobalTotal(1000), 42);
    assert_eq!(tagged_qnames(&out).len(), 6, "all 6 unique qnames selected");
    let mut reader = bam::Reader::from_path(&out).unwrap();
    let pair1_tagged: usize = reader
        .records()
        .filter(|r| {
            let r = r.as_ref().unwrap();
            r.qname() == b"pair1" && has_ys(r)
        })
        .count();
    assert_eq!(pair1_tagged, 2, "both mates of pair1 must be tagged");
}
