//! End-to-end tests for the two-pass subsample pipeline.
//!
//! A small BAM fixture is built by converting a SAM string to BAM with
//! rust-htslib itself (no external `samtools` dependency). The fixture
//! includes a paired read (shared qname on two records), three unmapped reads,
//! and known per-reference counts, so it exercises the bias fix, unmapped
//! subsampling (the `*` bucket), the per-reference reservoir, the keep-selected
//! vs keep-all output modes, and reproducibility.

#![allow(clippy::unwrap_used)]

use rust_htslib::bam::{self, Format, Read};
use sam_subsampler::{bam_io, config, config::SubsamplePlan, selection};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// Fixture SAM. `pair1` appears twice (a proper pair) to test the bias fix.
/// `unmapped`, `um2`, `um3` are unmapped reads → the `*` bucket.
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
um2\t4\t*\t0\t0\t*\t*\t0\t0\tNNNNNNNNNN\tIIIIIIIIII
um3\t4\t*\t0\t0\t*\t*\t0\t0\tNNNNNNNNNN\tIIIIIIIIII
";

/// Build a BAM from an arbitrary SAM string (SAM reader → BAM writer).
/// Returns the temp dir (keep alive) and the BAM path.
fn write_bam_from_sam_str(sam: &str) -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let sam_path = dir.path().join("in.sam");
    let bam_path = dir.path().join("in.bam");
    std::fs::write(&sam_path, sam).unwrap();

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

/// Build a BAM from the default fixture `SAM`.
fn write_bam_from_sam() -> (tempfile::TempDir, PathBuf) {
    write_bam_from_sam_str(SAM)
}

/// Run the pipeline (pass1 → select → pass2) with explicit format + reference.
fn run_pipeline_full(
    input: &Path,
    output: &Path,
    plan: SubsamplePlan,
    seed: u64,
    mode: bam_io::OutputMode,
    fmt: Format,
    reference: Option<&Path>,
) {
    let (qnames_by_ref, total) = bam_io::read_unique_qnames_by_ref(input).unwrap();
    let selected = selection::select(qnames_by_ref, &plan, seed);
    bam_io::tag_and_write(bam_io::TagWrite {
        input,
        output,
        output_format: fmt,
        reference,
        selected: &selected,
        tag: b"YS",
        total_records: total,
        mode,
        show_progress: false,
    })
    .unwrap();
}

/// Run the pipeline to a BAM output with no reference — the common case.
fn run_pipeline(
    input: &Path,
    output: &Path,
    plan: SubsamplePlan,
    seed: u64,
    mode: bam_io::OutputMode,
) {
    run_pipeline_full(input, output, plan, seed, mode, Format::Bam, None);
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

/// True when `q` is one of the fixture's unmapped read names.
fn is_unmapped_qname(q: &[u8]) -> bool {
    matches!(q, b"unmapped" | b"um2" | b"um3")
}

fn dbg_qnames(set: &HashSet<Vec<u8>>) -> Vec<String> {
    let mut v: Vec<String> = set
        .iter()
        .map(|q| String::from_utf8_lossy(q).into_owned())
        .collect();
    v.sort();
    v
}

// --- output modes: KeepSelected (default) vs TagInPlace (--keep-all) ---

#[test]
fn keep_selected_drops_unselected_records() {
    // Default mode writes only selected reads → output is smaller than input.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    let n_in = count_records(&bam_path);
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let n_out = count_records(&out);
    assert!(
        n_out < n_in,
        "keep-selected must drop records: {n_out} vs {n_in}"
    );
    assert!(n_out > 0, "some records must be kept");
}

#[test]
fn all_output_records_tagged_in_keep_selected() {
    // In keep-selected mode every written record was selected, so every record
    // carries the tag.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let n_out = count_records(&out);
    assert!(n_out > 0);
    assert_eq!(
        tagged_record_count(&out),
        n_out,
        "every kept record must carry the tag"
    );
}

#[test]
fn keep_all_preserves_record_count() {
    // --keep-all restores the pre-0.3 tag-in-place behavior: every record kept.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    let n_in = count_records(&bam_path);
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::TagInPlace,
    );
    assert_eq!(
        count_records(&out),
        n_in,
        "--keep-all must preserve every record"
    );
}

#[test]
fn unselected_reads_carry_no_tag_under_keep_all() {
    // Under --keep-all the unselected reads are still present but untagged.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    let total = count_records(&bam_path);
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::TagInPlace,
    );
    assert!(
        tagged_record_count(&out) < total,
        "under --keep-all some records must be untagged"
    );
}

// --- per-reference selection (Global / PerRef / Default) ---
//
// Fixture buckets: chr1 -> {r1,r2,r3,pair1} (4), chr2 -> {r4,r5} (2),
// '*' -> {unmapped,um2,um3} (3); 9 unique qnames, 10 records.

#[test]
fn selects_exactly_count_unique_qnames_per_ref() {
    // Global(2): chr1 -> 2, chr2 -> 2, '*' -> 2 => 6 unique qnames tagged.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let tagged = tagged_qnames(&out);
    assert_eq!(tagged.len(), 6, "tagged qnames: {:?}", dbg_qnames(&tagged));
}

#[test]
fn unmapped_reads_subsampled_under_count() {
    // The '*' bucket is sampled like a reference: Global(2) keeps 2 of 3
    // unmapped reads.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let tagged = tagged_qnames(&out);
    let unmapped_tagged = tagged.iter().filter(|q| is_unmapped_qname(q)).count();
    assert_eq!(unmapped_tagged, 2, "expected 2 of 3 unmapped tagged");
}

#[test]
fn star_bucket_in_config_controls_unmapped() {
    // A `*,N` row in the config CSV sets the unmapped subsample count.
    let (_dir, bam_path) = write_bam_from_sam();
    let csv = bam_path.with_file_name("refs.csv");
    std::fs::write(
        &csv,
        "seq_name,subsample_count\nchr1,1000\nchr2,1000\n*,1\n",
    )
    .unwrap();
    let map = config::load_config_csv(&csv).unwrap();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::PerRef(map),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let tagged = tagged_qnames(&out);
    // chr1 -> 4 (all), chr2 -> 2 (all), '*' -> 1 = 7; exactly 1 unmapped.
    assert_eq!(tagged.len(), 7, "tagged qnames: {:?}", dbg_qnames(&tagged));
    let unmapped_tagged = tagged.iter().filter(|q| is_unmapped_qname(q)).count();
    assert_eq!(unmapped_tagged, 1);
}

#[test]
fn paired_read_is_one_selection_unit_and_tags_both_mates() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    // Global(1000) selects every unique qname on every ref.
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(1000),
        42,
        bam_io::OutputMode::KeepSelected,
    );

    // pair1 is one selection unit, so it is selected exactly once as a qname…
    let mut reader = bam::Reader::from_path(&out).unwrap();
    let pair1_tagged_records: usize = reader
        .records()
        .filter(|r| {
            let r = r.as_ref().unwrap();
            r.qname() == b"pair1" && has_ys(r)
        })
        .count();
    // …but both of its records must be kept and tagged (pair-preserving bias fix).
    assert_eq!(pair1_tagged_records, 2, "both mates of pair1 must be kept");
}

#[test]
fn selected_tag_value_is_i32_one() {
    use rust_htslib::bam::record::Aux;
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
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
    run_pipeline(
        &bam_path,
        &out1,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    run_pipeline(
        &bam_path,
        &out2,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(tagged_qnames(&out1), tagged_qnames(&out2));
}

// --- global mode (--total-count / --ratio): reference-agnostic selection ---
//
// 9 unique qnames pooled across chr1, chr2, and the '*' unmapped bucket.

#[test]
fn global_total_tags_exactly_target_unique_qnames() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::GlobalTotal(3),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let tagged = tagged_qnames(&out);
    assert_eq!(tagged.len(), 3, "tagged qnames: {:?}", dbg_qnames(&tagged));
}

#[test]
fn global_total_differs_from_per_reference_count() {
    // Same N, different semantics. Per-ref Global(3) -> 3 chr1 + 2 chr2 + 3 '*'
    // = 8; GlobalTotal(3) -> exactly 3 (pooled, ignores reference).
    let (_dir, bam_path) = write_bam_from_sam();
    let per_ref_out = bam_path.with_file_name("per_ref.bam");
    let global_out = bam_path.with_file_name("global.bam");
    run_pipeline(
        &bam_path,
        &per_ref_out,
        SubsamplePlan::Global(3),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    run_pipeline(
        &bam_path,
        &global_out,
        SubsamplePlan::GlobalTotal(3),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(tagged_qnames(&per_ref_out).len(), 8);
    assert_eq!(tagged_qnames(&global_out).len(), 3);
}

#[test]
fn global_total_keep_all_preserves_record_count() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    let n_in = count_records(&bam_path);
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::GlobalTotal(3),
        42,
        bam_io::OutputMode::TagInPlace,
    );
    assert_eq!(
        count_records(&out),
        n_in,
        "--keep-all must preserve every record"
    );
}

#[test]
fn global_total_zero_tags_nothing() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::GlobalTotal(0),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(tagged_qnames(&out).len(), 0);
    assert_eq!(tagged_record_count(&out), 0);
    assert_eq!(count_records(&out), 0, "0 selected ⇒ empty output");
}

#[test]
fn global_ratio_tags_rounded_count() {
    // 9 unique pooled; round(9 * 0.5) = 5; round(9 * 0.34) = round(3.06) = 3.
    let (_dir, bam_path) = write_bam_from_sam();
    let half = bam_path.with_file_name("half.bam");
    run_pipeline(
        &bam_path,
        &half,
        SubsamplePlan::GlobalRatio(0.5),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(tagged_qnames(&half).len(), 5);
    let thirdish = bam_path.with_file_name("thirdish.bam");
    run_pipeline(
        &bam_path,
        &thirdish,
        SubsamplePlan::GlobalRatio(0.34),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(tagged_qnames(&thirdish).len(), 3);
}

#[test]
fn global_mode_reproducible_across_runs() {
    let (_dir, bam_path) = write_bam_from_sam();
    let out1 = bam_path.with_file_name("a.bam");
    let out2 = bam_path.with_file_name("b.bam");
    run_pipeline(
        &bam_path,
        &out1,
        SubsamplePlan::GlobalTotal(3),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    run_pipeline(
        &bam_path,
        &out2,
        SubsamplePlan::GlobalTotal(3),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(tagged_qnames(&out1), tagged_qnames(&out2));
}

#[test]
fn global_mode_selects_all_and_tags_both_pair_mates() {
    // GlobalTotal(large) selects every unique qname; pair1 must then have both
    // records kept (the qname-dedup bias fix holds in global mode too).
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::GlobalTotal(1000),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(
        tagged_qnames(&out).len(),
        9,
        "all 9 unique qnames selected (6 mapped + 3 unmapped)"
    );
    let mut reader = bam::Reader::from_path(&out).unwrap();
    let pair1_tagged: usize = reader
        .records()
        .filter(|r| {
            let r = r.as_ref().unwrap();
            r.qname() == b"pair1" && has_ys(r)
        })
        .count();
    assert_eq!(pair1_tagged, 2, "both mates of pair1 must be kept");
}

#[test]
fn global_mode_includes_unmapped_in_pool() {
    // The '*' bucket is part of the global pool: GlobalTotal(large) selects all
    // unmapped reads too.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::GlobalTotal(1000),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let tagged = tagged_qnames(&out);
    assert!(tagged.contains(b"unmapped".as_slice()));
    assert!(tagged.contains(b"um2".as_slice()));
    assert!(tagged.contains(b"um3".as_slice()));
}

// --- additional helpers for the complex scenarios ---

/// Write an indexed reference FASTA (single-line contigs) for CRAM tests.
/// Returns the FASTA path; `ref.fa.fai` is written beside it. Sequence content
/// is all `A` (only length matters for CRAM; mismatches are stored as diffs).
fn write_indexed_reference(dir: &Path, contigs: &[(&str, usize)]) -> PathBuf {
    let ref_path = dir.join("ref.fa");
    let mut bytes: Vec<u8> = Vec::new();
    let mut fai = String::new();
    let mut offset: u64 = 0;
    for (name, len) in contigs {
        let header = format!(">{name}\n");
        let seq_line = format!("{}\n", "A".repeat(*len));
        let seq_offset = offset + header.len() as u64;
        bytes.extend_from_slice(header.as_bytes());
        bytes.extend_from_slice(seq_line.as_bytes());
        fai.push_str(&format!(
            "{name}\t{len}\t{seq_offset}\t{len}\t{}\n",
            len + 1
        ));
        offset += (header.len() + seq_line.len()) as u64;
    }
    std::fs::write(&ref_path, &bytes).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), &fai).unwrap();
    ref_path
}

/// Build a BAM whose records are those of `src` in reversed order.
fn write_reversed_bam(src: &Path) -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("rev.bam");
    let mut reader = bam::Reader::from_path(src).unwrap();
    let header = bam::Header::from_template(reader.header());
    let mut recs: Vec<bam::Record> = reader.records().map(|r| r.unwrap()).collect();
    recs.reverse();
    let mut writer = bam::Writer::from_path(&out, &header, Format::Bam).unwrap();
    for rec in &recs {
        writer.write(rec).unwrap();
    }
    (dir, out)
}

/// True when `small` is a (non-contiguous) subsequence of `large` — order-aware.
fn is_subsequence<T: PartialEq>(small: &[T], large: &[T]) -> bool {
    let mut it = large.iter();
    small.iter().all(|s| it.any(|l| l == s))
}

// --- complex scenarios ---

#[test]
fn supplementary_and_secondary_kept_as_one_unit() {
    // A read with a primary, a supplementary (0x800) and a secondary (0x100)
    // alignment is ONE selection unit; when selected, ALL its records are kept
    // and tagged (the qname-dedup bias fix extended beyond mate pairs).
    let sam = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:1000
supread\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
supread\t2048\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII
supread\t256\tchr1\t200\t60\t3M\t*\t0\t0\tACG\tIII
other\t0\tchr1\t300\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII
";
    let (_dir, bam_path) = write_bam_from_sam_str(sam);
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(1000),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let mut reader = bam::Reader::from_path(&out).unwrap();
    let sup_tagged = reader
        .records()
        .filter(|r| {
            let r = r.as_ref().unwrap();
            r.qname() == b"supread" && has_ys(r)
        })
        .count();
    assert_eq!(sup_tagged, 3, "all 3 records of supread kept and tagged");
}

#[test]
fn cross_reference_supplementary_records_all_tagged() {
    // xread has a primary on chr1 and a supplementary on chr2. Global mode dedups
    // it to one selection unit; when selected, records on BOTH references are kept.
    let sam = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:1000
xread\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
xread\t2048\tchr2\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
";
    let (_dir, bam_path) = write_bam_from_sam_str(sam);
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::GlobalTotal(1000),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let tagged = tagged_qnames(&out);
    assert_eq!(tagged.len(), 1);
    assert!(tagged.contains(b"xread".as_slice()));
    let mut reader = bam::Reader::from_path(&out).unwrap();
    let xread_records = reader
        .records()
        .filter(|r| r.as_ref().unwrap().qname() == b"xread")
        .count();
    assert_eq!(
        xread_records, 2,
        "both chr1 primary and chr2 supplementary kept"
    );
}

#[test]
fn unmapped_placed_at_mate_position_goes_to_star_bucket() {
    // An unmapped read (flag 0x4) that still carries its mate's RNAME/POS must be
    // pooled under '*' (is_unmapped wins over tid), not the mate's reference.
    // Config chr1→0 keeps the mate out; *→1000 selects the unmapped read.
    let sam = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:1000
mate1\t0\tchr1\t1\t60\t10M\t=\t100\t200\tACGTACGTAC\tIIIIIIIIII
placed_unmapped\t4\tchr1\t100\t0\t*\t=\t1\t-200\tNNNNNNNNNN\tIIIIIIIIII
";
    let (_dir, bam_path) = write_bam_from_sam_str(sam);
    let csv = bam_path.with_file_name("refs.csv");
    std::fs::write(&csv, "seq_name,subsample_count\nchr1,0\n*,1000\n").unwrap();
    let map = config::load_config_csv(&csv).unwrap();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::PerRef(map),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let tagged = tagged_qnames(&out);
    assert!(
        tagged.contains(b"placed_unmapped".as_slice()),
        "unmapped-with-RNAME must land in '*' and be selected"
    );
    assert!(
        !tagged.contains(b"mate1".as_slice()),
        "chr1 count 0 ⇒ mapped mate not selected"
    );
}

#[test]
fn keep_selected_preserves_input_record_order() {
    // The output must be an order-preserving subsequence of the input.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let in_order: Vec<Vec<u8>> = bam::Reader::from_path(&bam_path)
        .unwrap()
        .records()
        .map(|r| r.unwrap().qname().to_vec())
        .collect();
    let out_order: Vec<Vec<u8>> = bam::Reader::from_path(&out)
        .unwrap()
        .records()
        .map(|r| r.unwrap().qname().to_vec())
        .collect();
    assert!(
        is_subsequence(&out_order, &in_order),
        "output order must follow input order"
    );
}

#[test]
fn selection_invariant_under_reversed_input_order() {
    // Selection is a pure function of (input SET, plan, seed): reversing the
    // input record order must yield an identical selected set.
    let (_dir, bam_path) = write_bam_from_sam();
    let (_rdir, rev_path) = write_reversed_bam(&bam_path);
    let out1 = bam_path.with_file_name("a.bam");
    let out2 = bam_path.with_file_name("b.bam");
    run_pipeline(
        &bam_path,
        &out1,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    run_pipeline(
        &rev_path,
        &out2,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(tagged_qnames(&out1), tagged_qnames(&out2));
}

#[test]
fn empty_input_produces_empty_output() {
    // A header-only file selects nothing and writes nothing; must not panic.
    let (_dir, bam_path) =
        write_bam_from_sam_str("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n");
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(1000),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    assert_eq!(count_records(&out), 0);
    assert_eq!(tagged_qnames(&out).len(), 0);
}

#[test]
fn header_preserved_after_filtering() {
    // Filtering does not trim the header: all @SQ targets survive.
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.bam");
    run_pipeline(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
    );
    let in_hdr = bam::Reader::from_path(&bam_path).unwrap();
    let out_hdr = bam::Reader::from_path(&out).unwrap();
    assert_eq!(
        in_hdr.header().target_count(),
        out_hdr.header().target_count(),
        "@SQ count must be preserved"
    );
}

#[test]
fn sam_output_round_trip() {
    // The pipeline can write SAM output (format selected by extension).
    let (_dir, bam_path) = write_bam_from_sam();
    let out = bam_path.with_file_name("out.sam");
    run_pipeline_full(
        &bam_path,
        &out,
        SubsamplePlan::Global(2),
        42,
        bam_io::OutputMode::KeepSelected,
        Format::Sam,
        None,
    );
    let n_in = count_records(&bam_path);
    let n_out = count_records(&out);
    assert!(n_out > 0, "sam output non-empty");
    assert!(n_out < n_in, "sam output subsampled");
}

#[test]
fn cram_output_round_trip_in_keep_selected_mode() {
    // CRAM path: set_reference on both write and read-back; the keep-selected
    // mode and aux tagging survive a real CRAM encode/decode cycle.
    let sam = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:50
@SQ\tSN:chr2\tLN:50
cr1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
cr2\t0\tchr1\t2\t60\t10M\t*\t0\t0\tTTTTACGTAC\tIIIIIIIIII
cr3\t0\tchr2\t1\t60\t10M\t*\t0\t0\tGGGGACGTAC\tIIIIIIIIII
";
    let (dir, bam_path) = write_bam_from_sam_str(sam);
    let ref_path = write_indexed_reference(dir.path(), &[("chr1", 50), ("chr2", 50)]);
    let out = dir.path().join("out.cram");
    run_pipeline_full(
        &bam_path,
        &out,
        SubsamplePlan::Global(1000),
        42,
        bam_io::OutputMode::KeepSelected,
        Format::Cram,
        Some(&ref_path),
    );
    let mut reader = bam::Reader::from_path(&out).unwrap();
    reader.set_reference(&ref_path).unwrap();
    let recs: Vec<bam::Record> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(recs.len(), 3, "all 3 reads kept through CRAM round-trip");
    assert_eq!(
        recs.iter().filter(|r| r.aux(b"YS").is_ok()).count(),
        3,
        "all kept CRAM records tagged"
    );
}

#[test]
fn reservoir_selection_is_uniform_over_seeds() {
    // 8 unique qnames on chr1; GlobalTotal(4) ⇒ each selected with P = 0.5.
    // Over many seeds every qname must be selected roughly equally (fairness).
    let sam = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:1000
r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r2\t0\tchr1\t2\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r3\t0\tchr1\t3\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r4\t0\tchr1\t4\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r5\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r6\t0\tchr1\t6\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r7\t0\tchr1\t7\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r8\t0\tchr1\t8\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
";
    let (_dir, bam_path) = write_bam_from_sam_str(sam);
    let (qnames_by_ref, _) = bam_io::read_unique_qnames_by_ref(&bam_path).unwrap();
    let names: Vec<Vec<u8>> = (1..=8).map(|i| format!("r{i}").into_bytes()).collect();
    let mut counts = [0u32; 8];
    for seed in 0..1000u64 {
        let selected =
            selection::select(qnames_by_ref.clone(), &SubsamplePlan::GlobalTotal(4), seed);
        for (i, name) in names.iter().enumerate() {
            if selected.contains(name) {
                counts[i] += 1;
            }
        }
    }
    // Expected 500 each (P = 0.5 × 1000 seeds); ±100 is ~6σ ⇒ not flaky.
    for (i, &c) in counts.iter().enumerate() {
        assert!(
            (400..=600).contains(&c),
            "r{} selected {c} times over 1000 seeds; expected ~500",
            i + 1
        );
    }
}
