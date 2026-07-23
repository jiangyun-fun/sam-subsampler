//! Pure per-reference reservoir sampling.
//!
//! No BAM I/O here — everything is deterministic and unit-testable.
//!
//! [`select_per_reference`] expects the caller (pass 1 of the pipeline) to
//! have deduplicated qnames per reference, so each read name is one selection
//! unit regardless of how many records (mate, secondary/supplementary
//! alignments) carry it. This corrects the per-record selection bias of the
//! original `bam_subsampler`, where a paired read had roughly twice the
//! probability of being chosen.

use crate::config::SubsamplePlan;
use rand::prelude::*;
use rand::rngs::StdRng;
use std::collections::HashSet;

/// Vitter's reservoir sampling — Algorithm R.
///
/// Returns up to `count` items chosen uniformly at random without replacement.
/// If `items.len() <= count`, all items are returned unchanged.
///
/// The output depends on the *order* of `items`; callers that need
/// cross-process reproducibility should sort the input first (see
/// [`select_per_reference`]).
pub fn reservoir_select<T>(mut items: Vec<T>, count: usize, rng: &mut impl Rng) -> Vec<T> {
    if items.len() <= count {
        return items;
    }
    for i in count..items.len() {
        // Uniform integer in the closed range [0, i].
        let j = rng.random_range(0..=i);
        if j < count {
            items.swap(j, i);
        }
    }
    items.truncate(count);
    items
}

/// Select the set of qnames to tag, dispatching on the [`SubsamplePlan`].
///
/// Single public entry point for selection. Per-reference plans
/// ([`SubsamplePlan::Global`]/[`PerRef`]/[`Default`]) draw one reservoir per
/// reference via [`select_per_reference`]. The reference-agnostic global plans
/// ([`SubsamplePlan::GlobalTotal`]/[`GlobalRatio`]) pool every unique qname
/// across references into one reservoir, ignoring which reference each read
/// mapped to.
///
/// The result is a pure function of (input, plan, seed): references are
/// processed in sorted order, pooled qnames are sorted, and a single RNG seeded
/// once drives every draw.
pub fn select(
    qnames_by_ref: crate::QnamesByRef,
    plan: &SubsamplePlan,
    seed: u64,
) -> HashSet<Vec<u8>> {
    match plan {
        SubsamplePlan::GlobalTotal(n) => {
            let pool = pooled_sorted_unique(qnames_by_ref);
            global_reservoir(pool, *n as usize, seed)
        }
        SubsamplePlan::GlobalRatio(r) => {
            let pool = pooled_sorted_unique(qnames_by_ref);
            // Base the ratio on the true deduped pool length, not the sum of
            // per-reference set sizes (which double-counts qnames that appear
            // on more than one reference).
            let target = ratio_to_target(pool.len(), *r);
            global_reservoir(pool, target, seed)
        }
        _ => select_per_reference(qnames_by_ref, plan, seed),
    }
}

/// Flatten every per-reference unique-qname set into one deduplicated, sorted
/// pool — the input to a single global reservoir draw.
///
/// Dedup is *across references* (a read whose records span multiple references
/// is one selection unit); sorting makes the draw reproducible despite
/// `HashSet`'s randomized iteration order.
fn pooled_sorted_unique(qnames_by_ref: crate::QnamesByRef) -> Vec<Vec<u8>> {
    let pooled: HashSet<Vec<u8>> = qnames_by_ref.into_values().flatten().collect();
    let mut uniq: Vec<Vec<u8>> = pooled.into_iter().collect();
    uniq.sort_unstable();
    uniq
}

/// Draw `target` qnames from a sorted pool with a fresh seeded RNG.
fn global_reservoir(pool: Vec<Vec<u8>>, target: usize, seed: u64) -> HashSet<Vec<u8>> {
    let mut rng = StdRng::seed_from_u64(seed);
    reservoir_select(pool, target, &mut rng)
        .into_iter()
        .collect()
}

/// Absolute target for a ratio plan: `round(total · ratio)` clamped to
/// `[1, total]` when `total > 0`; `0` when nothing is available to sample.
///
/// `ratio` is assumed finite and in `(0.0, 1.0]` (enforced by the CLI).
fn ratio_to_target(total: usize, ratio: f64) -> usize {
    if total == 0 {
        return 0;
    }
    let raw = (total as f64 * ratio).round() as usize;
    raw.clamp(1, total)
}

/// Select the set of qnames to tag — one reservoir per reference.
///
/// `qnames_by_ref` must already hold unique qnames per reference (dedup is the
/// caller's responsibility — see the module docs). Within each reference the
/// qnames are sorted before sampling so the result is reproducible across
/// process runs despite `HashSet`'s randomized iteration order. A single RNG,
/// seeded once, is drawn from for every reference (processed in sorted order),
/// which makes the full selection a pure function of (input, plan, seed).
pub fn select_per_reference(
    qnames_by_ref: crate::QnamesByRef,
    plan: &SubsamplePlan,
    seed: u64,
) -> HashSet<Vec<u8>> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut refs: Vec<String> = qnames_by_ref.keys().cloned().collect();
    refs.sort_unstable();

    let mut selected: HashSet<Vec<u8>> = HashSet::new();
    for name in refs {
        let target = plan.count_for(&name);
        let mut uniq: Vec<Vec<u8>> = qnames_by_ref
            .get(&name)
            .map(|set| set.iter().cloned().collect())
            .unwrap_or_default();
        uniq.sort_unstable();
        for q in reservoir_select(uniq, target, &mut rng) {
            selected.insert(q);
        }
    }
    selected
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use std::collections::HashMap;

    fn rng_from(seed: u64) -> StdRng {
        StdRng::seed_from_u64(seed)
    }

    // --- reservoir_select ---

    #[test]
    fn reservoir_returns_all_when_count_at_least_len() {
        let mut rng = rng_from(1);
        let out = reservoir_select(vec![10, 20, 30], 3, &mut rng);
        assert_eq!(out.len(), 3);
        // count > len
        let out = reservoir_select(vec![10, 20, 30], 10, &mut rng);
        assert_eq!(out.len(), 3);
    }

    #[test]
    fn reservoir_returns_exactly_count_items() {
        let mut rng = rng_from(7);
        let out = reservoir_select((0..100).collect::<Vec<_>>(), 5, &mut rng);
        assert_eq!(out.len(), 5);
    }

    #[test]
    fn reservoir_output_is_subset_of_input() {
        let mut rng = rng_from(42);
        let input: Vec<i32> = (0..50).collect();
        let out = reservoir_select(input.clone(), 8, &mut rng);
        assert!(out.iter().all(|x| input.contains(x)));
    }

    #[test]
    fn reservoir_output_has_no_duplicates() {
        let mut rng = rng_from(99);
        let out = reservoir_select((0..200).collect::<Vec<_>>(), 15, &mut rng);
        let unique: HashSet<_> = out.iter().collect();
        assert_eq!(unique.len(), out.len(), "reservoir produced duplicates");
    }

    #[test]
    fn reservoir_is_deterministic_for_same_seed() {
        let a = reservoir_select((0..1000).collect::<Vec<_>>(), 50, &mut rng_from(123));
        let b = reservoir_select((0..1000).collect::<Vec<_>>(), 50, &mut rng_from(123));
        assert_eq!(a, b);
    }

    #[test]
    fn reservoir_different_seeds_usually_differ() {
        // 50 of 1000 chosen; identical sets across 100 independent draws would
        // be astronomically improbable, so a single comparison is sufficient.
        let a = reservoir_select((0..1000).collect::<Vec<_>>(), 50, &mut rng_from(1));
        let b = reservoir_select((0..1000).collect::<Vec<_>>(), 50, &mut rng_from(2));
        assert_ne!(a, b);
    }

    #[test]
    fn reservoir_distribution_is_roughly_uniform() {
        // Every item in a small universe should be selected at least once over
        // many draws (a basic sanity check on uniformity).
        let mut counts = vec![0u32; 20];
        for seed in 0..400u64 {
            let out = reservoir_select((0..20).collect::<Vec<_>>(), 5, &mut rng_from(seed));
            for x in out {
                counts[x as usize] += 1;
            }
        }
        assert!(
            counts.iter().all(|&c| c > 0),
            "some items never selected: {counts:?}"
        );
    }

    // --- select_per_reference ---

    fn set_of(qnames: &[&[u8]]) -> HashSet<Vec<u8>> {
        qnames.iter().map(|q| q.to_vec()).collect()
    }

    #[test]
    fn per_ref_respects_global_count() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b", b"c", b"d", b"e"]));
        map.insert("chr2".into(), set_of(&[b"x", b"y", b"z"]));

        let plan = SubsamplePlan::Global(2);
        let selected = select_per_reference(map, &plan, 42);

        // chr1 contributes 2, chr2 has only 3 so contributes 2 → 4 total.
        assert_eq!(selected.len(), 4);
    }

    #[test]
    fn per_ref_respects_per_ref_plan() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b", b"c", b"d", b"e"]));
        map.insert("chr2".into(), set_of(&[b"w", b"x", b"y", b"z"]));

        let mut cfg = HashMap::new();
        cfg.insert("chr1".into(), 3u32);
        // chr2 unlisted → default (1000) → all 4 taken.
        let plan = SubsamplePlan::PerRef(cfg);

        let selected = select_per_reference(map, &plan, 42);
        assert_eq!(selected.len(), 7); // 3 + 4
    }

    #[test]
    fn per_ref_takes_all_when_fewer_than_target() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"only"])); // 1 unique
        let selected = select_per_reference(map, &SubsamplePlan::Global(100), 1);
        assert_eq!(selected.len(), 1);
    }

    #[test]
    fn per_ref_is_reproducible_regardless_of_hashmap_construction_order() {
        let mk = || {
            let mut m = HashMap::new();
            m.insert(
                "chr1".into(),
                set_of(&[b"a", b"b", b"c", b"d", b"e", b"f", b"g"]),
            );
            m
        };
        let a = select_per_reference(mk(), &SubsamplePlan::Global(3), 42);
        let b = select_per_reference(mk(), &SubsamplePlan::Global(3), 42);
        assert_eq!(a, b);
    }

    #[test]
    fn per_ref_unique_qname_is_one_selection_unit() {
        // Bias-fix guard: even though a read with many records is one unit,
        // the input set is unique by construction. A reference with N unique
        // qnames and target N-1 yields exactly N-1 selected qnames.
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"r1", b"r2", b"r3", b"r4", b"r5"]));
        let selected = select_per_reference(map, &SubsamplePlan::Global(4), 42);
        assert_eq!(selected.len(), 4);
    }

    #[test]
    fn per_ref_different_seed_likely_different_set() {
        let mk = || {
            let mut m = HashMap::new();
            m.insert(
                "chr1".into(),
                set_of(&[b"a", b"b", b"c", b"d", b"e", b"f", b"g", b"h", b"i", b"j"]),
            );
            m
        };
        let a = select_per_reference(mk(), &SubsamplePlan::Global(3), 1);
        let b = select_per_reference(mk(), &SubsamplePlan::Global(3), 2);
        assert_ne!(a, b);
    }

    // --- select (global dispatch) ---

    #[test]
    fn global_total_pools_across_references() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b", b"c", b"d", b"e"])); // 5
        map.insert("chr2".into(), set_of(&[b"w", b"x", b"y", b"z"])); // 4
        // 9 unique total; GlobalTotal(3) picks 3 regardless of reference.
        let selected = select(map, &SubsamplePlan::GlobalTotal(3), 42);
        assert_eq!(selected.len(), 3);
    }

    #[test]
    fn global_total_dedups_cross_reference_qname() {
        // "shared" appears on both references but is one selection unit.
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"shared", b"only1"]));
        map.insert("chr2".into(), set_of(&[b"shared", b"only2"]));
        // 3 unique total (shared, only1, only2); GlobalTotal(large) takes all 3.
        let selected = select(map, &SubsamplePlan::GlobalTotal(100), 7);
        assert_eq!(selected.len(), 3);
        assert!(selected.contains(b"shared".as_slice()));
    }

    #[test]
    fn global_total_zero_selects_nothing() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b"]));
        assert_eq!(select(map, &SubsamplePlan::GlobalTotal(0), 42).len(), 0);
    }

    #[test]
    fn global_ratio_derives_count() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b", b"c", b"d", b"e"])); // 5
        map.insert("chr2".into(), set_of(&[b"f", b"g", b"h", b"i", b"j"])); // 5
        // 10 unique; ratio 0.5 -> 5.
        assert_eq!(select(map, &SubsamplePlan::GlobalRatio(0.5), 42).len(), 5);
    }

    #[test]
    fn global_ratio_rounds_and_clamps_to_one() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b", b"c"])); // 3 unique
        // round(3 * 0.1) = round(0.3) = 0, clamped to [1, 3] -> 1.
        assert_eq!(select(map, &SubsamplePlan::GlobalRatio(0.1), 42).len(), 1);
    }

    #[test]
    fn global_ratio_one_takes_all() {
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b"]));
        map.insert("chr2".into(), set_of(&[b"c"]));
        // ratio 1.0 -> round(3 * 1.0) = 3 (all).
        assert_eq!(select(map, &SubsamplePlan::GlobalRatio(1.0), 42).len(), 3);
    }

    #[test]
    fn global_is_reproducible_per_seed() {
        // A large pool split across two references (1000 unique qnames). Selecting
        // 50 of 1000 makes a cross-seed collision astronomically unlikely, so a
        // single `assert_ne` is sufficient (mirrors reservoir_different_seeds_*).
        let mk = || {
            let mut m = HashMap::new();
            let chr1: HashSet<Vec<u8>> = (0..500).map(|i| format!("a{i}").into_bytes()).collect();
            let chr2: HashSet<Vec<u8>> =
                (500..1000).map(|i| format!("a{i}").into_bytes()).collect();
            m.insert("chr1".into(), chr1);
            m.insert("chr2".into(), chr2);
            m
        };
        let a = select(mk(), &SubsamplePlan::GlobalTotal(50), 1);
        let b = select(mk(), &SubsamplePlan::GlobalTotal(50), 1);
        assert_eq!(a, b);
        // A different seed very likely yields a different set of the same size.
        let c = select(mk(), &SubsamplePlan::GlobalTotal(50), 2);
        assert_eq!(c.len(), 50);
        assert_ne!(a, c);
    }

    #[test]
    fn select_dispatches_per_ref_for_global_count() {
        // Regression: SubsamplePlan::Global(n) still means per-reference, not pooled.
        let mut map = HashMap::new();
        map.insert("chr1".into(), set_of(&[b"a", b"b", b"c", b"d", b"e"])); // 5
        map.insert("chr2".into(), set_of(&[b"x", b"y", b"z"])); // 3
        // Global(2) -> 2 from chr1 + 2 from chr2 = 4 (NOT min(2, 8) = 2).
        assert_eq!(select(map, &SubsamplePlan::Global(2), 42).len(), 4);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;
    use std::collections::HashMap;

    proptest! {
        #[test]
        fn reservoir_output_always_subset_and_bounded(n in 0usize..200, k in 0usize..200, seed in 0u64..1000) {
            let input: Vec<usize> = (0..n).collect();
            let mut rng = StdRng::seed_from_u64(seed);
            let out = reservoir_select(input.clone(), k, &mut rng);
            prop_assert!(out.len() == n.min(k));
            let input_set: HashSet<usize> = input.into_iter().collect();
            prop_assert!(out.iter().all(|x| input_set.contains(x)));
            let unique: HashSet<usize> = out.iter().copied().collect();
            prop_assert_eq!(unique.len(), out.len());
        }

        #[test]
        fn global_ratio_output_is_bounded_subset(
            n in 0usize..100,
            ratio in 0.01f64..=1.0,
            seed in 0u64..1000,
        ) {
            let mut map: HashMap<String, HashSet<Vec<u8>>> = HashMap::new();
            let pool: Vec<Vec<u8>> = (0..n).map(|i| format!("q{i}").into_bytes()).collect();
            let universe: HashSet<Vec<u8>> = pool.iter().cloned().collect();
            if n > 0 {
                map.insert("chr1".into(), pool.iter().cloned().collect());
            }
            let out = select(map, &SubsamplePlan::GlobalRatio(ratio), seed);
            prop_assert!(out.iter().all(|q| universe.contains(q)));
            prop_assert_eq!(out.len(), ratio_to_target(n, ratio));
        }
    }
}
