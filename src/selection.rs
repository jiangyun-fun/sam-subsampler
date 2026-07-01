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
use std::collections::{HashMap, HashSet};

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

/// Select the set of qnames to tag — one reservoir per reference.
///
/// `qnames_by_ref` must already hold unique qnames per reference (dedup is the
/// caller's responsibility — see the module docs). Within each reference the
/// qnames are sorted before sampling so the result is reproducible across
/// process runs despite `HashSet`'s randomized iteration order. A single RNG,
/// seeded once, is drawn from for every reference (processed in sorted order),
/// which makes the full selection a pure function of (input, plan, seed).
pub fn select_per_reference(
    qnames_by_ref: HashMap<String, HashSet<Vec<u8>>>,
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
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

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
    }
}
