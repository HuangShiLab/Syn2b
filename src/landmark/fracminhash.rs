//! FracMinHash landmark selection — the sketch-based alternative to enzyme digestion.
//!
//! # Why this exists
//!
//! Everything downstream of landmark extraction in [`crate::synteny::scoring`] —
//! the repeat filter, the shared-tag restriction, per-contig closure, the
//! direction-free adjacency sets, SCJ, the orientation channel,
//! `observable_fraction`, the >=2-landmark relocation rule, and the error model —
//! consumes only a list of `(canonical identity, position, contig, orientation)`.
//! It has never depended on those landmarks coming from a restriction digest. This
//! module supplies the same four things from a hash sketch instead, so the same
//! structural mathematics applies to either.
//!
//! # The selection rule
//!
//! FracMinHash keeps a k-mer when `h(canonical(kmer)) < u64::MAX / scale`. Two
//! properties matter here and neither is shared by minimizers:
//!
//! - **Context-free.** Whether a k-mer is kept depends on the k-mer alone, never on
//!   its neighbours. A substitution 10 bp away cannot change the decision. Minimizer
//!   selection is window-relative, so one edit re-selects a whole neighbourhood and
//!   landmark identity is not stable across genomes — which is exactly what an
//!   adjacency-based structural metric cannot tolerate.
//! - **Genome-independent.** The threshold is a fixed constant, unlike bottom-s
//!   MinHash where the cutoff is the s-th smallest hash *of that genome* and
//!   therefore moves with genome size and content.
//!
//! Both properties are what make the selected set of one genome comparable to the
//! selected set of another, which is the precondition for every metric downstream.
//!
//! # Density
//!
//! Expected spacing is `scale` bp, so `scale` is a continuous knob on landmark
//! count `m` — where the enzyme path offers only the discrete steps of a 1-, 2-, 4-
//! or 16-enzyme panel. The error model `Var(err) = 1.504*p(1-p)/m + 0.0205^2` is a
//! function of `m`, so this is the only landmark source that can sweep it with
//! everything else held fixed.
//!
//! # Measured against the enzyme panel on E. coli K-12
//!
//! | property | BcgI | 4-enzyme panel | FracMinHash |
//! |---|---|---|---|
//! | landmarks | 2,872 | 6,079 | 6,034 (k=31, s=750) |
//! | multi-copy instances | 2.19% | 3.13% | 2.55% |
//! | multi-copy families | 13 | 38 | 39 |
//! | **unique landmarks 1 substitution from a family** | **7 (0.249%)** | **20 (0.340%)** | **0 (0.000%)** |
//! | cross-species leakage to B. subtilis | 0.000% | — | 0.000% |
//!
//! The last row is the operationally important one, and note that it is *not* about
//! repeat content: FracMinHash carries just as many genuine multi-copy families,
//! because repeats are a property of the genome rather than of the selection rule.
//! The difference is that none of its unique landmarks sits one substitution away
//! from one.
//!
//! That is the `sub_2` mechanism (see `docs/MATH_REVIEW.md`). A multi-copy family is
//! dropped by the per-genome uniqueness filter; in a diverged genome, once enough of
//! its copies are destroyed the survivor becomes unique and collides with a locus
//! elsewhere, so the metric reads a landmark that teleported. Enzyme landmarks must
//! contain a recognition motif, so they occupy a far smaller region of sequence
//! space and near-collisions are correspondingly likelier. FracMinHash k-mers are
//! drawn from the whole space with no shared constraint, and the >=2-landmark
//! relocation rule that exists to reject this has nothing to reject here.

use crate::enzyme::enzyme::EnzymeType;
use crate::tgt::tag::{Strand, Tag};
use anyhow::{bail, Result};

/// Largest k that fits the 32-byte tag buffer as ASCII.
pub const MAX_K: usize = 32;

/// FracMinHash landmark selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FracMinHash {
    /// k-mer length, 1..=32.
    pub k: usize,
    /// Compression factor. One k-mer in `scale` is kept, so expected landmark
    /// spacing is `scale` bp.
    pub scale: u64,
}

impl FracMinHash {
    pub fn new(k: usize, scale: u64) -> Result<Self> {
        if k == 0 || k > MAX_K {
            bail!("FracMinHash k must be 1..={} (got {}), because the tag buffer holds {} ASCII bases", MAX_K, k, MAX_K);
        }
        if scale == 0 {
            bail!("FracMinHash scale must be >= 1 (got 0)");
        }
        Ok(Self { k, scale })
    }

    /// Keep a k-mer when its canonical hash falls below this.
    ///
    /// `u64::MAX / scale` rather than `(u64::MAX + 1) / scale`, so at `scale = 1`
    /// the threshold is `u64::MAX` and the single hash equal to it is dropped —
    /// one k-mer in 2^64, which is not worth a wider integer to recover.
    pub fn threshold(&self) -> u64 {
        u64::MAX / self.scale
    }

    /// Select landmarks from one contig.
    ///
    /// `offset` is added to every position, matching
    /// [`crate::enzyme::digest::digest_genome_contig`], so a multi-contig genome
    /// carries genome-wide coordinates.
    ///
    /// Runs of non-ACGT reset the rolling encoder: a k-mer overlapping an `N`
    /// has no defined canonical form and is skipped, exactly as the enzyme path
    /// skips windows failing `is_pure_atcg`.
    pub fn landmarks(&self, sequence: &[u8], contig_id: u16, offset: u64) -> Vec<Tag> {
        let k = self.k;
        if sequence.len() < k {
            return Vec::new();
        }
        let thresh = self.threshold();
        let kmask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        let top_shift = 2 * (k - 1);

        let mut out = Vec::new();
        let (mut fwd, mut rev) = (0u64, 0u64);
        let mut valid = 0usize; // consecutive ACGT bases ending at the current one

        for (i, &b) in sequence.iter().enumerate() {
            let Some(c) = base_code(b) else {
                valid = 0;
                fwd = 0;
                rev = 0;
                continue;
            };
            fwd = ((fwd << 2) | c) & kmask;
            rev = (rev >> 2) | ((3 - c) << top_shift);
            valid += 1;
            if valid < k {
                continue;
            }
            // Select on the canonical k-mer, never on the forward one: the same
            // locus must be chosen whichever strand the contig was deposited on.
            // This is what makes the sketch reverse-complement symmetric, and the
            // whole orientation channel depends on it.
            let canonical = fwd.min(rev);
            if splitmix64(canonical) >= thresh {
                continue;
            }
            let start = i + 1 - k;
            let mut seq = [0u8; 32];
            for (dst, &src) in seq[..k].iter_mut().zip(&sequence[start..start + k]) {
                // Uppercase: soft-masked assemblies must not compare unequal to
                // unmasked ones. The rolling encoder already accepts either case.
                *dst = src.to_ascii_uppercase();
            }
            // No recognition site exists, so `strand` records which strand carried
            // the canonical form — the same information the enzyme path records by
            // which pattern matched. Scoring derives orientation from the stored
            // sequence via `is_revcomp_of_canonical`, not from this field.
            let strand = if fwd <= rev { Strand::Forward } else { Strand::Reverse };
            out.push(Tag::new(
                seq,
                offset + start as u64,
                EnzymeType::FracMinHash,
                strand,
                contig_id,
            ));
        }
        out
    }
}

#[inline]
fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// splitmix64. Chosen over a keyed hash so a TGT written by one build is
/// comparable to one written by any other: the sketch must be reproducible across
/// machines and versions, which rules out anything seeded from the environment.
#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn revcomp(s: &str) -> String {
        s.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', other => other,
            })
            .collect()
    }

    /// Deterministic pseudo-random DNA, so the tests do not depend on a fixture.
    fn dna(n: usize, seed: u64) -> String {
        let mut x = seed | 1;
        (0..n)
            .map(|_| {
                x = splitmix64(x);
                b"ACGT"[(x >> 33) as usize % 4] as char
            })
            .collect()
    }

    #[test]
    fn rejects_impossible_parameters() {
        assert!(FracMinHash::new(0, 100).is_err());
        assert!(FracMinHash::new(33, 100).is_err(), "k must fit the 32-byte buffer");
        assert!(FracMinHash::new(31, 0).is_err());
        assert!(FracMinHash::new(32, 1).is_ok(), "k = 32 must not overflow the mask");
    }

    #[test]
    fn selection_is_reverse_complement_symmetric() {
        // The property the whole orientation channel rests on: a genome and its
        // reverse complement must select the same loci, so that an inverted segment
        // yields the same landmarks read from either strand.
        let seq = dna(20_000, 7);
        let rc = revcomp(&seq);
        let fmh = FracMinHash::new(21, 50).unwrap();

        let fwd = fmh.landmarks(seq.as_bytes(), 0, 0);
        let rev = fmh.landmarks(rc.as_bytes(), 0, 0);
        assert!(!fwd.is_empty(), "sanity: the sketch must select something");
        assert_eq!(fwd.len(), rev.len(), "same loci must be selected on both strands");

        let mut a: Vec<[u8; 32]> = fwd.iter().map(|t| t.canonical_sequence()).collect();
        let mut b: Vec<[u8; 32]> = rev.iter().map(|t| t.canonical_sequence()).collect();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "canonical identities must be identical, not merely equinumerous");

        // And the positions must mirror: a landmark at `p` in the forward sequence
        // starts at `len - k - p` in the reverse complement.
        let n = seq.len() as u64;
        let k = 21u64;
        let mut pa: Vec<u64> = fwd.iter().map(|t| n - k - t.position).collect();
        let mut pb: Vec<u64> = rev.iter().map(|t| t.position).collect();
        pa.sort_unstable();
        pb.sort_unstable();
        assert_eq!(pa, pb, "mirrored positions must match");
    }

    #[test]
    fn density_tracks_the_scale_factor() {
        // Expected spacing is `scale` bp. Checked loosely — this is a statistical
        // claim about a finite sequence, not an identity.
        let seq = dna(400_000, 11);
        for scale in [100u64, 500, 2000] {
            let fmh = FracMinHash::new(31, scale).unwrap();
            let n = fmh.landmarks(seq.as_bytes(), 0, 0).len() as f64;
            let expected = (seq.len() - 30) as f64 / scale as f64;
            let ratio = n / expected;
            assert!(
                (0.75..1.25).contains(&ratio),
                "scale {scale}: got {n} landmarks, expected ~{expected:.0} (ratio {ratio:.3})"
            );
        }
    }

    #[test]
    fn selection_is_context_free() {
        // The property minimizers lack. A substitution changes only the k-mers that
        // physically contain it; every landmark whose k-mer is untouched must
        // survive with its position intact. This is why adjacency is stable across
        // genomes, and it is the reason this source can be swapped in at all.
        let seq = dna(50_000, 3);
        let fmh = FracMinHash::new(21, 40).unwrap();
        let before = fmh.landmarks(seq.as_bytes(), 0, 0);

        let mut mutated: Vec<u8> = seq.clone().into_bytes();
        let site = 25_000usize;
        mutated[site] = if mutated[site] == b'A' { b'C' } else { b'A' };
        let after = fmh.landmarks(&mutated, 0, 0);

        // Landmarks that do not overlap the edited base must be bit-identical.
        let untouched = |t: &Tag| (t.position as usize + 21) <= site || t.position as usize > site;
        let keep_b: Vec<_> = before.iter().filter(|t| untouched(t)).collect();
        let keep_a: Vec<_> = after.iter().filter(|t| untouched(t)).collect();
        assert_eq!(
            keep_b.len(),
            keep_a.len(),
            "a substitution must not re-select landmarks that do not contain it"
        );
        for (x, y) in keep_b.iter().zip(keep_a.iter()) {
            assert_eq!(x.position, y.position);
            assert_eq!(x.sequence, y.sequence);
        }
    }

    #[test]
    fn ambiguous_bases_are_skipped_without_shifting_positions() {
        let mut seq = dna(10_000, 5).into_bytes();
        for b in seq.iter_mut().take(5_050).skip(5_000) {
            *b = b'N';
        }
        let fmh = FracMinHash::new(21, 40).unwrap();
        let tags = fmh.landmarks(&seq, 0, 0);
        for t in &tags {
            let s = t.position as usize;
            assert!(
                s + 21 <= 5_000 || s >= 5_050,
                "a landmark at {s} overlaps the N-run"
            );
            assert!(
                t.sequence[..21].iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')),
                "stored sequence must be pure ACGT"
            );
        }
        // The encoder must resume after the N-run, not stop at it.
        assert!(
            tags.iter().any(|t| t.position as usize >= 5_050),
            "selection must resume after ambiguous bases"
        );
    }

    #[test]
    fn positions_and_sequences_agree_with_the_source() {
        let seq = dna(30_000, 13);
        let fmh = FracMinHash::new(25, 60).unwrap();
        for t in fmh.landmarks(seq.as_bytes(), 3, 1_000_000) {
            let start = (t.position - 1_000_000) as usize;
            assert_eq!(t.contig_id, 3);
            assert_eq!(
                &t.sequence[..25],
                &seq.as_bytes()[start..start + 25],
                "stored sequence must be the forward window at `position`"
            );
            assert_eq!(t.enzyme, EnzymeType::FracMinHash);
        }
    }

    #[test]
    fn contig_offsets_shift_positions_only() {
        let seq = dna(20_000, 17);
        let fmh = FracMinHash::new(21, 50).unwrap();
        let a = fmh.landmarks(seq.as_bytes(), 1, 0);
        let b = fmh.landmarks(seq.as_bytes(), 2, 7_777);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.sequence, y.sequence);
            assert_eq!(y.position, x.position + 7_777);
        }
    }

    #[test]
    fn case_is_normalised() {
        // Soft-masked assemblies must not compare unequal to unmasked ones.
        let seq = dna(20_000, 19);
        let lower = seq.to_lowercase();
        let fmh = FracMinHash::new(21, 50).unwrap();
        let a = fmh.landmarks(seq.as_bytes(), 0, 0);
        let b = fmh.landmarks(lower.as_bytes(), 0, 0);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.sequence, y.sequence, "case must not change the stored tag");
            assert_eq!(x.position, y.position);
        }
    }
}
