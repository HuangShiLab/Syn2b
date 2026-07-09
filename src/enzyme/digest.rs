//! In silico digestion engine.
//!
//! Provides functions to simulate restriction enzyme digestion of genomic
//! sequences and extract 2bRAD tags.
//!
//! Algorithm (based on Fast2bRAD-M / 2bRAD-M paper):
//! For each enzyme, slide a window of `tag_length` across the sequence.
//! If all anchors in any pattern match within the window, and the window
//! contains only A/T/C/G bases (no degenerate bases), the window is a valid tag.

use crate::enzyme::enzyme::{Enzyme, EnzymeType};
use crate::tgt::tag::{Strand, Tag};

// ── Lookup table: only A/T/C/G/a/t/c/g are true ────────────────────────────

const ATCG_TABLE: [bool; 256] = {
    let mut table = [false; 256];
    table[b'A' as usize] = true; table[b'a' as usize] = true;
    table[b'T' as usize] = true; table[b't' as usize] = true;
    table[b'C' as usize] = true; table[b'c' as usize] = true;
    table[b'G' as usize] = true; table[b'g' as usize] = true;
    table
};

#[inline]
/// Check if all bytes in the window are A/T/C/G (no degenerate bases).
pub fn is_pure_atcg(window: &[u8]) -> bool {
    window.iter().all(|&b| ATCG_TABLE[b as usize])
}

// ── Internal helper ─────────────────────────────────────────────────────────

fn digest_sequence(sequence: &[u8], enzyme: EnzymeType, contig_id: u16, offset: u64) -> Vec<Tag> {
    let props = Enzyme::properties(enzyme);
    let tag_len = props.tag_length as usize;
    let mut tags = Vec::new();

    for pattern in props.patterns {
        // Find the first anchor (smallest offset) to use as skip-anchor
        let Some(first_anchor) = pattern.anchors.iter().min_by_key(|a| a.offset) else {
            // Fallback for patterns with no anchors (shouldn't happen in practice)
            for pos in 0..sequence.len().saturating_sub(tag_len - 1) {
                if pos + tag_len > sequence.len() {
                    break;
                }
                let window = &sequence[pos..pos + tag_len];
                if pattern.matches(window) && is_pure_atcg(window) {
                    let mut tag_seq = [0u8; 32];
                    let copy_len = tag_len.min(32);
                    tag_seq[..copy_len].copy_from_slice(&window[..copy_len]);
                    tags.push(Tag::new(
                        tag_seq,
                        offset + pos as u64,
                        enzyme,
                        Strand::Forward,
                        contig_id,
                    ));
                }
            }
            continue;
        };

        let first_byte = first_anchor.motif[0];
        for candidate_pos in memchr::memchr_iter(first_byte, sequence) {
            if candidate_pos < first_anchor.offset {
                continue;
            }
            let window_start = candidate_pos - first_anchor.offset;
            if window_start + tag_len > sequence.len() {
                continue;
            }
            let window = &sequence[window_start..window_start + tag_len];
            if pattern.matches(window) && is_pure_atcg(window) {
                let mut tag_seq = [0u8; 32];
                let copy_len = tag_len.min(32);
                tag_seq[..copy_len].copy_from_slice(&window[..copy_len]);
                tags.push(Tag::new(
                    tag_seq,
                    offset + window_start as u64,
                    enzyme,
                    Strand::Forward,
                    contig_id,
                ));
            }
        }
    }

    tags.sort_by_key(|t| t.position);
    tags.dedup_by_key(|t| t.position);
    tags
}

// ── Public API ────────────────────────────────────────────────────────────

/// Digest a single contig sequence with cumulative offset and contig_id.
/// 
/// `offset` is the cumulative length of all preceding contigs (added to position).
/// `contig_id` is a 1-based identifier for this contig (0 = not specified).
pub fn digest_genome_contig(sequence: &[u8], enzyme: EnzymeType, contig_id: u16, offset: u64) -> Vec<Tag> {
    digest_sequence(sequence, enzyme, contig_id, offset)
}

/// Digest a genome sequence with a single enzyme and extract tags.
///
/// For each pattern defined by the enzyme, slides a `tag_length` window
/// across the sequence. If every anchor in the pattern matches, and the
/// window contains only A/T/C/G bases, the window is emitted as a `Tag`.
/// Tags are sorted by position and deduplicated.
pub fn digest_genome(sequence: &[u8], enzyme: EnzymeType) -> Vec<Tag> {
    digest_sequence(sequence, enzyme, 0, 0)
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_digest_empty_sequence() {
        let tags = digest_genome(b"", EnzymeType::BcgI);
        assert!(tags.is_empty());
    }

    #[test]
    fn test_digest_no_sites() {
        let tags = digest_genome(b"ATATATATATATATATATATATATATATAT", EnzymeType::BcgI);
        assert!(tags.is_empty());
    }

    #[test]
    fn test_digest_bcgI_real_site() {
        // BcgI: tag_length=32, forward anchors at offset 10=CGA, offset 19=TGC
        // Build a 32-bp window: 10 A's + CGA + 6 A's + TGC + 10 A's
        let seq = b"AAAAAAAAAACGAAAAAAATGCAAAAAAAAAA";
        assert_eq!(seq.len(), 32);
        let tags = digest_genome(seq, EnzymeType::BcgI);
        assert_eq!(tags.len(), 1);
        assert_eq!(tags[0].position, 0);
        // Verify reverse pattern matches reverse-complement sequence
        let rev_seq = b"TTTTTTTTTTGCATTTTTTTCGTTTTTTTTTT";
        let tags_rev = digest_genome(rev_seq, EnzymeType::BcgI);
        assert_eq!(tags_rev.len(), 1);
    }

    #[test]
    fn test_digest_bcgI_with_n_rejected() {
        // Same as above but with an N — should be rejected by is_pure_atcg
        let seq = b"AAAAAAAAAACGAAAAAAATGCAAAANAAA";
        let tags = digest_genome(seq, EnzymeType::BcgI);
        assert!(tags.is_empty(), "N-containing window should be rejected");
    }

    #[test]
    fn test_digest_cjeI_correct_site() {
        // CjeI: tag_length=37, forward anchors at offset 8=CCA, offset 17=GT
        // Build a 37-bp window: 8 A's + CCA + 6 A's + GT + 18 A's
        let seq = b"AAAAAAAACCAAAAAAAGTAAAAAAAAAAAAAAAAAA";
        assert_eq!(seq.len(), 37);
        let tags = digest_genome(seq, EnzymeType::CjeI);
        assert_eq!(tags.len(), 1);
    }

    #[test]
    fn test_digest_cjePI_correct_site() {
        // CjePI: tag_length=38, forward anchors at offset 7=CCA, offset 17=TC
        // Build a 38-bp window: 7 A's + CCA + 7 A's + TC + 19 A's
        let seq = b"AAAAAAACCAAAAAAAATCAAAAAAAAAAAAAAAAAAA";
        assert_eq!(seq.len(), 38);
        let tags = digest_genome(seq, EnzymeType::CjePI);
        assert_eq!(tags.len(), 1);
    }
}
