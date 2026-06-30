//! In silico digestion engine.
//!
//! Provides functions to simulate restriction enzyme digestion of genomic
//! sequences and extract 2bRAD tags.

use crate::enzyme::enzyme::{Enzyme, EnzymeType};
use crate::tgt::tag::{Strand, Tag};

/// Find all recognition sites of an enzyme in a DNA sequence.
fn find_sites(enzyme: &Enzyme, sequence: &[u8]) -> Vec<usize> {
    let site = enzyme.recognition_site;
    let mut positions = Vec::new();
    if site.is_empty() || sequence.len() < site.len() {
        return positions;
    }
    for i in 0..=sequence.len() - site.len() {
        if matches_site(sequence, i, site) {
            positions.push(i);
        }
    }
    positions
}

/// Check if the sequence at `offset` matches the recognition site.
/// Handles degenerate bases in the site (R=A|G, Y=C|T, N=any).
fn matches_site(seq: &[u8], offset: usize, site: &[u8]) -> bool {
    for (i, &base) in site.iter().enumerate() {
        let seq_base = seq[offset + i];
        let match_ok = match base {
            b'A' | b'T' | b'G' | b'C' => seq_base == base,
            b'N' => true,
            b'R' => seq_base == b'A' || seq_base == b'G',
            b'Y' => seq_base == b'C' || seq_base == b'T',
            b'M' => seq_base == b'A' || seq_base == b'C',
            b'K' => seq_base == b'G' || seq_base == b'T',
            b'S' => seq_base == b'G' || seq_base == b'C',
            b'W' => seq_base == b'A' || seq_base == b'T',
            b'H' => seq_base == b'A' || seq_base == b'C' || seq_base == b'T',
            b'B' => seq_base == b'C' || seq_base == b'G' || seq_base == b'T',
            b'V' => seq_base == b'A' || seq_base == b'C' || seq_base == b'G',
            b'D' => seq_base == b'A' || seq_base == b'G' || seq_base == b'T',
            _ => false,
        };
        if !match_ok {
            return false;
        }
    }
    true
}

/// Digest a genome sequence with a single enzyme and extract tags.
///
/// # Arguments
/// * `sequence` — The genomic DNA sequence (uppercase A/C/G/T bytes)
/// * `enzyme` — The enzyme type to use for digestion
///
/// # Returns
/// A vector of `Tag` structs, each representing an extracted 2bRAD tag.
pub fn digest_genome(sequence: &[u8], enzyme: EnzymeType) -> Vec<Tag> {
    let props = Enzyme::properties(enzyme);
    let sites = find_sites(&props, sequence);
    let tag_len = props.tag_length as usize;
    let mut tags = Vec::new();

    for &site_pos in &sites {
        // Compute the 5' and 3' cut positions
        let cut_5 = (site_pos as i64 + props.cut_offset_5 as i64).max(0) as usize;
        let cut_3 = (site_pos as i64 + props.cut_offset_3 as i64).max(0) as usize;

        if cut_3 <= cut_5 || cut_3 - cut_5 != tag_len {
            continue;
        }
        if cut_3 > sequence.len() {
            continue;
        }

        // Extract the tag sequence between the two cut sites
        let mut tag_seq = [0u8; 32];
        let extracted_len = cut_3 - cut_5;
        if extracted_len <= 32 {
            tag_seq[..extracted_len].copy_from_slice(&sequence[cut_5..cut_3]);
        } else {
            // Shouldn't happen if tag_length matches, but clamp to 32
            tag_seq.copy_from_slice(&sequence[cut_5..cut_5 + 32]);
        }

        tags.push(Tag::new(
            tag_seq,
            cut_5 as u64,
            enzyme,
            Strand::Forward,
        ));
    }

    tags.sort_by_key(|t| t.position);
    tags
}

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
}
