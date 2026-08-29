//! Tag data structure — 32bp 2bRAD tag with metadata

use crate::enzyme::enzyme::EnzymeType;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Strand orientation of a tag
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Convert to u8 (0 = Forward, 1 = Reverse)
    pub fn to_u8(&self) -> u8 {
        match self {
            Strand::Forward => 0,
            Strand::Reverse => 1,
        }
    }

    /// Convert from u8 (0 = Forward, anything else = Reverse)
    pub fn from_u8(v: u8) -> Self {
        if v == 0 { Strand::Forward } else { Strand::Reverse }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// A 32bp 2bRAD tag with position, enzyme, strand, and contig metadata
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Tag {
    pub sequence: [u8; 32],
    pub position: u64,
    pub enzyme: EnzymeType,
    pub strand: Strand,
    pub contig_id: u16,   // 0 = not specified / single contig; 1+ = index into contig_names
}

impl Tag {
    /// Create a new Tag
    pub fn new(sequence: [u8; 32], position: u64, enzyme: EnzymeType, strand: Strand, contig_id: u16) -> Self {
        Self {
            sequence,
            position,
            enzyme,
            strand,
            contig_id,
        }
    }

    /// Return sequence as ASCII string (trims trailing null bytes from fixed 32-byte buffer)
    pub fn sequence_str(&self) -> String {
        self.sequence
            .iter()
            .take_while(|&&b| b != 0)
            .map(|&b| b as char)
            .collect()
    }

    /// Number of real bases in `sequence`, which is a zero-padded fixed buffer.
    pub fn seq_len(&self) -> usize {
        self.sequence.iter().take_while(|&&b| b != 0).count()
    }

    /// Strand-canonical sequence: `min(sequence, reverse_complement(sequence))`,
    /// zero-padded back to 32 bytes.
    ///
    /// # Why this is required for structural comparison
    ///
    /// A tag inside an inverted segment is read from the opposite strand, so the
    /// same locus yields the reverse complement of the sequence seen in the
    /// un-inverted genome. Keying tag identity on the raw forward window
    /// therefore makes every tag in an inverted region *vanish* from the shared
    /// set. That looks like a strong inversion signal, but it is the same
    /// mechanism by which substitutions destroy tags, so the two causes become
    /// inseparable: measured on E. coli K-12, a genome with a 400 kb inversion
    /// and zero substitutions scored 0.8234 while a genome with 0.1%
    /// substitutions and no inversion scored 0.8678 — indistinguishable.
    ///
    /// Canonicalising makes an inverted tag match its homolog, which removes the
    /// substitution confound. The inversion is then visible in the tag *order*
    /// instead, which is what [`crate::synteny::scoring::structural_synteny`]
    /// measures.
    pub fn canonical_sequence(&self) -> [u8; 32] {
        let n = self.seq_len();
        let rc = self.revcomp_sequence();
        if rc[..n] < self.sequence[..n] {
            rc
        } else {
            self.sequence
        }
    }

    /// Reverse complement of the stored sequence, in the same zero-padded buffer.
    pub fn revcomp_sequence(&self) -> [u8; 32] {
        let n = self.seq_len();
        let mut rc = [0u8; 32];
        for i in 0..n {
            rc[n - 1 - i] = match self.sequence[i] {
                b'A' | b'a' => b'T',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'T' | b't' => b'A',
                other => other,
            };
        }
        rc
    }

    /// True when the stored window is the reverse complement of its canonical
    /// representative — the bit that flips when this locus is inverted.
    ///
    /// The digester always stores the window as read off the forward strand of
    /// the assembly, so if a locus is inverted between two genomes the two
    /// stored sequences are reverse complements, they share a canonical form
    /// (and so still match), and exactly one of them reports `true` here. That
    /// makes this an inversion indicator that survives canonicalisation.
    pub fn is_revcomp_of_canonical(&self) -> bool {
        let n = self.seq_len();
        self.revcomp_sequence()[..n] < self.sequence[..n]
    }

    /// True when the tag is its own reverse complement, so
    /// [`Self::is_revcomp_of_canonical`] is `false` in every genome and the
    /// landmark cannot report orientation at all.
    pub fn is_palindromic(&self) -> bool {
        let n = self.seq_len();
        n > 0 && self.revcomp_sequence()[..n] == self.sequence[..n]
    }

    /// Hamming distance between two tag sequences
    pub fn hamming_distance(&self, other: &Tag) -> u8 {
        self.sequence
            .iter()
            .zip(other.sequence.iter())
            .map(|(&a, &b)| if a != b { 1 } else { 0 })
            .sum()
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let seq_str = self.sequence_str();
        let strand_char = match self.strand {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        };
        write!(
            f,
            "{}@{}{}:{}:{}",
            seq_str, self.position, strand_char, self.enzyme as u8, self.contig_id
        )
    }
}
