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
