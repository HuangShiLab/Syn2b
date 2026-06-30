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

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

/// A 32bp 2bRAD tag with position, enzyme, and strand metadata
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Tag {
    pub sequence: [u8; 32],
    pub position: u64,
    pub enzyme: EnzymeType,
    pub strand: Strand,
}

impl Tag {
    /// Create a new Tag
    pub fn new(sequence: [u8; 32], position: u64, enzyme: EnzymeType, strand: Strand) -> Self {
        Self {
            sequence,
            position,
            enzyme,
            strand,
        }
    }

    /// Return sequence as ASCII string
    pub fn sequence_str(&self) -> String {
        self.sequence
            .iter()
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
            "{}@{}{}:{}",
            seq_str, self.position, strand_char, self.enzyme as u8
        )
    }
}
