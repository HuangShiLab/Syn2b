//! TGT record — a genome's complete tag-gap-tag representation

use crate::tgt::gap::Gap;
use crate::tgt::tag::Tag;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fmt;

/// A single genome's TGT record containing ordered tags and gaps
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TgtRecord {
    pub genome_id: String,
    pub tags: Vec<Tag>,
    pub gaps: Vec<Gap>,
    pub total_length: u64,
    #[serde(skip)]
    pub sequence_set: HashSet<[u8; 32]>,
    #[serde(skip)]
    pub adjacency_set: HashSet<([u8; 32], [u8; 32])>,
    pub contig_names: Vec<String>,   // contig_id 1+ index here; contig_id 0 means single contig / not specified
    pub contig_offsets: Vec<u64>,  // cumulative length of contigs up to index i
}

impl TgtRecord {
    /// Create a new empty TGT record
    pub fn new(genome_id: &str, total_length: u64) -> Self {
        Self {
            genome_id: genome_id.to_string(),
            tags: Vec::new(),
            gaps: Vec::new(),
            total_length,
            sequence_set: HashSet::new(),
            adjacency_set: HashSet::new(),
            contig_names: Vec::new(),
            contig_offsets: Vec::new(),
        }
    }

    /// Append a tag, auto-computing the gap from the previous tag
    pub fn add_tag(&mut self, tag: Tag) {
        if let Some(last_tag) = self.tags.last() {
            let gap_size = if tag.position > last_tag.position {
                (tag.position - last_tag.position) as u32
            } else {
                0
            };
            self.gaps.push(Gap::new(gap_size));
            // Update adjacency set for fast pairwise scoring
            let (a, b) = if last_tag.sequence <= tag.sequence {
                (last_tag.sequence, tag.sequence)
            } else {
                (tag.sequence, last_tag.sequence)
            };
            self.adjacency_set.insert((a, b));
        }
        self.sequence_set.insert(tag.sequence);
        self.tags.push(tag);
    }

    /// Return the number of distinct enzymes used in this record
    pub fn enzyme_count(&self) -> usize {
        use std::collections::HashSet;
        let enzymes: HashSet<_> = self.tags.iter().map(|t| t.enzyme).collect();
        enzymes.len()
    }

    /// Return the number of tags
    pub fn tag_count(&self) -> usize {
        self.tags.len()
    }

    /// Return the mean gap size
    pub fn mean_gap(&self) -> f64 {
        if self.gaps.is_empty() {
            return 0.0;
        }
        let sum: u64 = self.gaps.iter().map(|g| g.size as u64).sum();
        sum as f64 / self.gaps.len() as f64
    }

    /// Return the median gap size
    pub fn median_gap(&self) -> f64 {
        if self.gaps.is_empty() {
            return 0.0;
        }
        let mut sizes: Vec<u32> = self.gaps.iter().map(|g| g.size).collect();
        sizes.sort();
        let mid = sizes.len() / 2;
        if sizes.len() % 2 == 0 {
            (sizes[mid - 1] as f64 + sizes[mid] as f64) / 2.0
        } else {
            sizes[mid] as f64
        }
    }

    /// Return the maximum gap size
    pub fn max_gap(&self) -> u32 {
        self.gaps.iter().map(|g| g.size).max().unwrap_or(0)
    }

    /// Estimated genome coverage fraction
    pub fn coverage_fraction(&self) -> f64 {
        if self.total_length == 0 {
            return 0.0;
        }
        let tag_bp = self.tags.len() * 32;
        tag_bp as f64 / self.total_length as f64
    }
}

impl fmt::Display for TgtRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, ">{}|length={}", self.genome_id, self.total_length)?;
        // Write contig metadata if available
        if !self.contig_names.is_empty() && !self.contig_offsets.is_empty() {
            let mut contig_parts = Vec::new();
            for i in 0..self.contig_names.len() {
                let start = self.contig_offsets[i];
                let end = if i + 1 < self.contig_offsets.len() {
                    self.contig_offsets[i + 1]
                } else {
                    self.total_length
                };
                let len = end - start;
                contig_parts.push(format!("{}:{}", self.contig_names[i], len));
            }
            writeln!(f, "#contigs={}", contig_parts.join(";"))?;
        }
        if self.tags.is_empty() {
            return Ok(());
        }
        // Write all tags in genome order with position and gaps
        for (i, tag) in self.tags.iter().enumerate() {
            if i > 0 && i <= self.gaps.len() {
                write!(f, " {}", self.gaps[i - 1])?;
            }
            if i > 0 {
                write!(f, " ")?;
            }
            if tag.contig_id > 0 && (tag.contig_id as usize - 1) < self.contig_names.len() {
                write!(f, "{}:{}@{}:{}", tag.enzyme, tag.sequence_str(), tag.position, self.contig_names[tag.contig_id as usize - 1])?;
            } else if tag.contig_id > 0 {
                write!(f, "{}:{}@{}:{}", tag.enzyme, tag.sequence_str(), tag.position, tag.contig_id)?;
            } else {
                write!(f, "{}:{}@{}", tag.enzyme, tag.sequence_str(), tag.position)?;
            }
        }
        writeln!(f)?;
        Ok(())
    }
}
