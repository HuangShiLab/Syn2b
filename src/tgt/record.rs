//! TGT record — a genome's complete tag-gap-tag representation

use crate::tgt::gap::Gap;
use crate::tgt::tag::Tag;
use serde::{Deserialize, Serialize};
use std::fmt;

/// A single genome's TGT record containing ordered tags and gaps
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TgtRecord {
    pub genome_id: String,
    pub tags: Vec<Tag>,
    pub gaps: Vec<Gap>,
    pub total_length: u64,
}

impl TgtRecord {
    /// Create a new empty TGT record
    pub fn new(genome_id: &str, total_length: u64) -> Self {
        Self {
            genome_id: genome_id.to_string(),
            tags: Vec::new(),
            gaps: Vec::new(),
            total_length,
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
        }
        self.tags.push(tag);
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

    /// Number of distinct enzymes represented among this record's tags.
    pub fn enzyme_count(&self) -> usize {
        use std::collections::HashSet;
        self.tags
            .iter()
            .map(|t| t.enzyme)
            .collect::<HashSet<_>>()
            .len()
    }
}

impl fmt::Display for TgtRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, ">{}|length={}", self.genome_id, self.total_length)?;
        if self.tags.is_empty() {
            return Ok(());
        }
        // Group tags by enzyme
        use std::collections::HashMap;
        let mut by_enzyme: HashMap<String, Vec<(usize, &Tag)>> = HashMap::new();
        for (i, tag) in self.tags.iter().enumerate() {
            let key = format!("{:?}", tag.enzyme);
            by_enzyme.entry(key).or_default().push((i, tag));
        }
        for (enzyme, tag_refs) in by_enzyme {
            write!(f, "{}:", enzyme)?;
            for (j, (orig_idx, tag)) in tag_refs.iter().enumerate() {
                if j > 0 && *orig_idx > 0 {
                    // Find gap between this and previous tag in full record
                    if *orig_idx <= self.gaps.len() {
                        write!(f, " {}", self.gaps[*orig_idx - 1])?;
                    }
                }
                write!(f, " {}", tag.sequence_str())?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}
