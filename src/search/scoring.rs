//! Placement scoring and conservative confidence assignment.

use serde::{Deserialize, Serialize};

/// Initial placement weights from the Syn2b-Contig design.
pub const W_SHARED_TAGS: f64 = 0.15;
pub const W_ADJACENCIES: f64 = 0.30;
pub const W_CHAIN_COVERAGE: f64 = 0.30;
pub const W_ORIENTATION: f64 = 0.10;
pub const W_UNIQUENESS: f64 = 0.10;
pub const W_INTERIOR: f64 = 0.05;

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct PlacementComponents {
    pub shared_tag_fraction: f64,
    pub conserved_adjacency_fraction: f64,
    pub chained_block_coverage: f64,
    pub orientation_consistency: f64,
    pub uniqueness_penalty: f64,
    pub contig_interior_support: f64,
}

pub fn placement_score(components: &PlacementComponents) -> f64 {
    let raw = W_SHARED_TAGS * components.shared_tag_fraction
        + W_ADJACENCIES * components.conserved_adjacency_fraction
        + W_CHAIN_COVERAGE * components.chained_block_coverage
        + W_ORIENTATION * components.orientation_consistency
        + W_UNIQUENESS * components.uniqueness_penalty
        + W_INTERIOR * components.contig_interior_support;
    raw.clamp(0.0, 1.0)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum PlacementConfidence {
    High,
    Medium,
    Low,
    Exploratory,
}

impl PlacementConfidence {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::High => "high",
            Self::Medium => "medium",
            Self::Low => "low",
            Self::Exploratory => "exploratory",
        }
    }
}

/// Length-dependent support thresholds from the handoff, used only to assign
/// a coarse confidence label. Statistical calibration remains future work.
pub fn support_thresholds(query_length: u64) -> (usize, usize) {
    match query_length {
        500_000.. => (25, 15),
        100_000..=499_999 => (15, 8),
        50_000..=99_999 => (10, 5),
        10_000..=49_999 => (5, 2),
        _ => (0, 0),
    }
}

pub fn confidence(
    score: f64,
    shared_tags: usize,
    shared_adjacencies: usize,
    query_length: u64,
) -> PlacementConfidence {
    let (min_tags, min_adj) = support_thresholds(query_length);
    let enough_support = shared_tags >= min_tags && shared_adjacencies >= min_adj;
    if score >= 0.80 && enough_support {
        PlacementConfidence::High
    } else if score >= 0.65 {
        PlacementConfidence::Medium
    } else if score >= 0.50 && query_length >= 10_000 {
        PlacementConfidence::Low
    } else {
        PlacementConfidence::Exploratory
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_components_score_one() {
        let c = PlacementComponents {
            shared_tag_fraction: 1.0,
            conserved_adjacency_fraction: 1.0,
            chained_block_coverage: 1.0,
            orientation_consistency: 1.0,
            uniqueness_penalty: 1.0,
            contig_interior_support: 1.0,
        };
        assert_eq!(placement_score(&c), 1.0);
    }

    #[test]
    fn thresholds_follow_query_length() {
        assert_eq!(support_thresholds(600_000), (25, 15));
        assert_eq!(support_thresholds(120_000), (15, 8));
        assert_eq!(support_thresholds(20_000), (5, 2));
        assert_eq!(support_thresholds(1_000), (0, 0));
    }
}
