//! Synteny detection engine — tag adjacency graph, scoring, and block extraction
//!
//! This module implements the core synteny detection pipeline:
//!
//! 1. **Graph construction** (`graph.rs`): Build a tag adjacency graph from
//!    multiple genomes' TGT records. Shared tags become nodes; consecutive tags
//!    in a genome become edges.
//!
//! 2. **Graph simplification** (`graph.rs`): Remove low-support edges and
//!    isolated nodes using ntSynt-inspired anchored node classification.
//!
//! 3. **Backbone extraction** (`graph.rs`): Extract linear paths (degree-2
//!    chains) from the simplified graph. These are the synteny backbones.
//!
//! 4. **Scoring** (`scoring.rs`): Evaluate synteny conservation using edge
//!    weights, Jaccard similarity, Kendall's tau, and breakpoint counting.
//!
//! 5. **Block extraction** (`blocks.rs`): Convert linear paths into synteny
//!    blocks with genomic coordinates, detect indels, and filter by size.

pub mod blocks;
pub mod graph;
pub mod scoring;

pub use blocks::{
    extract_synteny_blocks, detect_indels, filter_blocks_by_size, SyntenyBlock,
};
pub use graph::{AdjacencyEdge, TagAdjacencyGraph, TagNode};
pub use scoring::{
    adjacency_jaccard, breakpoint_count, kendall_tag_order, pairwise_synteny_matrix, position_correlation, windowed_position_correlation, windowed_density_correlation, hamming_distance, count_approximate_common_tags, windowed_synteny_score, windowed_synteny_score_approx,
    synteny_score,
};
