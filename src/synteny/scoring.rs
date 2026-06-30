//! Synteny scoring functions
//!
//! Provides metrics for evaluating synteny conservation between genomes,
//! including Jaccard similarity on tag adjacencies, Kendall's tau rank
//! correlation on tag order, breakpoint counting, and path-level scoring.

use crate::synteny::graph::TagAdjacencyGraph;
use crate::tgt::record::TgtRecord;
use crate::tgt::tag::Tag;
use std::collections::{HashMap, HashSet};

// ─────────────────────────────────────────────────────────────────────────────
// Synteny score for a linear path
// ─────────────────────────────────────────────────────────────────────────────

/// Score synteny conservation for a linear path.
///
/// The score is the average edge weight along the path divided by the
/// maximum possible weight (the number of genomes in the graph), adjusted
/// by a length bonus that rewards longer conserved segments:
///
///   score = (mean_weight / num_genomes) * min(1.0, sqrt(path_len / 10.0))
///
/// Higher score = more conserved synteny.
/// - A score of 1.0 means all genomes agree on every adjacency in the path
///   and the path is sufficiently long.
/// - A score of 0.0 means no synteny conservation.
pub fn synteny_score(path: &[u64], graph: &TagAdjacencyGraph) -> f64 {
    if path.len() < 2 {
        return 0.0;
    }

    let mut total_weight = 0u32;
    let mut edge_count = 0usize;

    // Sum edge weights along the path
    for window in path.windows(2) {
        let src = window[0];
        let tgt = window[1];
        // Try both edge directions (undirected path)
        let weight = graph
            .edges
            .get(&(src, tgt))
            .or_else(|| graph.edges.get(&(tgt, src)))
            .map_or(0, |e| e.weight);
        total_weight += weight;
        edge_count += 1;
    }

    if edge_count == 0 || graph.num_genomes == 0 {
        return 0.0;
    }

    let mean_weight = total_weight as f64 / edge_count as f64;
    let max_weight = graph.num_genomes as f64;
    let weight_ratio = mean_weight / max_weight;

    // Length bonus: longer paths score slightly higher, saturating at len=10
    let path_len = path.len() as f64;
    let length_factor = (path_len / 10.0).sqrt().min(1.0);

    weight_ratio * length_factor
}

// ─────────────────────────────────────────────────────────────────────────────
// Pairwise synteny matrix
// ─────────────────────────────────────────────────────────────────────────────

/// Compute pairwise synteny matrix for all genomes in the graph.
///
/// For each pair of genomes (g_i, g_j), computes the Jaccard similarity of
/// their adjacency sets: the set of (tag_a, tag_b) pairs that are adjacent
/// in each genome.
///
/// Returns a map from (genome_i, genome_j) to synteny score in [0, 1].
pub fn pairwise_synteny_matrix(
    graph: &TagAdjacencyGraph,
) -> HashMap<(String, String), f64> {
    let mut matrix = HashMap::new();

    if graph.num_genomes < 2 {
        return matrix;
    }

    // Build adjacency sets for each genome
    let genome_ids: Vec<String> = graph
        .nodes
        .values()
        .flat_map(|n| n.positions.keys().cloned())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    for i in 0..genome_ids.len() {
        for j in (i + 1)..genome_ids.len() {
            let gi = &genome_ids[i];
            let gj = &genome_ids[j];

            let score = pairwise_score(graph, gi, gj);
            matrix.insert((gi.clone(), gj.clone()), score);
            matrix.insert((gj.clone(), gi.clone()), score);
        }
    }

    matrix
}

/// Compute the pairwise synteny score between two genomes using their
/// shared adjacencies in the graph.
fn pairwise_score(graph: &TagAdjacencyGraph, genome_a: &str, genome_b: &str) -> f64 {
    // Collect tag adjacencies for genome A
    let mut adj_a: HashSet<(u64, u64)> = HashSet::new();
    let mut adj_b: HashSet<(u64, u64)> = HashSet::new();

    // Build adjacency sets by traversing the genome order
    if let Some(order_a) = genome_order_in_graph(graph, genome_a) {
        for w in order_a.windows(2) {
            let (u, v) = (w[0], w[1]);
            // Normalize: smaller id first
            adj_a.insert(if u < v { (u, v) } else { (v, u) });
        }
    }

    if let Some(order_b) = genome_order_in_graph(graph, genome_b) {
        for w in order_b.windows(2) {
            let (u, v) = (w[0], w[1]);
            adj_b.insert(if u < v { (u, v) } else { (v, u) });
        }
    }

    if adj_a.is_empty() && adj_b.is_empty() {
        return 0.0;
    }

    let intersection: HashSet<_> = adj_a.intersection(&adj_b).collect();
    let union: HashSet<_> = adj_a.union(&adj_b).collect();

    intersection.len() as f64 / union.len() as f64
}

/// Reconstruct the tag ID ordering for a genome from the graph's node data.
fn genome_order_in_graph(graph: &TagAdjacencyGraph, genome_id: &str) -> Option<Vec<u64>> {
    // Collect all (tag_id, position) pairs for this genome
    let mut tagged: Vec<(u64, u64)> = graph
        .nodes
        .values()
        .filter_map(|node| {
            node.positions
                .get(genome_id)
                .map(|&(pos, _)| (node.tag_id, pos))
        })
        .collect();

    if tagged.is_empty() {
        return None;
    }

    // Sort by genomic position
    tagged.sort_by_key(|&(_, pos)| pos);
    Some(tagged.into_iter().map(|(id, _)| id).collect())
}

// ─────────────────────────────────────────────────────────────────────────────
// Adjacency Jaccard
// ─────────────────────────────────────────────────────────────────────────────

/// Jaccard similarity of adjacent tag pairs between two TGT records.
///
/// For each record, builds a set of (tag_sequence_i, tag_sequence_{i+1}) pairs.
/// The Jaccard similarity is |A ∩ B| / |A ∪ B|.
///
/// This operates on tag *sequences* (not IDs), so it can compare records
/// that have not been inserted into a graph.
///
/// Returns a value in [0, 1] where 1.0 means identical adjacencies.
pub fn adjacency_jaccard(record_a: &TgtRecord, record_b: &TgtRecord) -> f64 {
    let set_a = adjacency_set(record_a);
    let set_b = adjacency_set(record_b);

    if set_a.is_empty() && set_b.is_empty() {
        return 0.0;
    }

    let intersection: HashSet<_> = set_a.intersection(&set_b).collect();
    let union: HashSet<_> = set_a.union(&set_b).collect();

    intersection.len() as f64 / union.len() as f64
}

/// Build a set of (seq_a, seq_b) pairs from adjacent tags in a record.
fn adjacency_set(record: &TgtRecord) -> HashSet<([u8; 32], [u8; 32])> {
    let mut set = HashSet::new();
    for w in record.tags.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        // Normalize order for undirected comparison
        if a.sequence <= b.sequence {
            set.insert((a.sequence, b.sequence));
        } else {
            set.insert((b.sequence, a.sequence));
        }
    }
    set
}

// ─────────────────────────────────────────────────────────────────────────────
// Kendall's tau
// ─────────────────────────────────────────────────────────────────────────────

/// Kendall's tau correlation on tag presence order between two records.
///
/// 1. Find common tags between the two records (by exact sequence match).
/// 2. For each common tag, record its rank (position index) in each genome.
/// 3. Compute Kendall's tau_b: (concordant_pairs - discordant_pairs) / total_pairs.
///
/// Returns a value in [-1, 1] where 1.0 means identical relative order,
/// -1.0 means completely reversed order, and 0.0 means uncorrelated.
pub fn kendall_tag_order(record_a: &TgtRecord, record_b: &TgtRecord) -> f64 {
    // Build rank maps: sequence → position index
    let rank_a = build_rank_map(&record_a.tags);
    let rank_b = build_rank_map(&record_b.tags);

    // Find common tags
    let common: Vec<&[u8; 32]> = rank_a
        .keys()
        .filter(|seq| rank_b.contains_key(*seq))
        .collect();

    let n = common.len();
    if n < 2 {
        return 0.0;
    }

    // Get ranks in both records
    let ranks_a: Vec<usize> = common.iter().map(|seq| rank_a[*seq]).collect();
    let ranks_b: Vec<usize> = common.iter().map(|seq| rank_b[*seq]).collect();

    // Count concordant and discordant pairs
    let mut concordant = 0usize;
    let mut discordant = 0usize;

    for i in 0..n {
        for j in (i + 1)..n {
            let order_a = (ranks_a[i] as i64) - (ranks_a[j] as i64);
            let order_b = (ranks_b[i] as i64) - (ranks_b[j] as i64);

            if order_a == 0 || order_b == 0 {
                // Tied rank — skip in simple Kendall's tau
                continue;
            }

            if order_a.signum() == order_b.signum() {
                concordant += 1;
            } else {
                discordant += 1;
            }
        }
    }

    let total_pairs = concordant + discordant;
    if total_pairs == 0 {
        return 0.0;
    }

    (concordant as f64 - discordant as f64) / total_pairs as f64
}

/// Build a map from tag sequence to its position index (rank) in the genome.
fn build_rank_map(tags: &[Tag]) -> HashMap<[u8; 32], usize> {
    tags.iter()
        .enumerate()
        .map(|(idx, tag)| (tag.sequence, idx))
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Breakpoint count
// ─────────────────────────────────────────────────────────────────────────────

/// Breakpoint count: number of positions where tag adjacency changes.
///
/// Two genomes share an adjacency if they have the same pair of consecutive
/// tags (by sequence). A breakpoint is a position in one genome where the
/// next tag differs from what the other genome has at the corresponding
/// location.
///
/// Algorithm:
/// 1. Find common tags between the two records.
/// 2. Extract the subsequence of common tags in each genome (preserving order).
/// 3. Count how many adjacent pairs in genome A are *not* adjacent in genome B.
///
/// This is a symmetric measure: breakpoint_count(A, B) == breakpoint_count(B, A).
pub fn breakpoint_count(record_a: &TgtRecord, record_b: &TgtRecord) -> usize {
    // Get ranks of common tags in each genome
    let rank_a = build_rank_map(&record_a.tags);
    let rank_b = build_rank_map(&record_b.tags);

    let common_seqs: Vec<[u8; 32]> = rank_a
        .keys()
        .filter(|seq| rank_b.contains_key(*seq))
        .copied()
        .collect();

    if common_seqs.len() < 2 {
        return 0;
    }

    // Build adjacency sets of common tags for each genome
    let mut adj_a: HashSet<([u8; 32], [u8; 32])> = HashSet::new();
    let mut adj_b: HashSet<([u8; 32], [u8; 32])> = HashSet::new();

    // For genome A: scan through tags and record adjacencies between common tags
    for w in record_a.tags.windows(2) {
        let (s1, s2) = (w[0].sequence, w[1].sequence);
        if rank_a.contains_key(&s1)
            && rank_a.contains_key(&s2)
            && rank_b.contains_key(&s1)
            && rank_b.contains_key(&s2)
        {
            let key = normalize_pair(s1, s2);
            adj_a.insert(key);
        }
    }

    for w in record_b.tags.windows(2) {
        let (s1, s2) = (w[0].sequence, w[1].sequence);
        if rank_a.contains_key(&s1)
            && rank_a.contains_key(&s2)
            && rank_b.contains_key(&s1)
            && rank_b.contains_key(&s2)
        {
            let key = normalize_pair(s1, s2);
            adj_b.insert(key);
        }
    }

    // Breakpoints = adjacencies present in one but not the other
    let symmetric_diff: HashSet<_> = adj_a.symmetric_difference(&adj_b).collect();
    symmetric_diff.len()
}

/// Normalize a tag pair for undirected comparison.
fn normalize_pair(a: [u8; 32], b: [u8; 32]) -> ([u8; 32], [u8; 32]) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzyme::enzyme::EnzymeType;
    use crate::tgt::record::TgtRecord;
    use crate::tgt::tag::{Strand, Tag};
    use crate::synteny::graph::TagAdjacencyGraph;

    /// Helper: create a tag with index-encoded sequence
    fn make_tag(idx: u8, position: u64, enzyme: EnzymeType, strand: Strand) -> Tag {
        let mut seq = [b'A'; 32];
        seq[0] = idx;
        seq[1] = idx.wrapping_add(1);
        Tag::new(seq, position, enzyme, strand)
    }

    /// Helper: create a record with evenly spaced tags
    fn make_record(genome_id: &str, tag_count: usize) -> TgtRecord {
        let mut record = TgtRecord::new(genome_id, 1_000_000);
        for i in 0..tag_count {
            let tag = make_tag(i as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward);
            record.add_tag(tag);
        }
        record
    }

    // ─────────────────────────────────────────────────────────────────────────
    // synteny_score
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_synteny_score_perfect() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 10);
        let record_b = make_record("genome_a", 10);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let paths = graph.linear_paths();
        assert!(!paths.is_empty(), "Should have at least one path");

        // With 2 identical genomes, edges have weight 2 (out of max 2)
        // Score should be close to 1.0
        let score = synteny_score(&paths[0], &graph);
        assert!(
            score >= 0.9 && score <= 1.0,
            "Two identical genomes should have near-perfect synteny score, got {}",
            score
        );
    }

    #[test]
    fn test_synteny_score_empty_path() {
        let graph = TagAdjacencyGraph::new();
        let score = synteny_score(&[], &graph);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_synteny_score_singleton_path() {
        let graph = TagAdjacencyGraph::new();
        let score = synteny_score(&[42], &graph);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_synteny_score_single_genome() {
        let mut graph = TagAdjacencyGraph::new();
        let record = make_record("genome_a", 10);
        graph.add_genome("genome_a", &record);
        graph.build_edges();
        graph.simplify(1);

        let paths = graph.linear_paths();
        // With only 1 genome, edge weight = 1, max = 1, ratio = 1.0
        let score = synteny_score(&paths[0], &graph);
        assert!(
            score >= 0.9 && score <= 1.0,
            "Single genome should have perfect score, got {}",
            score
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // adjacency_jaccard
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_adjacency_jaccard_identical() {
        let record_a = make_record("a", 5);
        let record_b = make_record("a", 5);
        let jaccard = adjacency_jaccard(&record_a, &record_b);
        assert!(
            (jaccard - 1.0).abs() < 1e-9,
            "Identical records should have Jaccard = 1.0, got {}",
            jaccard
        );
    }

    #[test]
    fn test_adjacency_jaccard_completely_different() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags 0, 1, 2
        for i in 0..3 {
            record_a.add_tag(make_tag(i, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: tags 10, 11, 12 (no overlap with A)
        for i in 10..13 {
            record_b.add_tag(make_tag(i, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        let jaccard = adjacency_jaccard(&record_a, &record_b);
        assert_eq!(jaccard, 0.0, "No common tags means Jaccard = 0");
    }

    #[test]
    fn test_adjacency_jaccard_partial() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags 0, 1, 2, 3 (adjacencies: (0,1), (1,2), (2,3))
        for i in 0..4 {
            record_a.add_tag(make_tag(i, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: tags 0, 1, 2, 4 (shared adjacencies: (0,1), (1,2); unique: (2,4))
        record_b.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(1, 200, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(2, 300, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(4, 400, EnzymeType::BcgI, Strand::Forward));

        // intersection = {(0,1), (1,2)} → 2
        // union = {(0,1), (1,2), (2,3), (2,4)} → 4
        // Jaccard = 2/4 = 0.5
        let jaccard = adjacency_jaccard(&record_a, &record_b);
        assert!(
            (jaccard - 0.5).abs() < 1e-9,
            "Partial overlap should give Jaccard = 0.5, got {}",
            jaccard
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // kendall_tag_order
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_kendall_identical_order() {
        let record_a = make_record("a", 5);
        let record_b = make_record("a", 5);
        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            (tau - 1.0).abs() < 1e-9,
            "Identical order should give tau = 1.0, got {}",
            tau
        );
    }

    #[test]
    fn test_kendall_reversed_order() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags in order 0, 1, 2, 3
        for i in 0..4 {
            record_a.add_tag(make_tag(i, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: same tags in reversed order 3, 2, 1, 0
        for i in (0..4).rev() {
            record_b.add_tag(make_tag(i, ((3 - i) * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            (tau - (-1.0)).abs() < 1e-9,
            "Reversed order should give tau = -1.0, got {}",
            tau
        );
    }

    #[test]
    fn test_kendall_no_common_tags() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        record_a.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(10, 100, EnzymeType::BcgI, Strand::Forward));

        let tau = kendall_tag_order(&record_a, &record_b);
        assert_eq!(tau, 0.0, "No common tags should give tau = 0");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // breakpoint_count
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_breakpoint_count_identical() {
        let record_a = make_record("a", 5);
        let record_b = make_record("a", 5);
        let bp = breakpoint_count(&record_a, &record_b);
        assert_eq!(bp, 0, "Identical genomes should have 0 breakpoints");
    }

    #[test]
    fn test_breakpoint_count_with_rearrangement() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags 0, 1, 2, 3, 4 (adjacencies: 0-1, 1-2, 2-3, 3-4)
        for i in 0..5 {
            record_a.add_tag(make_tag(i, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: tags 0, 1, 3, 2, 4 (breakpoint at 1-3 and 3-2 vs 2-3)
        // Adjacencies in B: 0-1, 1-3, 3-2, 2-4
        // Common adjacencies: 0-1
        // Breakpoints: (1-2 in A but not B), (2-3 in A, 3-2 in B = same undirected),
        //              (3-4 in A but not B)
        record_b.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(1, 200, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(3, 300, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(2, 400, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(4, 500, EnzymeType::BcgI, Strand::Forward));

        let bp = breakpoint_count(&record_a, &record_b);
        assert!(bp > 0, "Rearranged genome should have breakpoints, got {}", bp);
    }

    #[test]
    fn test_breakpoint_count_empty() {
        let record_a = TgtRecord::new("a", 1_000_000);
        let record_b = TgtRecord::new("b", 1_000_000);
        let bp = breakpoint_count(&record_a, &record_b);
        assert_eq!(bp, 0, "Empty records should have 0 breakpoints");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // pairwise_synteny_matrix
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_pairwise_matrix_identical_genomes() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 10);
        let record_b = make_record("genome_a", 10);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let matrix = pairwise_synteny_matrix(&graph);
        let score = matrix
            .get(&("genome_a".to_string(), "genome_b".to_string()))
            .copied()
            .unwrap_or(0.0);
        assert!(
            score >= 0.9,
            "Identical genomes should have high pairwise score, got {}",
            score
        );
    }

    #[test]
    fn test_pairwise_matrix_single_genome() {
        let mut graph = TagAdjacencyGraph::new();
        let record = make_record("genome_a", 5);
        graph.add_genome("genome_a", &record);
        graph.build_edges();

        let matrix = pairwise_synteny_matrix(&graph);
        assert!(matrix.is_empty(), "Single genome should produce empty matrix");
    }
}
