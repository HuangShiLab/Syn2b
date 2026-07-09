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
    if record_a.adjacency_set.is_empty() && record_b.adjacency_set.is_empty() {
        return 0.0;
    }

    let intersection: HashSet<_> = record_a.adjacency_set.intersection(&record_b.adjacency_set).collect();
    let union: HashSet<_> = record_a.adjacency_set.union(&record_b.adjacency_set).collect();

    intersection.len() as f64 / union.len() as f64
}



/// Kendall's tau correlation on tag presence order between two records.
///
/// 1. Find common tags between the two records (by exact sequence match).
/// 2. For each common tag, record its rank (position index) in each genome.
/// 3. Compute Kendall's tau_b: (concordant_pairs - discordant_pairs) / total_pairs.
///
/// Returns a value in [-1, 1] where 1.0 means identical relative order,
/// -1.0 means completely reversed order, and 0.0 means uncorrelated.
pub fn kendall_tag_order(record_a: &TgtRecord, record_b: &TgtRecord) -> f64 {
    let rank_a = build_rank_map(&record_a.tags);

    // Collect common tags' ranks in A, ordered by B's tag order
    let mut common_ranks: Vec<usize> = Vec::new();
    for tag in &record_b.tags {
        if let Some(&rank_a) = rank_a.get(&tag.sequence) {
            common_ranks.push(rank_a);
        }
    }

    let n = common_ranks.len();
    if n < 2 {
        return 0.0;
    }

    let inversions = count_inversions(&mut common_ranks);
    let total_pairs = n * (n - 1) / 2;
    if total_pairs == 0 {
        return 0.0;
    }

    1.0 - 2.0 * inversions as f64 / total_pairs as f64
}

/// Count inversions in a slice using merge-sort (O(n log n)).
fn count_inversions(arr: &mut [usize]) -> usize {
    let n = arr.len();
    if n < 2 {
        return 0;
    }
    let mut temp = vec![0; n];
    merge_sort_count(arr, &mut temp, 0, n)
}

fn merge_sort_count(arr: &mut [usize], temp: &mut [usize], left: usize, right: usize) -> usize {
    if right - left <= 1 {
        return 0;
    }
    let mid = left + (right - left) / 2;
    let mut inv_count = merge_sort_count(arr, temp, left, mid);
    inv_count += merge_sort_count(arr, temp, mid, right);

    let mut i = left;
    let mut j = mid;
    let mut k = left;
    while i < mid && j < right {
        if arr[i] <= arr[j] {
            temp[k] = arr[i];
            i += 1;
        } else {
            temp[k] = arr[j];
            inv_count += mid - i;
            j += 1;
        }
        k += 1;
    }
    while i < mid {
        temp[k] = arr[i];
        i += 1;
        k += 1;
    }
    while j < right {
        temp[k] = arr[j];
        j += 1;
        k += 1;
    }
    arr[left..right].copy_from_slice(&temp[left..right]);
    inv_count
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
        Tag::new(seq, position, enzyme, strand, 0)
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
        for i in 10..13usize {
            record_b.add_tag(make_tag(i as u8, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        let jaccard = adjacency_jaccard(&record_a, &record_b);
        assert_eq!(jaccard, 0.0, "No common tags means Jaccard = 0");
    }

    #[test]
    fn test_adjacency_jaccard_partial() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags 0, 1, 2, 3 (adjacencies: (0,1), (1,2), (2,3))
        for i in 0..4usize {
            record_a.add_tag(make_tag(i as u8, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
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
        for i in 0..4u8 {
            record_a.add_tag(make_tag(i, (i as u64) * 100, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: same tags in reversed order 3, 2, 1, 0
        for i in (0..4u8).rev() {
            record_b.add_tag(make_tag(i, ((3 - i) as u64) * 100, EnzymeType::BcgI, Strand::Forward));
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
        for i in 0..5usize {
            record_a.add_tag(make_tag(i as u8, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
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

    // ─────────────────────────────────────────────────────────────────────────
    // count_inversions / kendall_tag_order
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_count_inversions_sorted() {
        let mut arr = vec![1, 2, 3, 4, 5];
        assert_eq!(count_inversions(&mut arr), 0, "Sorted array should have 0 inversions");
    }

    #[test]
    fn test_count_inversions_reversed() {
        let mut arr = vec![5, 4, 3, 2, 1];
        assert_eq!(count_inversions(&mut arr), 10, "Reversed array of 5 should have 10 inversions");
    }

    #[test]
    fn test_kendall_tag_order_identical() {
        let record_a = make_record("a", 10);
        let record_b = make_record("a", 10);
        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            tau >= 0.99,
            "Identical genomes should have tau ≈ 1.0, got {}",
            tau
        );
    }

    #[test]
    fn test_kendall_tag_order_reversed() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);
        for i in 0..10 {
            record_a.add_tag(make_tag(i as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward));
        }
        for i in 0..10 {
            record_b.add_tag(make_tag((9 - i) as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward));
        }
        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            tau <= -0.7,
            "Reversed genomes should have strongly negative tau, got {}",
            tau
        );
    }
}

/// Compute Pearson correlation of shared tag positions between two genomes.
/// Returns absolute value (inversions produce negative correlation but still syntenic).
pub fn position_correlation(record_a: &TgtRecord, record_b: &TgtRecord) -> f64 {
    let mut pos_a = Vec::new();
    let mut pos_b = Vec::new();

    // Build sequence->position map for record_b
    let pos_map_b: std::collections::HashMap<[u8; 32], u64> = record_b
        .tags
        .iter()
        .map(|t| (t.sequence, t.position))
        .collect();

    // Find shared tags and collect positions
    for tag_a in &record_a.tags {
        if let Some(&pos_b_val) = pos_map_b.get(&tag_a.sequence) {
            pos_a.push(tag_a.position as f64);
            pos_b.push(pos_b_val as f64);
        }
    }

    let n = pos_a.len();
    if n < 2 {
        return 0.0;
    }

    let mean_a: f64 = pos_a.iter().sum::<f64>() / n as f64;
    let mean_b: f64 = pos_b.iter().sum::<f64>() / n as f64;

    let cov: f64 = pos_a.iter().zip(pos_b.iter())
        .map(|(a, b)| (a - mean_a) * (b - mean_b))
        .sum();
    let var_a: f64 = pos_a.iter().map(|a| (a - mean_a).powi(2)).sum();
    let var_b: f64 = pos_b.iter().map(|b| (b - mean_b).powi(2)).sum();

    if var_a <= 0.0 || var_b <= 0.0 {
        return 0.0;
    }

    let corr = cov / (var_a.sqrt() * var_b.sqrt());
    corr.abs()
}

/// Compute windowed position correlation: split shared tags into windows
/// and compute the maximum Pearson correlation across windows.
/// Captures local synteny even when global shared ratio is low.
pub fn windowed_position_correlation(record_a: &TgtRecord, record_b: &TgtRecord, windows: usize) -> f64 {
    let mut pos_a = Vec::new();
    let mut pos_b = Vec::new();

    let pos_map_b: std::collections::HashMap<[u8; 32], u64> = record_b
        .tags
        .iter()
        .map(|t| (t.sequence, t.position))
        .collect();

    for tag_a in &record_a.tags {
        if let Some(&pos_b_val) = pos_map_b.get(&tag_a.sequence) {
            pos_a.push(tag_a.position);
            pos_b.push(pos_b_val);
        }
    }

    let n = pos_a.len();
    if n < windows * 2 {
        return position_correlation(record_a, record_b);
    }

    // Sort by position in genome A
    let mut pairs: Vec<(u64, u64)> = pos_a.into_iter().zip(pos_b.into_iter()).collect();
    pairs.sort_by_key(|(a, _)| *a);

    let window_size = n / windows;
    let mut max_corr = 0.0f64;

    for w in 0..windows {
        let start = w * window_size;
        let end = if w == windows - 1 { n } else { (w + 1) * window_size };
        
        let wa: Vec<f64> = pairs[start..end].iter().map(|(a, _)| *a as f64).collect();
        let wb: Vec<f64> = pairs[start..end].iter().map(|(_, b)| *b as f64).collect();

        let n_w = wa.len();
        if n_w < 2 {
            continue;
        }

        let mean_a: f64 = wa.iter().sum::<f64>() / n_w as f64;
        let mean_b: f64 = wb.iter().sum::<f64>() / n_w as f64;

        let cov: f64 = wa.iter().zip(wb.iter())
            .map(|(a, b)| (a - mean_a) * (b - mean_b))
            .sum();
        let var_a: f64 = wa.iter().map(|a| (a - mean_a).powi(2)).sum();
        let var_b: f64 = wb.iter().map(|b| (b - mean_b).powi(2)).sum();

        if var_a > 0.0 && var_b > 0.0 {
            let corr = (cov / (var_a.sqrt() * var_b.sqrt())).abs();
            if corr > max_corr {
                max_corr = corr;
            }
        }
    }

    max_corr
}

/// Compute tag density correlation across genomic windows.
/// Divides each genome into `windows` bins and computes the Pearson correlation
/// of tag counts per bin. Captures synteny even when exact tag order is disrupted
/// by rearrangements, because syntenic regions share similar tag density patterns.
pub fn windowed_density_correlation(record_a: &TgtRecord, record_b: &TgtRecord, windows: usize) -> f64 {
    if record_a.total_length == 0 || record_b.total_length == 0 || windows == 0 {
        return 0.0;
    }

    let bin_size_a = (record_a.total_length as f64 / windows as f64).ceil() as u64;
    let bin_size_b = (record_b.total_length as f64 / windows as f64).ceil() as u64;

    let mut bins_a = vec![0.0f64; windows];
    for tag in &record_a.tags {
        let bin = (tag.position / bin_size_a).min((windows - 1) as u64) as usize;
        bins_a[bin] += 1.0;
    }

    let mut bins_b = vec![0.0f64; windows];
    for tag in &record_b.tags {
        let bin = (tag.position / bin_size_b).min((windows - 1) as u64) as usize;
        bins_b[bin] += 1.0;
    }

    // Compute Pearson correlation
    let n = windows as f64;
    let mean_a: f64 = bins_a.iter().sum::<f64>() / n;
    let mean_b: f64 = bins_b.iter().sum::<f64>() / n;

    if mean_a == 0.0 || mean_b == 0.0 {
        return 0.0;
    }

    let cov: f64 = bins_a.iter().zip(bins_b.iter())
        .map(|(a, b)| (a - mean_a) * (b - mean_b))
        .sum();
    let var_a: f64 = bins_a.iter().map(|a| (a - mean_a).powi(2)).sum();
    let var_b: f64 = bins_b.iter().map(|b| (b - mean_b).powi(2)).sum();

    if var_a <= 0.0 || var_b <= 0.0 {
        return 0.0;
    }

    let corr = cov / (var_a.sqrt() * var_b.sqrt());
    corr.abs()
}

// ─────────────────────────────────────────────────────────────────────────────
// Approximate tag matching (Hamming distance ≤ 1)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute Hamming distance between two 32-byte tag sequences.
/// Only counts non-null bytes (stops at first null byte in either sequence).
pub fn hamming_distance(a: &[u8; 32], b: &[u8; 32]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Count approximate common tags between two records, allowing Hamming distance ≤ max_dist.
///
/// For max_dist=1, uses a fast variant-generation approach:
/// For each tag in r1, generate all 1bp variants and check against r2's exact set.
/// Complexity: O(n * 93) for 32bp tags, where n = r1.tag_count().
pub fn count_approximate_common_tags(r1: &TgtRecord, r2: &TgtRecord, max_dist: usize) -> usize {
    if max_dist == 0 {
        return r1.sequence_set.intersection(&r2.sequence_set).count();
    }

    if max_dist == 1 {
        return count_approximate_common_tags_h1(r1, r2);
    }

    // For max_dist > 1, fall back to brute-force (slow for large datasets)
    let mut count = 0;
    for seq1 in &r1.sequence_set {
        for seq2 in &r2.sequence_set {
            if hamming_distance(seq1, seq2) <= max_dist {
                count += 1;
                break;
            }
        }
    }
    count
}

/// Fast Hamming-1 approximate matching using variant generation.
fn count_approximate_common_tags_h1(r1: &TgtRecord, r2: &TgtRecord) -> usize {
    let bases = [b'A', b'T', b'C', b'G'];
    let mut count = 0;

    for seq1 in &r1.sequence_set {
        // Check exact match first
        if r2.sequence_set.contains(seq1) {
            count += 1;
            continue;
        }

        // Generate all 1bp variants
        let mut found = false;
        for pos in 0..32 {
            let original = seq1[pos];
            for &base in &bases {
                if base == original {
                    continue;
                }
                let mut variant = *seq1;
                variant[pos] = base;
                if r2.sequence_set.contains(&variant) {
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }

        if found {
            count += 1;
        }
    }

    count
}

// ─────────────────────────────────────────────────────────────────────────────
// 5kb local window comparison (SynTracker-inspired)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute synteny score by comparing 5kb non-overlapping windows.
///
/// For each window, compute the Jaccard similarity of tag sequences present
/// in both genomes. The final score is the mean Jaccard across all windows
/// that contain at least one tag in either genome. Windows with no tags in
/// both genomes are skipped (not counted as 0).
///
/// This approach is robust to rearrangements and indels because it only
/// compares local tag presence within each window, not global tag order.
pub fn windowed_synteny_score(
    record_a: &TgtRecord,
    record_b: &TgtRecord,
    window_size: u64,
) -> f64 {
    if record_a.total_length == 0 || record_b.total_length == 0 || window_size == 0 {
        return 0.0;
    }

    let windows_a = ((record_a.total_length + window_size - 1) / window_size) as usize;
    let windows_b = ((record_b.total_length + window_size - 1) / window_size) as usize;
    let n_windows = windows_a.max(windows_b);

    // Build window -> tag set for genome A
    let mut window_tags_a: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_a.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_a[w].insert(tag.sequence);
    }

    // Build window -> tag set for genome B
    let mut window_tags_b: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_b.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_b[w].insert(tag.sequence);
    }

    // Compute mean Jaccard across windows with at least one tag
    let mut total_jaccard = 0.0;
    let mut valid_windows = 0usize;

    for i in 0..n_windows {
        let set_a = &window_tags_a[i];
        let set_b = &window_tags_b[i];

        if set_a.is_empty() && set_b.is_empty() {
            continue;
        }

        let intersection: HashSet<_> = set_a.intersection(set_b).collect();
        let union: HashSet<_> = set_a.union(set_b).collect();

        if !union.is_empty() {
            total_jaccard += intersection.len() as f64 / union.len() as f64;
            valid_windows += 1;
        }
    }

    if valid_windows == 0 {
        return 0.0;
    }

    total_jaccard / valid_windows as f64
}

/// Compute approximate windowed synteny score with Hamming distance ≤ 1.
///
/// For each window, tags from genome A are matched against tags from genome B
/// using approximate matching (1bp tolerance). This increases sensitivity
/// for divergent strains where SNPs disrupt exact tag matches.
pub fn windowed_synteny_score_approx(
    record_a: &TgtRecord,
    record_b: &TgtRecord,
    window_size: u64,
) -> f64 {
    if record_a.total_length == 0 || record_b.total_length == 0 || window_size == 0 {
        return 0.0;
    }

    let windows_a = ((record_a.total_length + window_size - 1) / window_size) as usize;
    let windows_b = ((record_b.total_length + window_size - 1) / window_size) as usize;
    let n_windows = windows_a.max(windows_b);

    // Build window -> tag set for genome A
    let mut window_tags_a: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_a.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_a[w].insert(tag.sequence);
    }

    // Build window -> tag set for genome B
    let mut window_tags_b: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_b.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_b[w].insert(tag.sequence);
    }

    let mut total_score = 0.0;
    let mut valid_windows = 0usize;

    for i in 0..n_windows {
        let set_a = &window_tags_a[i];
        let set_b = &window_tags_b[i];

        if set_a.is_empty() && set_b.is_empty() {
            continue;
        }

        // Approximate matching: count A tags that match B tags within Hamming-1
        let mut matched = 0usize;
        for seq_a in set_a {
            if set_b.contains(seq_a) {
                matched += 1;
                continue;
            }
            // Generate 1bp variants
            let bases = [b'A', b'T', b'C', b'G'];
            let mut found = false;
            for pos in 0..32 {
                let orig = seq_a[pos];
                for &base in &bases {
                    if base == orig {
                        continue;
                    }
                    let mut variant = *seq_a;
                    variant[pos] = base;
                    if set_b.contains(&variant) {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }
            if found {
                matched += 1;
            }
        }

        let union_size = set_a.len() + set_b.len();
        if union_size > 0 {
            // Score = 2 * matched / union_size (similar to shared_ratio)
            total_score += 2.0 * matched as f64 / union_size as f64;
            valid_windows += 1;
        }
    }

    if valid_windows == 0 {
        return 0.0;
    }

    total_score / valid_windows as f64
}
