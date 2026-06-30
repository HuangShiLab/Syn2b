//! Tag adjacency graph — core data structure for synteny detection
//!
//! Inspired by ntSynt's minimizer graph approach. Each node represents a unique
//! 2bRAD tag sequence, and edges connect tags that are adjacent in at least one
//! genome. Linear paths through the graph correspond to syntenic backbones.

use crate::tgt::record::TgtRecord;
use crate::tgt::tag::{Strand, Tag};
use std::collections::{HashMap, HashSet};

// ─────────────────────────────────────────────────────────────────────────────
// Data structures
// ─────────────────────────────────────────────────────────────────────────────

/// A node in the tag adjacency graph — represents a tag shared across genomes.
/// Each node stores the unique tag sequence and, for every genome that contains
/// this exact sequence, its genomic position and strand.
pub struct TagNode {
    pub tag_id: u64,
    pub sequence: [u8; 32],
    /// genome_id → (position, strand)
    pub positions: HashMap<String, (u64, Strand)>,
}

impl TagNode {
    /// Create a new TagNode with the given id and sequence.
    fn new(tag_id: u64, sequence: [u8; 32]) -> Self {
        Self {
            tag_id,
            sequence,
            positions: HashMap::new(),
        }
    }

    /// Return true if this node is shared by at least two genomes.
    pub fn is_shared(&self) -> bool {
        self.positions.len() >= 2
    }
}

impl std::fmt::Debug for TagNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TagNode")
            .field("tag_id", &self.tag_id)
            .field("sequence", &std::str::from_utf8(&self.sequence).unwrap_or("???"))
            .field("genomes", &self.positions.len())
            .finish()
    }
}

/// An edge connecting two adjacent tags in at least one genome.
/// The weight equals the number of distinct genomes where the source tag
/// and target tag appear consecutively (in that order).
pub struct AdjacencyEdge {
    pub source: u64, // Source tag ID
    pub target: u64, // Target tag ID
    pub weight: u32, // Number of genomes supporting this adjacency
    /// genome_id → inter-tag distance (gap size in bp)
    pub distances: HashMap<String, u32>,
}

impl AdjacencyEdge {
    fn new(source: u64, target: u64) -> Self {
        Self {
            source,
            target,
            weight: 0,
            distances: HashMap::new(),
        }
    }

    /// Add a genome's support to this edge.
    fn add_genome(&mut self, genome_id: String, distance: u32) {
        if !self.distances.contains_key(&genome_id) {
            self.weight += 1;
        }
        self.distances.insert(genome_id, distance);
    }
}

impl std::fmt::Debug for AdjacencyEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AdjacencyEdge")
            .field("source", &self.source)
            .field("target", &self.target)
            .field("weight", &self.weight)
            .finish()
    }
}

/// Tag adjacency graph for synteny detection.
///
/// The graph is built incrementally:
/// 1. `add_genome` ingests each genome's TGT record, creating nodes for unique tag sequences.
/// 2. `build_edges` creates edges between consecutive tags within each genome.
/// 3. `simplify` removes low-support noise.
/// 4. `linear_paths` extracts the syntenic backbones.
pub struct TagAdjacencyGraph {
    pub nodes: HashMap<u64, TagNode>,
    pub edges: HashMap<(u64, u64), AdjacencyEdge>,
    pub num_genomes: usize,
    /// genome_id → ordered list of tag IDs as they appear in that genome
    genome_order: HashMap<String, Vec<u64>>,
    /// sequence → tag_id lookup for deduplication
    sequence_to_id: HashMap<[u8; 32], u64>,
    next_tag_id: u64,
}

impl TagAdjacencyGraph {
    /// Create a new empty graph.
    pub fn new() -> Self {
        Self {
            nodes: HashMap::new(),
            edges: HashMap::new(),
            num_genomes: 0,
            genome_order: HashMap::new(),
            sequence_to_id: HashMap::new(),
            next_tag_id: 0,
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Genome ingestion
    // ─────────────────────────────────────────────────────────────────────────

    /// Add a genome's TGT record to the graph.
    ///
    /// For each tag in the record:
    /// - If the exact sequence already exists as a node, add this genome's position.
    /// - Otherwise, create a new TagNode.
    ///
    /// Tags are stored in the order they appear in the record (which must be
    /// sorted by position).
    pub fn add_genome(&mut self, genome_id: &str, record: &TgtRecord) {
        self.num_genomes += 1;
        let mut ordered_ids: Vec<u64> = Vec::with_capacity(record.tags.len());

        for tag in &record.tags {
            let tag_id = *self.sequence_to_id.entry(tag.sequence).or_insert_with(|| {
                let id = self.next_tag_id;
                self.next_tag_id += 1;
                let node = TagNode::new(id, tag.sequence);
                self.nodes.insert(id, node);
                id
            });

            // Add this genome's position to the node
            if let Some(node) = self.nodes.get_mut(&tag_id) {
                node.positions
                    .insert(genome_id.to_string(), (tag.position, tag.strand));
            }

            ordered_ids.push(tag_id);
        }

        self.genome_order
            .insert(genome_id.to_string(), ordered_ids);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Edge construction
    // ─────────────────────────────────────────────────────────────────────────

    /// Build adjacency edges from the ordering of tags within each genome.
    ///
    /// For each genome, iterate through its tags in position order. For each
    /// adjacent pair (tag_i, tag_{i+1}), create or update a directed edge.
    /// The edge weight counts how many genomes share this exact adjacency.
    /// The edge also stores the inter-tag distance for each supporting genome.
    pub fn build_edges(&mut self) {
        self.edges.clear();

        for (genome_id, ordered_ids) in &self.genome_order {
            for window in ordered_ids.windows(2) {
                let src = window[0];
                let tgt = window[1];
                let edge_key = (src, tgt);

                // Compute inter-tag distance using stored positions
                let distance = if let (Some(src_node), Some(tgt_node)) =
                    (self.nodes.get(&src), self.nodes.get(&tgt))
                {
                    if let (Some((src_pos, _)), Some((tgt_pos, _))) =
                        (src_node.positions.get(genome_id), tgt_node.positions.get(genome_id))
                    {
                        (*tgt_pos - *src_pos) as u32
                    } else {
                        0
                    }
                } else {
                    0
                };

                self.edges
                    .entry(edge_key)
                    .or_insert_with(|| AdjacencyEdge::new(src, tgt))
                    .add_genome(genome_id.clone(), distance);
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Graph simplification
    // ─────────────────────────────────────────────────────────────────────────

    /// Simplify the graph using ntSynt-inspired anchored node classification:
    ///
    /// 1. Remove all edges with weight < `min_weight`.
    /// 2. Remove all nodes with degree 0 (isolated — no incident edges remain).
    ///
    /// Nodes that are shared across multiple genomes but fall on low-support
    /// edges are effectively "unanchored" and removed, cleaning up the graph
    /// for backbone extraction.
    pub fn simplify(&mut self, min_weight: u32) {
        // Step 1: Remove edges below weight threshold
        self.edges
            .retain(|_, edge| edge.weight >= min_weight);

        // Step 2: Compute node degrees from remaining edges
        let mut degree: HashMap<u64, usize> = HashMap::new();
        for ((src, tgt), _) in &self.edges {
            *degree.entry(*src).or_insert(0) += 1;
            *degree.entry(*tgt).or_insert(0) += 1;
        }

        // Step 3: Remove isolated nodes (degree 0)
        let isolated: Vec<u64> = self
            .nodes
            .keys()
            .filter(|id| *degree.get(*id).unwrap_or(&0) == 0)
            .copied()
            .collect();

        for tag_id in &isolated {
            self.nodes.remove(tag_id);
            self.sequence_to_id
                .retain(|_, id| id != tag_id);
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Linear path extraction
    // ─────────────────────────────────────────────────────────────────────────

    /// Extract linear paths (degree-2 chains) = synteny backbones.
    ///
    /// A linear path is a maximal chain of nodes where each interior node
    /// has exactly two incident edges (one in, one out). These correspond to
    /// conserved syntenic regions across genomes.
    ///
    /// Returns each path as an ordered list of tag IDs.
    pub fn linear_paths(&self) -> Vec<Vec<u64>> {
        if self.nodes.is_empty() || self.edges.is_empty() {
            return Vec::new();
        }

        // Build adjacency list (undirected for path finding)
        let mut adj: HashMap<u64, Vec<u64>> = HashMap::new();
        for ((src, tgt), _) in &self.edges {
            adj.entry(*src).or_default().push(*tgt);
            adj.entry(*tgt).or_default().push(*src);
        }

        let mut visited: HashSet<u64> = HashSet::new();
        let mut paths: Vec<Vec<u64>> = Vec::new();

        // Find paths starting from endpoints (degree != 2) or unvisited degree-2 nodes
        let mut start_candidates: Vec<u64> = adj
            .keys()
            .filter(|&node| {
                let deg = adj.get(node).map_or(0, |v| v.len());
                deg != 2 || visited.contains(node)
            })
            .copied()
            .collect();

        // If all nodes are degree-2 (pure cycle component), pick arbitrary starts
        if start_candidates.is_empty() && !adj.is_empty() {
            start_candidates = adj.keys().copied().take(1).collect();
        }

        // Extend from each endpoint
        for &start in &start_candidates {
            if visited.contains(&start) {
                continue;
            }

            let mut path = vec![start];
            visited.insert(start);

            // Extend forward
            let mut current = start;
            loop {
                let neighbors = adj.get(&current).map_or(&[][..], |v| v.as_slice());
                let mut extended = false;
                for &next in neighbors {
                    if !visited.contains(&next) {
                        // For degree-2 nodes, continue the chain
                        let next_deg = adj.get(&next).map_or(0, |v| v.len());
                        path.push(next);
                        visited.insert(next);
                        current = next;
                        extended = true;
                        if next_deg != 2 {
                            break;
                        }
                    }
                }
                if !extended {
                    break;
                }
            }

            // Also try extending backward from start if start had degree 1
            let start_deg = adj.get(&start).map_or(0, |v| v.len());
            if start_deg == 1 {
                // Already went forward; the path is complete
            } else if start_deg > 2 {
                // Start is a branch point; we only captured one branch
            }

            if path.len() >= 2 {
                paths.push(path);
            } else {
                // Singleton — don't include as a synteny backbone
            }
        }

        // Handle remaining unvisited degree-2 nodes (cyclic components)
        for &node in adj.keys() {
            if visited.contains(&node) {
                continue;
            }
            let deg = adj.get(&node).map_or(0, |v| v.len());
            if deg == 2 {
                // Extract the cycle
                let mut cycle = Vec::new();
                let mut current = node;
                loop {
                    visited.insert(current);
                    cycle.push(current);
                    let neighbors = adj.get(&current).unwrap();
                    let unvisited_next = neighbors.iter().find(|&&n| !visited.contains(&n));
                    match unvisited_next {
                        Some(&next) => current = next,
                        None => break,
                    }
                }
                if cycle.len() >= 2 {
                    paths.push(cycle);
                }
            }
        }

        paths
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Common tag queries
    // ─────────────────────────────────────────────────────────────────────────

    /// Find tags shared across **all** genomes (exact sequence match).
    ///
    /// Returns the tag IDs of nodes whose `positions` map contains every
    /// ingested genome.
    pub fn find_common_tags(&self) -> Vec<u64> {
        self.nodes
            .values()
            .filter(|node| node.positions.len() == self.num_genomes)
            .map(|node| node.tag_id)
            .collect()
    }

    /// Find tags shared across genomes with Hamming tolerance.
    ///
    /// Tags are grouped by sequence similarity: two tags are in the same group
    /// if their 32bp sequences differ by at most `tolerance` positions. Only
    /// groups containing tags from all genomes are returned.
    ///
    /// This is an O(N²) operation over unique tag sequences; for large graphs
    /// a more efficient clustering approach (e.g., locality-sensitive hashing)
    /// would be warranted.
    pub fn find_common_tags_tolerance(&self, tolerance: u8) -> Vec<u64> {
        if self.nodes.is_empty() || self.num_genomes == 0 {
            return Vec::new();
        }

        // Collect all (sequence, tag_id, genome_count) tuples
        let node_list: Vec<(u64, [u8; 32], usize)> = self
            .nodes
            .values()
            .map(|n| (n.tag_id, n.sequence, n.positions.len()))
            .collect();

        // Union-find to cluster similar sequences
        let n = node_list.len();
        let mut parent: Vec<usize> = (0..n).collect();
        let mut rank: Vec<usize> = vec![0; n];

        fn find(parent: &mut [usize], x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        fn union(parent: &mut [usize], rank: &mut [usize], x: usize, y: usize) {
            let rx = find(parent, x);
            let ry = find(parent, y);
            if rx == ry {
                return;
            }
            match rank[rx].cmp(&rank[ry]) {
                std::cmp::Ordering::Less => parent[rx] = ry,
                std::cmp::Ordering::Greater => parent[ry] = rx,
                std::cmp::Ordering::Equal => {
                    parent[ry] = rx;
                    rank[rx] += 1;
                }
            }
        }

        // Compare all pairs of tag sequences
        for i in 0..n {
            for j in (i + 1)..n {
                let dist = hamming_distance_32(&node_list[i].1, &node_list[j].1);
                if dist <= tolerance {
                    union(&mut parent, &mut rank, i, j);
                }
            }
        }

        // Group by cluster and check if cluster spans all genomes
        let mut cluster_genomes: HashMap<usize, HashSet<String>> = HashMap::new();
        let mut cluster_representative: HashMap<usize, u64> = HashMap::new();

        // Build index: tag_id → index in node_list
        let id_to_idx: HashMap<u64, usize> = node_list
            .iter()
            .enumerate()
            .map(|(i, (id, _, _))| (*id, i))
            .collect();

        for (node_idx, &(tag_id, _, _)) in node_list.iter().enumerate() {
            let root = find(&mut parent, node_idx);
            cluster_representative.entry(root).or_insert(tag_id);

            // Collect all genome IDs for this node
            if let Some(node) = self.nodes.get(&tag_id) {
                let entry = cluster_genomes.entry(root).or_insert_with(HashSet::new);
                for genome_id in node.positions.keys() {
                    entry.insert(genome_id.clone());
                }
            }
        }

        // Return one representative tag_id per cluster that spans all genomes
        cluster_genomes
            .into_iter()
            .filter(|(_, genomes)| genomes.len() == self.num_genomes)
            .filter_map(|(root, _)| cluster_representative.get(&root).copied())
            .collect()
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Helpers
    // ─────────────────────────────────────────────────────────────────────────

    /// Compute the degree of each node in the current graph.
    pub(crate) fn node_degrees(&self) -> HashMap<u64, usize> {
        let mut degree: HashMap<u64, usize> = HashMap::new();
        for &(src, tgt) in self.edges.keys() {
            *degree.entry(src).or_insert(0) += 1;
            *degree.entry(tgt).or_insert(0) += 1;
        }
        degree
    }

    /// Get the ordered tag IDs for a given genome.
    pub(crate) fn genome_order(&self, genome_id: &str) -> Option<&Vec<u64>> {
        self.genome_order.get(genome_id)
    }

    /// Get the number of unique tag nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Get the number of adjacency edges.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }
}

impl Default for TagAdjacencyGraph {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility functions
// ─────────────────────────────────────────────────────────────────────────────

/// Hamming distance between two 32-byte sequences.
fn hamming_distance_32(a: &[u8; 32], b: &[u8; 32]) -> u8 {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count() as u8
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

    /// Helper: create a tag with the given index (encodes index into sequence bytes)
    fn make_tag(idx: u8, position: u64, enzyme: EnzymeType, strand: Strand) -> Tag {
        let mut seq = [b'A'; 32];
        seq[0] = idx; // Encode the index into first byte for dedup
        seq[1] = idx.wrapping_add(1);
        Tag::new(seq, position, enzyme, strand)
    }

    /// Helper: create a TGT record with evenly spaced tags
    fn make_record(genome_id: &str, tag_count: usize) -> TgtRecord {
        let mut record = TgtRecord::new(genome_id, 1_000_000);
        for i in 0..tag_count {
            let tag = make_tag(i as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward);
            record.add_tag(tag);
        }
        record
    }

    #[test]
    fn test_new_graph_is_empty() {
        let graph = TagAdjacencyGraph::new();
        assert_eq!(graph.node_count(), 0);
        assert_eq!(graph.edge_count(), 0);
        assert_eq!(graph.num_genomes, 0);
    }

    #[test]
    fn test_add_single_genome() {
        let mut graph = TagAdjacencyGraph::new();
        let record = make_record("genome_a", 5);
        graph.add_genome("genome_a", &record);

        assert_eq!(graph.num_genomes, 1);
        assert_eq!(graph.node_count(), 5);
    }

    #[test]
    fn test_add_two_genomes_shared_tags() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5); // Same tags

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);

        // 5 unique tag sequences, each shared by 2 genomes
        assert_eq!(graph.node_count(), 5);
        assert_eq!(graph.num_genomes, 2);

        // Each node should have 2 positions
        for node in graph.nodes.values() {
            assert_eq!(node.positions.len(), 2);
        }
    }

    #[test]
    fn test_build_edges_single_genome() {
        let mut graph = TagAdjacencyGraph::new();
        let record = make_record("genome_a", 5);
        graph.add_genome("genome_a", &record);
        graph.build_edges();

        // 5 tags in order → 4 edges
        assert_eq!(graph.edge_count(), 4);
        for edge in graph.edges.values() {
            assert_eq!(edge.weight, 1);
        }
    }

    #[test]
    fn test_build_edges_two_identical_genomes() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();

        // Same tag order in both genomes → edges should have weight 2
        assert_eq!(graph.edge_count(), 4);
        for edge in graph.edges.values() {
            assert_eq!(edge.weight, 2);
        }
    }

    #[test]
    fn test_simplify_removes_low_weight_edges() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();

        // All edges have weight 2; min_weight=3 should remove all
        graph.simplify(3);
        assert_eq!(graph.edge_count(), 0);
        // Nodes with degree 0 should be removed
        assert_eq!(graph.node_count(), 0);
    }

    #[test]
    fn test_simplify_keeps_high_weight_edges() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();

        // min_weight=2 should keep all edges (weight=2)
        graph.simplify(2);
        assert_eq!(graph.edge_count(), 4);
        assert_eq!(graph.node_count(), 5);
    }

    #[test]
    fn test_linear_paths_two_genomes() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let paths = graph.linear_paths();
        // Should have one linear path with all 5 nodes
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].len(), 5);
    }

    #[test]
    fn test_find_common_tags() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);

        let common = graph.find_common_tags();
        // All 5 tags are shared across both genomes
        assert_eq!(common.len(), 5);
    }

    #[test]
    fn test_find_common_tags_not_all_shared() {
        let mut graph = TagAdjacencyGraph::new();
        // Genome A: tags 0, 1, 2, 3, 4
        let record_a = make_record("genome_a", 5);
        graph.add_genome("genome_a", &record_a);

        // Genome B: only tags 0, 1 (different sequence indices)
        let mut record_b = TgtRecord::new("genome_b", 1_000_000);
        record_b.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(1, 200, EnzymeType::BcgI, Strand::Forward));
        graph.add_genome("genome_b", &record_b);

        let common = graph.find_common_tags();
        // Only tags 0 and 1 are in both genomes
        assert_eq!(common.len(), 2);
    }

    #[test]
    fn test_find_common_tags_tolerance_exact() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);

        // tolerance=0 means exact match — same as find_common_tags
        let common = graph.find_common_tags_tolerance(0);
        assert_eq!(common.len(), 5);
    }

    #[test]
    fn test_empty_graph_linear_paths() {
        let graph = TagAdjacencyGraph::new();
        let paths = graph.linear_paths();
        assert!(paths.is_empty());
    }

    /// Test that a known inversion (reversed tag order) is detected as
    /// a break in linear paths.
    #[test]
    fn test_inversion_detection() {
        let mut graph = TagAdjacencyGraph::new();

        // Genome A: tags in order 0,1,2,3,4
        let record_a = make_record("genome_a", 5);
        graph.add_genome("genome_a", &record_a);

        // Genome B: tags in order 0,1 (then inverted segment 4,3,2)
        let mut record_b = TgtRecord::new("genome_b", 1_000_000);
        record_b.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(1, 200, EnzymeType::BcgI, Strand::Forward));
        // Inverted segment: 4, 3, 2
        record_b.add_tag(make_tag(4, 300, EnzymeType::BcgI, Strand::Reverse));
        record_b.add_tag(make_tag(3, 400, EnzymeType::BcgI, Strand::Reverse));
        record_b.add_tag(make_tag(2, 500, EnzymeType::BcgI, Strand::Reverse));
        graph.add_genome("genome_b", &record_b);

        graph.build_edges();
        graph.simplify(1);

        // The inversion breaks the adjacency at positions 1→2 and 2→3.
        // Expected edges (both genomes agree):
        //   0→1 (weight 2)
        // Edges only in genome A: 1→2, 2→3, 3→4
        // Edges only in genome B: 1→4, 4→3, 3→2
        // After simplify(1), all edges with weight ≥ 1 remain.
        //
        // The graph should show a break: path [0,1] then a separate path [2,3,4]
        // (or similar depending on edge directions).
        let paths = graph.linear_paths();

        // There should be multiple paths because the inversion breaks the backbone
        assert!(
            paths.len() >= 2,
            "Inversion should produce at least 2 separate paths, got {:?}",
            paths
        );
    }
}
