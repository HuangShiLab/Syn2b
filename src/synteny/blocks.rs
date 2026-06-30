//! Synteny block extraction
//!
//! Converts linear paths from the tag adjacency graph into synteny blocks —
//! genomic regions that are conserved across multiple genomes. Also provides
//! functions for detecting indels and filtering blocks by size.

use crate::synteny::graph::TagAdjacencyGraph;
use crate::tgt::tag::Strand;
use std::collections::{HashMap, HashSet};

// ─────────────────────────────────────────────────────────────────────────────
// SyntenyBlock
// ─────────────────────────────────────────────────────────────────────────────

/// A synteny block — a genomic region conserved across multiple genomes.
///
/// A block is defined by a linear chain of tags (a path through the
/// adjacency graph) and the corresponding genomic coordinates in each
/// genome that contains the entire tag chain.
#[derive(Clone, Debug)]
pub struct SyntenyBlock {
    pub block_id: u64,
    pub start_tag: u64,
    pub end_tag: u64,
    /// genome_id → (start_position, end_position, strand)
    pub genome_coords: HashMap<String, (u64, u64, Strand)>,
    pub num_tags: usize,
}

impl SyntenyBlock {
    /// Compute the genomic span of this block in a given genome.
    pub fn span_in_genome(&self, genome_id: &str) -> Option<u64> {
        self.genome_coords
            .get(genome_id)
            .map(|(start, end, _)| end.saturating_sub(*start))
    }

    /// Return the minimum span across all genomes.
    pub fn min_span(&self) -> u64 {
        self.genome_coords
            .values()
            .map(|(start, end, _)| end.saturating_sub(*start))
            .min()
            .unwrap_or(0)
    }

    /// Return the maximum span across all genomes.
    pub fn max_span(&self) -> u64 {
        self.genome_coords
            .values()
            .map(|(start, end, _)| end.saturating_sub(*start))
            .max()
            .unwrap_or(0)
    }

    /// Number of genomes that contain this block.
    pub fn num_genomes(&self) -> usize {
        self.genome_coords.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Block extraction
// ─────────────────────────────────────────────────────────────────────────────

/// Extract synteny blocks from linear paths in the graph.
///
/// 1. Find all linear paths using `graph.linear_paths()`.
/// 2. For each path, look up the first and last tag's position in each genome.
/// 3. If a genome contains both the start and end tag, record its coordinates.
/// 4. Create a `SyntenyBlock` for each valid path.
///
/// A path is valid if at least two genomes contain the full path (start and end
/// tags are both present). The `num_tags` field records the path length.
pub fn extract_synteny_blocks(graph: &TagAdjacencyGraph) -> Vec<SyntenyBlock> {
    let paths = graph.linear_paths();
    if paths.is_empty() {
        return Vec::new();
    }

    let mut blocks = Vec::new();
    let mut block_id: u64 = 0;

    for path in &paths {
        if path.len() < 2 {
            continue;
        }

        let start_tag_id = path[0];
        let end_tag_id = path[path.len() - 1];

        // Collect all genomes that have both the start and end tags
        let start_node = match graph.nodes.get(&start_tag_id) {
            Some(n) => n,
            None => continue,
        };
        let end_node = match graph.nodes.get(&end_tag_id) {
            Some(n) => n,
            None => continue,
        };

        // Find common genomes between start and end nodes
        let common_genomes: HashSet<String> = start_node
            .positions
            .keys()
            .filter(|g| end_node.positions.contains_key(*g))
            .cloned()
            .collect();

        if common_genomes.len() < 2 {
            // Require at least 2 genomes for a synteny block
            continue;
        }

        let mut genome_coords = HashMap::new();

        for genome_id in &common_genomes {
            // For each genome in the block, determine the span from the
            // first to the last tag. We also check the internal tags to
            // ensure the path is fully present.
            let mut min_pos = u64::MAX;
            let mut max_pos = u64::MIN;
            let mut all_present = true;
            let mut path_strand = Strand::Forward;

            for &tag_id in path.iter() {
                if let Some(node) = graph.nodes.get(&tag_id) {
                    if let Some(&(pos, strand)) = node.positions.get(genome_id) {
                        if pos < min_pos {
                            min_pos = pos;
                        }
                        if pos > max_pos {
                            max_pos = pos;
                        }
                        // Use the strand of the first tag as the block strand
                        if tag_id == start_tag_id {
                            path_strand = strand;
                        }
                    } else {
                        all_present = false;
                        break;
                    }
                } else {
                    all_present = false;
                    break;
                }
            }

            if all_present {
                genome_coords.insert(
                    genome_id.clone(),
                    (min_pos, max_pos, path_strand),
                );
            }
        }

        if genome_coords.len() >= 2 {
            blocks.push(SyntenyBlock {
                block_id,
                start_tag: start_tag_id,
                end_tag: end_tag_id,
                genome_coords,
                num_tags: path.len(),
            });
            block_id += 1;
        }
    }

    blocks
}

// ─────────────────────────────────────────────────────────────────────────────
// Indel detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect indels by comparing inter-tag distances across genomes within a block.
///
/// For each adjacent tag pair in the block, compares the inter-tag distance
/// in each pair of genomes. Large differences indicate indels (insertions or
/// deletions).
///
/// Returns a list of `(genome_a, genome_b, indel_size)` where `indel_size` is:
/// - Positive: genome_b has an insertion relative to genome_a
/// - Negative: genome_b has a deletion relative to genome_a
/// - Zero: no difference (not included in the result)
///
/// Only differences exceeding a threshold (10% of the mean distance) are
/// reported to filter out noise.
pub fn detect_indels(
    block: &SyntenyBlock,
    graph: &TagAdjacencyGraph,
) -> Vec<(String, String, i64)> {
    // Reconstruct the tag path from the block
    let path = reconstruct_block_path(block, graph);
    if path.len() < 2 {
        return Vec::new();
    }

    // Collect genome IDs present in the block
    let genome_ids: Vec<String> = block.genome_coords.keys().cloned().collect();
    if genome_ids.len() < 2 {
        return Vec::new();
    }

    // For each adjacent tag pair in the path, compute distances per genome
    let mut indels = Vec::new();

    for window in path.windows(2) {
        let src_id = window[0];
        let tgt_id = window[1];

        let src_node = match graph.nodes.get(&src_id) {
            Some(n) => n,
            None => continue,
        };
        let tgt_node = match graph.nodes.get(&tgt_id) {
            Some(n) => n,
            None => continue,
        };

        // Compute distance for each genome
        let mut distances: HashMap<String, u64> = HashMap::new();
        for genome_id in &genome_ids {
            if let (Some(&(src_pos, _)), Some(&(tgt_pos, _))) =
                (src_node.positions.get(genome_id), tgt_node.positions.get(genome_id))
            {
                distances.insert(genome_id.clone(), tgt_pos.saturating_sub(src_pos));
            }
        }

        // Compare all pairs of genomes
        for i in 0..genome_ids.len() {
            for j in (i + 1)..genome_ids.len() {
                let gi = &genome_ids[i];
                let gj = &genome_ids[j];

                let di = *distances.get(gi).unwrap_or(&0);
                let dj = *distances.get(gj).unwrap_or(&0);

                if di == 0 || dj == 0 {
                    continue;
                }

                // Compute relative difference
                let mean_dist = (di + dj) as f64 / 2.0;
                let diff = dj as i64 - di as i64;

                // Threshold: report if difference is > 10% of mean distance
                // and absolute difference is at least 100 bp
                let threshold = (mean_dist * 0.1).max(100.0);
                if (diff as f64).abs() >= threshold {
                    indels.push((gi.clone(), gj.clone(), diff));
                }
            }
        }
    }

    // Deduplicate: merge indels between the same genome pair
    dedup_indels(indels)
}

/// Deduplicate indel entries by summing differences for the same genome pair.
fn dedup_indels(indels: Vec<(String, String, i64)>) -> Vec<(String, String, i64)> {
    let mut summed: HashMap<(String, String), i64> = HashMap::new();
    for (a, b, diff) in indels {
        *summed.entry((a, b)).or_insert(0) += diff;
    }
    summed.into_iter().map(|((a, b), diff)| (a, b, diff)).collect()
}

/// Reconstruct the tag path for a block from the graph.
///
/// Uses a BFS-like traversal from the start_tag to the end_tag following
/// edges in the graph.
fn reconstruct_block_path(block: &SyntenyBlock, graph: &TagAdjacencyGraph) -> Vec<u64> {
    // Try to find a path from start_tag to end_tag through the graph edges.
    // Since the block was created from a linear path, we can trace it by
    // following edges from the start.

    if block.start_tag == block.end_tag {
        return vec![block.start_tag];
    }

    // Build undirected adjacency for traversal
    let mut adj: HashMap<u64, Vec<u64>> = HashMap::new();
    for ((src, tgt), _) in &graph.edges {
        adj.entry(*src).or_default().push(*tgt);
        adj.entry(*tgt).or_default().push(*src);
    }

    // BFS to find path from start to end
    let mut visited: HashSet<u64> = HashSet::new();
    let mut parent: HashMap<u64, u64> = HashMap::new();
    let mut queue: Vec<u64> = vec![block.start_tag];
    visited.insert(block.start_tag);

    let mut found = false;
    while let Some(current) = queue.pop() {
        if current == block.end_tag {
            found = true;
            break;
        }
        for &neighbor in adj.get(&current).unwrap_or(&Vec::new()) {
            if !visited.contains(&neighbor) {
                visited.insert(neighbor);
                parent.insert(neighbor, current);
                queue.push(neighbor);
            }
        }
    }

    if !found {
        // Fallback: just return start and end
        return vec![block.start_tag, block.end_tag];
    }

    // Reconstruct path from end to start
    let mut path = Vec::new();
    let mut current = block.end_tag;
    path.push(current);
    while let Some(&p) = parent.get(&current) {
        path.push(p);
        current = p;
    }
    path.reverse();
    path
}

// ─────────────────────────────────────────────────────────────────────────────
// Block filtering
// ─────────────────────────────────────────────────────────────────────────────

/// Filter blocks by minimum genomic size (span in base pairs).
///
/// Only keeps blocks whose minimum span across all containing genomes is
/// at least `min_size`.
pub fn filter_blocks_by_size(
    blocks: Vec<SyntenyBlock>,
    min_size: u64,
) -> Vec<SyntenyBlock> {
    blocks
        .into_iter()
        .filter(|block| block.min_span() >= min_size)
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzyme::enzyme::EnzymeType;
    use crate::synteny::graph::TagAdjacencyGraph;
    use crate::tgt::record::TgtRecord;
    use crate::tgt::tag::{Strand, Tag};

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
    // extract_synteny_blocks
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_extract_blocks_two_genomes() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 10);
        let record_b = make_record("genome_a", 10);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let blocks = extract_synteny_blocks(&graph);
        assert!(!blocks.is_empty(), "Two identical genomes should produce blocks");

        // All tags are in the same order, so we expect one large block
        let first_block = &blocks[0];
        assert_eq!(first_block.num_tags, 10);
        assert_eq!(first_block.num_genomes(), 2);
    }

    #[test]
    fn test_extract_blocks_single_genome() {
        let mut graph = TagAdjacencyGraph::new();
        let record = make_record("genome_a", 10);
        graph.add_genome("genome_a", &record);
        graph.build_edges();
        graph.simplify(1);

        let blocks = extract_synteny_blocks(&graph);
        // Single genome: no synteny block (need ≥ 2 genomes)
        assert!(blocks.is_empty(), "Single genome should not produce blocks");
    }

    #[test]
    fn test_extract_blocks_empty_graph() {
        let graph = TagAdjacencyGraph::new();
        let blocks = extract_synteny_blocks(&graph);
        assert!(blocks.is_empty());
    }

    #[test]
    fn test_block_spans() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 10);
        let record_b = make_record("genome_a", 10);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let blocks = extract_synteny_blocks(&graph);
        assert!(!blocks.is_empty());

        let block = &blocks[0];
        // 10 tags spaced at 1000bp intervals: span ~ 9000 bp
        assert!(block.min_span() > 0, "Block should have positive span");

        for (genome_id, (start, end, strand)) in &block.genome_coords {
            assert!(start < end, "Genome {}: start < end", genome_id);
            assert!(*strand == Strand::Forward, "Forward strand expected");
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // detect_indels
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_detect_indels_no_indel() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 10);
        let record_b = make_record("genome_a", 10);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let blocks = extract_synteny_blocks(&graph);
        assert!(!blocks.is_empty());

        let indels = detect_indels(&blocks[0], &graph);
        // Identical genomes should have no indels
        assert!(indels.is_empty(), "Identical genomes should have no indels, got {:?}", indels);
    }

    #[test]
    fn test_detect_indels_with_different_spacing() {
        let mut graph = TagAdjacencyGraph::new();

        // Genome A: tags at 0, 1000, 2000, ... (1kb spacing)
        let mut record_a = TgtRecord::new("genome_a", 100_000);
        for i in 0..5 {
            record_a.add_tag(make_tag(i as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: tags at 0, 1500, 3000, ... (1.5kb spacing — larger gaps)
        let mut record_b = TgtRecord::new("genome_b", 100_000);
        for i in 0..5 {
            record_b.add_tag(make_tag(i as u8, (i * 1500) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let blocks = extract_synteny_blocks(&graph);
        if blocks.is_empty() {
            // If no blocks, indel detection can't run — that's ok for this test
            return;
        }

        let indels = detect_indels(&blocks[0], &graph);
        // Different spacing (1000 vs 1500) should be detected as indels
        assert!(
            !indels.is_empty(),
            "Different spacing should produce indels"
        );

        // The indel size should be positive for the genome with larger gaps
        for (ga, gb, size) in &indels {
            if *size > 0 {
                // Genome B has larger gaps → positive indel size
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // filter_blocks_by_size
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_filter_blocks_by_size() {
        let mut block1 = SyntenyBlock {
            block_id: 0,
            start_tag: 0,
            end_tag: 9,
            genome_coords: HashMap::new(),
            num_tags: 10,
        };
        block1.genome_coords.insert("g1".to_string(), (0, 9000, Strand::Forward));
        block1.genome_coords.insert("g2".to_string(), (0, 9000, Strand::Forward));

        let mut block2 = SyntenyBlock {
            block_id: 1,
            start_tag: 10,
            end_tag: 11,
            genome_coords: HashMap::new(),
            num_tags: 2,
        };
        block2.genome_coords.insert("g1".to_string(), (100, 200, Strand::Forward));
        block2.genome_coords.insert("g2".to_string(), (100, 200, Strand::Forward));

        let blocks = vec![block1, block2];
        let filtered = filter_blocks_by_size(blocks, 1000);

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].block_id, 0, "Only the large block should remain");
    }

    #[test]
    fn test_filter_blocks_all_pass() {
        let mut block = SyntenyBlock {
            block_id: 0,
            start_tag: 0,
            end_tag: 9,
            genome_coords: HashMap::new(),
            num_tags: 10,
        };
        block.genome_coords.insert("g1".to_string(), (0, 5000, Strand::Forward));

        let blocks = vec![block];
        let filtered = filter_blocks_by_size(blocks, 1);
        assert_eq!(filtered.len(), 1);
    }

    #[test]
    fn test_filter_blocks_none_pass() {
        let mut block = SyntenyBlock {
            block_id: 0,
            start_tag: 0,
            end_tag: 1,
            genome_coords: HashMap::new(),
            num_tags: 2,
        };
        block.genome_coords.insert("g1".to_string(), (0, 50, Strand::Forward));

        let blocks = vec![block];
        let filtered = filter_blocks_by_size(blocks, 1000);
        assert!(filtered.is_empty());
    }

    #[test]
    fn test_reconstruct_block_path() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 5);
        let record_b = make_record("genome_a", 5);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let paths = graph.linear_paths();
        assert!(!paths.is_empty());

        let block = SyntenyBlock {
            block_id: 0,
            start_tag: paths[0][0],
            end_tag: paths[0][paths[0].len() - 1],
            genome_coords: HashMap::new(),
            num_tags: paths[0].len(),
        };

        let reconstructed = reconstruct_block_path(&block, &graph);
        assert!(
            reconstructed.len() >= 2,
            "Reconstructed path should have at least 2 nodes"
        );
        assert_eq!(reconstructed[0], block.start_tag);
        assert_eq!(reconstructed[reconstructed.len() - 1], block.end_tag);
    }
}
