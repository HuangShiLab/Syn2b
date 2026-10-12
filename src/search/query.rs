//! Candidate lookup, chaining, scoring, and stable placement output model.

use crate::search::chaining::{chain_matches, ChainedBlock, TagMatch};
use crate::search::index::{canonical_sequence, contig_lengths, ContigIndex};
use crate::search::scoring::{
    confidence, placement_score, PlacementComponents, PlacementConfidence,
};
use crate::tgt::{Tag, TgtRecord};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Clone, Copy, Debug)]
pub struct SearchOptions {
    pub min_tags: usize,
    pub min_adjacencies: usize,
    pub max_targets: usize,
}

impl Default for SearchOptions {
    fn default() -> Self {
        Self {
            min_tags: 5,
            min_adjacencies: 3,
            max_targets: 25,
        }
    }
}

/// A scored placement for one query contig.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Placement {
    pub query_id: String,
    pub query_length: u64,
    pub query_tags: usize,
    pub query_adjacencies: usize,
    pub best_genome: String,
    pub best_species_cluster: String,
    pub best_contig: String,
    pub orientation: char,
    pub shared_tags: usize,
    pub shared_adjacencies: usize,
    pub chained_blocks: usize,
    pub query_coverage: f64,
    pub reference_support_fraction: f64,
    pub placement_score: f64,
    pub confidence: PlacementConfidence,
    pub secondary_genome: Option<String>,
    pub secondary_score: Option<f64>,
    pub chimera_warning: bool,
    pub best_block: Option<ChainedBlock>,
    /// Scored targets retained up to `max_targets`, in deterministic ranking.
    pub candidates: Vec<CandidateSummary>,
}

/// One ranked target retained for top-k benchmark reporting.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CandidateSummary {
    pub genome: String,
    pub species_cluster: String,
    pub score: f64,
    pub shared_tags: usize,
    pub shared_adjacencies: usize,
    pub orientation: char,
}

#[derive(Clone)]
struct QueryAdjacency {
    a_id: u32,
    b_id: u32,
}

#[derive(Default)]
struct CandidateSeeds {
    tags: BTreeSet<u32>,
    adjacencies: BTreeSet<(u32, u32)>,
}

#[derive(Debug)]
struct ScoredCandidate {
    genome_idx: u32,
    components: PlacementComponents,
    score: f64,
    block: ChainedBlock,
    shared_tags: usize,
    shared_adjacencies: usize,
    orientation_consistency: f64,
    reference_support: f64,
}

/// Search every contig in one query TGT record.
///
/// Results are ordered by query contig id and descending placement score.
pub fn search_record(
    record: &TgtRecord,
    index: &ContigIndex,
    options: SearchOptions,
) -> Vec<Placement> {
    let mut by_contig: BTreeMap<u16, Vec<&Tag>> = BTreeMap::new();
    for tag in &record.tags {
        by_contig.entry(tag.contig_id).or_default().push(tag);
    }

    let mut results = Vec::new();
    for (contig_id, query_tags) in by_contig {
        let mut tags = query_tags;
        tags.sort_by_key(|tag| (tag.enzyme.index(), tag.position));
        let mut placement = placement_for_tags(record, contig_id, &tags, index, options);
        placement.query_id = query_contig_name(record, contig_id);
        results.push(placement);
    }
    results
}

/// Search every contig in multiple query records.
pub fn search_records(
    records: &[TgtRecord],
    index: &ContigIndex,
    options: SearchOptions,
) -> Vec<Placement> {
    let mut results = Vec::new();
    for record in records {
        results.extend(search_record(record, index, options));
    }
    results
}

fn placement_for_tags(
    record: &TgtRecord,
    contig_id: u16,
    query_tags: &[&Tag],
    index: &ContigIndex,
    options: SearchOptions,
) -> Placement {
    let query_adjacencies = query_adjacency_set(index, query_tags);
    let query_length = query_contig_length(record, contig_id);
    let mut fallback = Placement {
        query_id: query_contig_name(record, contig_id),
        query_length,
        query_tags: query_tags.len(),
        query_adjacencies: query_adjacencies.len(),
        best_genome: "*".to_string(),
        best_species_cluster: "*".to_string(),
        best_contig: "*".to_string(),
        orientation: '+',
        shared_tags: 0,
        shared_adjacencies: 0,
        chained_blocks: 0,
        query_coverage: 0.0,
        reference_support_fraction: 0.0,
        placement_score: 0.0,
        confidence: PlacementConfidence::Exploratory,
        secondary_genome: None,
        secondary_score: None,
        chimera_warning: false,
        best_block: None,
        candidates: Vec::new(),
    };

    if query_tags.len() < options.min_tags {
        return fallback;
    }

    let seeds = candidate_seeds(index, query_tags, &query_adjacencies);
    let eligible: Vec<u32> = seeds
        .iter()
        .filter(|(_, seed)| {
            seed.tags.len() >= options.min_tags && seed.adjacencies.len() >= options.min_adjacencies
        })
        .map(|(&genome_idx, _)| genome_idx)
        .collect();
    if eligible.is_empty() {
        return fallback;
    }

    // Provisional score without uniqueness, then add a rank-aware uniqueness term.
    let mut scored: Vec<ScoredCandidate> = eligible
        .iter()
        .filter_map(|&genome_idx| {
            score_candidate(
                genome_idx,
                query_tags,
                &query_adjacencies,
                index,
                query_length,
            )
        })
        .collect();
    if scored.is_empty() {
        return fallback;
    }
    scored.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.genome_idx.cmp(&b.genome_idx))
    });
    for rank in 0..scored.len() {
        let uniqueness = if rank == 0 {
            if scored.len() == 1 {
                1.0
            } else {
                let best = scored[0].score;
                let second = scored[1].score;
                if best <= f64::EPSILON {
                    0.0
                } else {
                    (1.0 - (best - second) / best).clamp(0.0, 1.0)
                }
            }
        } else {
            0.0
        };
        scored[rank].components.uniqueness_penalty = uniqueness;
        scored[rank].score = placement_score(&scored[rank].components);
    }
    scored.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.block.n_adjacencies.cmp(&a.block.n_adjacencies))
            .then_with(|| b.block.n_tags.cmp(&a.block.n_tags))
            .then_with(|| a.genome_idx.cmp(&b.genome_idx))
            .then_with(|| {
                a.block
                    .reference_contig_id
                    .cmp(&b.block.reference_contig_id)
            })
    });
    scored.truncate(options.max_targets.max(1));

    let best = &scored[0];
    let genome = &index.genomes[best.genome_idx as usize];
    let contig_name = index
        .contigs
        .iter()
        .find(|c| c.genome_id == best.genome_idx && c.contig_id == best.block.reference_contig_id)
        .map(|c| c.contig_name.clone())
        .unwrap_or_else(|| best.block.reference_contig_id.to_string());

    fallback.best_genome = genome.genome_id.clone();
    fallback.best_species_cluster = genome
        .metadata
        .species_cluster
        .clone()
        .unwrap_or_else(|| genome.genome_id.clone());
    fallback.best_contig = contig_name;
    fallback.orientation = if best.block.orientation { '+' } else { '-' };
    fallback.shared_tags = best.shared_tags;
    fallback.shared_adjacencies = best.shared_adjacencies;
    fallback.chained_blocks = 1;
    fallback.query_coverage = best.components.chained_block_coverage;
    fallback.reference_support_fraction = best.reference_support;
    fallback.placement_score = best.score;
    fallback.confidence = confidence(
        best.score,
        best.shared_tags,
        best.shared_adjacencies,
        query_length,
    );
    fallback.secondary_genome = scored.get(1).map(|candidate| {
        index.genomes[candidate.genome_idx as usize]
            .genome_id
            .clone()
    });
    fallback.secondary_score = scored.get(1).map(|candidate| candidate.score);
    fallback.best_block = Some(best.block.clone());
    fallback.candidates = scored
        .iter()
        .map(|candidate| CandidateSummary {
            genome: index.genomes[candidate.genome_idx as usize]
                .genome_id
                .clone(),
            species_cluster: index.genomes[candidate.genome_idx as usize]
                .metadata
                .species_cluster
                .clone()
                .unwrap_or_else(|| {
                    index.genomes[candidate.genome_idx as usize]
                        .genome_id
                        .clone()
                }),
            score: candidate.score,
            shared_tags: candidate.shared_tags,
            shared_adjacencies: candidate.shared_adjacencies,
            orientation: if candidate.block.orientation {
                '+'
            } else {
                '-'
            },
        })
        .collect();
    fallback
}

fn candidate_seeds(
    index: &ContigIndex,
    query_tags: &[&Tag],
    query_adjacencies: &[QueryAdjacency],
) -> HashMap<u32, CandidateSeeds> {
    let mut candidates: HashMap<u32, CandidateSeeds> = HashMap::new();
    for tag in query_tags {
        let tag_id = index.sequence_id(canonical_sequence(&tag.sequence));
        let Some(tag_id) = tag_id else { continue };
        for posting in index.tag_postings.get(&tag_id).into_iter().flatten() {
            candidates
                .entry(posting.genome_id)
                .or_default()
                .tags
                .insert(tag_id);
        }
    }
    for adjacency in query_adjacencies {
        let key = if adjacency.a_id <= adjacency.b_id {
            (adjacency.a_id, adjacency.b_id)
        } else {
            (adjacency.b_id, adjacency.a_id)
        };
        for posting in index.adjacency_postings.get(&key).into_iter().flatten() {
            let seed = candidates.entry(posting.genome_id).or_default();
            seed.adjacencies.insert(key);
            seed.tags.insert(posting.tag_a_id);
            seed.tags.insert(posting.tag_b_id);
        }
    }
    candidates
}

fn score_candidate(
    genome_idx: u32,
    query_tags: &[&Tag],
    query_adjacencies: &[QueryAdjacency],
    index: &ContigIndex,
    query_length: u64,
) -> Option<ScoredCandidate> {
    // Repeated landmarks have many valid coordinate mappings and can create a
    // false long chain. Following the genome-synteny scorer, chain only
    // landmarks that are unique on both the query contig and candidate contig.
    let mut query_tag_counts: HashMap<u32, usize> = HashMap::new();
    let mut query_identity_count = 0usize;
    for tag in query_tags {
        let identity = canonical_sequence(&tag.sequence);
        if let Some(tag_id) = index.sequence_id(identity) {
            *query_tag_counts.entry(tag_id).or_default() += 1;
        }
        query_identity_count += 1;
    }
    let mut matches: Vec<TagMatch> = Vec::new();
    for tag in query_tags {
        // A landmark unique to the query/source may be absent from the target.
        // That is normal in source-excluded placement and must not discard the
        // entire candidate; only the missing anchor is skipped.
        let Some(tag_id) = index.sequence_id(canonical_sequence(&tag.sequence)) else {
            continue;
        };
        if query_tag_counts[&tag_id] != 1 {
            continue;
        }
        for posting in index.tag_postings.get(&tag_id).into_iter().flatten() {
            if posting.genome_id == genome_idx {
                matches.push(TagMatch {
                    query_position: tag.position,
                    reference_position: posting.position,
                    tag_id,
                    reference_contig_id: posting.contig_id,
                });
            }
        }
    }
    // Repeat filtering is per candidate reference contig. Genome-wide counts
    // would discard valid anchors merely because homologous loci occur on other
    // contigs; MAG contig fragmentation differs between source and target.
    let mut reference_counts: HashMap<(u32, u32), usize> = HashMap::new();
    for matched in &matches {
        *reference_counts
            .entry((matched.reference_contig_id, matched.tag_id))
            .or_default() += 1;
    }
    matches.retain(|matched| {
        reference_counts
            .get(&(matched.reference_contig_id, matched.tag_id))
            .copied()
            .unwrap_or(0)
            == 1
    });
    let mut best: Option<ScoredCandidate> = None;
    {
        for orientation in [true, false] {
            let chain = chain_matches(&matches, orientation);
            if chain.len() < 2 {
                continue;
            }
            let selected_ids: BTreeSet<u32> = chain.iter().map(|m| m.tag_id).collect();
            let mut shared_adjacencies = BTreeSet::new();
            let mut oriented = 0usize;
            let mut comparable = 0usize;
            for adjacency in query_adjacencies {
                let key = if adjacency.a_id <= adjacency.b_id {
                    (adjacency.a_id, adjacency.b_id)
                } else {
                    (adjacency.b_id, adjacency.a_id)
                };
                for posting in index.adjacency_postings.get(&key).into_iter().flatten() {
                    if posting.genome_id != genome_idx {
                        continue;
                    }
                    if !selected_ids.contains(&posting.tag_a_id)
                        || !selected_ids.contains(&posting.tag_b_id)
                    {
                        continue;
                    }
                    shared_adjacencies.insert(key);
                    comparable += 1;
                    let adjacency_forward =
                        posting.tag_a_id == adjacency.a_id && posting.tag_b_id == adjacency.b_id;
                    if adjacency_forward == orientation {
                        oriented += 1;
                    }
                }
            }

            // Report the contig contributing the largest part of the chain.
            let mut contig_counts: BTreeMap<u32, usize> = BTreeMap::new();
            for match_item in &chain {
                *contig_counts
                    .entry(match_item.reference_contig_id)
                    .or_default() += 1;
            }
            let best_reference_contig = contig_counts
                .iter()
                .max_by_key(|(contig_id, count)| (**count, std::cmp::Reverse(**contig_id)))
                .map(|(contig_id, _)| *contig_id)
                .unwrap_or_default();

            let query_start = chain.first().unwrap().query_position;
            let query_end = chain.last().unwrap().query_position;
            let reference_start = chain
                .iter()
                .map(|m| m.reference_position)
                .min()
                .unwrap_or_default();
            let reference_end = chain
                .iter()
                .map(|m| m.reference_position)
                .max()
                .unwrap_or_default();
            let block = ChainedBlock {
                query_start,
                query_end,
                reference_start,
                reference_end,
                orientation,
                reference_contig_id: best_reference_contig,
                n_tags: chain.len() as u64,
                n_adjacencies: shared_adjacencies.len() as u64,
            };
            let query_span = query_end.saturating_sub(query_start).saturating_add(1);
            let coverage = if query_length > 0 {
                query_span as f64 / query_length as f64
            } else {
                1.0
            };
            let orientation_consistency = if comparable == 0 {
                0.0
            } else {
                oriented as f64 / comparable as f64
            };
            let genome = &index.genomes[genome_idx as usize];
            let reference_support = if genome.total_tags == 0 {
                0.0
            } else {
                block.n_tags as f64 / genome.total_tags as f64
            };
            let components = PlacementComponents {
                shared_tag_fraction: block.n_tags as f64 / query_identity_count as f64,
                conserved_adjacency_fraction: if query_adjacencies.is_empty() {
                    0.0
                } else {
                    shared_adjacencies.len() as f64 / query_adjacencies.len() as f64
                },
                chained_block_coverage: coverage.clamp(0.0, 1.0),
                orientation_consistency,
                uniqueness_penalty: 0.0,
                contig_interior_support: if query_tags.len() <= 2 {
                    0.0
                } else {
                    ((block.n_tags.saturating_sub(2)) as f64
                        / (query_identity_count.saturating_sub(2)) as f64)
                        .clamp(0.0, 1.0)
                },
            };
            let candidate = ScoredCandidate {
                genome_idx,
                score: placement_score(&components),
                shared_tags: block.n_tags as usize,
                shared_adjacencies: shared_adjacencies.len(),
                orientation_consistency,
                reference_support,
                block,
                components,
            };
            let better = match &best {
                None => true,
                Some(current) => {
                    (
                        candidate.block.n_tags,
                        candidate.block.n_adjacencies,
                        ordered(candidate.orientation_consistency),
                        candidate.genome_idx,
                        candidate.block.reference_contig_id,
                    ) > (
                        current.block.n_tags,
                        current.block.n_adjacencies,
                        ordered(current.orientation_consistency),
                        current.genome_idx,
                        current.block.reference_contig_id,
                    )
                }
            };
            if better {
                best = Some(candidate);
            }
        }
    }
    best
}

// A private tuple helper keeps the comparison above explicit without exposing
// orientation as a comparable public numeric value.
fn ordered(value: f64) -> u64 {
    (value.clamp(0.0, 1.0) * 1_000_000.0).round() as u64
}

fn query_adjacency_set(index: &ContigIndex, query_tags: &[&Tag]) -> Vec<QueryAdjacency> {
    let mut groups: BTreeMap<(u16, u8), Vec<&Tag>> = BTreeMap::new();
    for tag in query_tags {
        groups
            .entry((tag.contig_id, tag.enzyme.index()))
            .or_default()
            .push(tag);
    }
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    for (_, mut group) in groups {
        group.sort_by_key(|tag| tag.position);
        let mut previous: Option<&Tag> = None;
        for tag in group {
            if let Some(previous_tag) = previous {
                let max_tag_len = std::cmp::max(
                    previous_tag.enzyme.properties().tag_length,
                    tag.enzyme.properties().tag_length,
                ) as u64;
                if tag.position > previous_tag.position
                    && tag.position - previous_tag.position > max_tag_len
                {
                    let a_id = index.sequence_id(canonical_sequence(&previous_tag.sequence));
                    let b_id = index.sequence_id(canonical_sequence(&tag.sequence));
                    if let (Some(a_id), Some(b_id)) = (a_id, b_id) {
                        // Ordered adjacency identity, not genomic instance, is
                        // the unit of synteny. Repeated operons therefore count
                        // once in both numerator and denominator.
                        if seen.insert((a_id, b_id)) {
                            out.push(QueryAdjacency { a_id, b_id });
                        }
                    }
                }
            }
            previous = Some(tag);
        }
    }
    out
}

fn query_contig_name(record: &TgtRecord, contig_id: u16) -> String {
    if contig_id == 0 {
        if record.contig_names.is_empty() {
            record.genome_id.clone()
        } else {
            record.contig_names[0].clone()
        }
    } else {
        record
            .contig_names
            .get(contig_id as usize - 1)
            .cloned()
            .unwrap_or_else(|| format!("contig_{}", contig_id))
    }
}

fn query_contig_length(record: &TgtRecord, contig_id: u16) -> u64 {
    let lengths = contig_lengths(record);
    if contig_id == 0 {
        lengths.first().copied().unwrap_or(record.total_length)
    } else {
        lengths.get(contig_id as usize - 1).copied().unwrap_or(0)
    }
}

/// Stable `placements.tsv` writer.
pub fn write_placements(
    path: impl AsRef<std::path::Path>,
    placements: &[Placement],
) -> anyhow::Result<()> {
    use std::io::Write;
    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(out, "query_id\tquery_length\tquery_tags\tquery_adjacencies\tbest_genome\tbest_species_cluster\tbest_contig\torientation\tshared_tags\tshared_adjacencies\tchained_blocks\tquery_coverage\treference_support_fraction\tplacement_score\tconfidence\tsecondary_genome\tsecondary_score\tchimera_warning")?;
    for p in placements {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            p.query_id,
            p.query_length,
            p.query_tags,
            p.query_adjacencies,
            p.best_genome,
            p.best_species_cluster,
            p.best_contig,
            p.orientation,
            p.shared_tags,
            p.shared_adjacencies,
            p.chained_blocks,
            p.query_coverage,
            p.reference_support_fraction,
            p.placement_score,
            p.confidence.as_str(),
            p.secondary_genome.as_deref().unwrap_or("*"),
            p.secondary_score
                .map_or("*".to_string(), |v| format!("{:.6}", v)),
            p.chimera_warning,
        )?;
    }
    Ok(())
}

/// Stable ranked candidate table for top-k benchmarking.
pub fn write_candidate_hits(
    path: impl AsRef<std::path::Path>,
    placements: &[Placement],
) -> anyhow::Result<()> {
    use std::io::Write;
    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        out,
        "query_id	rank	genome	species_cluster	score	shared_tags	shared_adjacencies	orientation"
    )?;
    for placement in placements {
        for (rank, candidate) in placement.candidates.iter().enumerate() {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}",
                placement.query_id,
                rank + 1,
                candidate.genome,
                candidate.species_cluster,
                candidate.score,
                candidate.shared_tags,
                candidate.shared_adjacencies,
                candidate.orientation,
            )?;
        }
    }
    Ok(())
}

/// Stable `placement_details.tsv` writer.
pub fn write_placement_details(
    path: impl AsRef<std::path::Path>,
    placements: &[Placement],
) -> anyhow::Result<()> {
    use std::io::Write;
    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(out, "query_id\tquery_block_id\tquery_start\tquery_end\treference_genome\treference_contig\treference_start\treference_end\torientation\tn_tags\tn_adjacencies\tblock_score")?;
    for placement in placements {
        let Some(block) = &placement.best_block else {
            continue;
        };
        writeln!(
            out,
            "{}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
            placement.query_id,
            block.query_start,
            block.query_end,
            placement.best_genome,
            placement.best_contig,
            block.reference_start,
            block.reference_end,
            if block.orientation { '+' } else { '-' },
            block.n_tags,
            block.n_adjacencies,
            placement.placement_score,
        )?;
    }
    Ok(())
}

/// Internal naming helper shared with chimera window output.
pub(crate) fn query_stub_name(record: &TgtRecord, contig_id: u16) -> String {
    query_contig_name(record, contig_id)
}
