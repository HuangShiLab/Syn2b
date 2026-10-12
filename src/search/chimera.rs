//! Sliding-window chimera segmentation.
//!
//! The first implementation reports every window's dominant genome/cluster and
//! marks adjacent-window transitions. Breakpoint localization is therefore
//! bounded by the requested window size; posterior calibration is future work.

use crate::search::query::{query_stub_name, search_record, Placement, SearchOptions};
use crate::tgt::TgtRecord;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug)]
pub struct ChimeraOptions {
    pub window: u64,
    pub min_window_tags: usize,
    pub max_targets: usize,
    /// Minimum min(best, second-best) score required to call a transition.
    /// Calibrated on simulated 250-500 kb + 250-500 kb OAPGC chimeras.
    pub min_transition_confidence: f64,
}

impl Default for ChimeraOptions {
    fn default() -> Self {
        Self {
            window: 50_000,
            min_window_tags: 10,
            max_targets: 25,
            min_transition_confidence: 0.20,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ChimeraSegment {
    pub query_id: String,
    pub segment: usize,
    pub start: u64,
    pub end: u64,
    pub best_cluster: String,
    pub best_genome: String,
    pub score: f64,
    pub second_cluster: String,
    pub second_score: f64,
    pub cluster_transition: bool,
    pub chimera_confidence: f64,
}

/// Scan each contig independently in fixed coordinate windows.
pub fn scan_record(
    record: &TgtRecord,
    index: &crate::search::index::ContigIndex,
    options: ChimeraOptions,
) -> Vec<ChimeraSegment> {
    if options.window == 0 {
        return Vec::new();
    }
    let mut by_contig: BTreeMap<u16, Vec<&crate::tgt::Tag>> = BTreeMap::new();
    for tag in &record.tags {
        by_contig.entry(tag.contig_id).or_default().push(tag);
    }

    let mut segments = Vec::new();
    for (contig_id, mut tags) in by_contig {
        tags.sort_by_key(|tag| tag.position);
        let contig_start = tags.first().map_or(0, |tag| tag.position);
        let contig_end = tags.last().map_or(0, |tag| tag.position);
        let mut segment = 1usize;
        let mut previous_cluster: Option<String> = None;
        let mut start = contig_start;
        while start <= contig_end {
            let end = start.saturating_add(options.window).saturating_sub(1);
            let window_tags: Vec<&crate::tgt::Tag> = tags
                .iter()
                .copied()
                .filter(|tag| tag.position >= start && tag.position <= end)
                .collect();
            let mut window_record =
                TgtRecord::new(&record.genome_id, end.saturating_sub(start) + 1);
            window_record.contig_names = record.contig_names.clone();
            window_record.contig_offsets = record.contig_offsets.clone();
            window_record.contig_circular = record.contig_circular.clone();
            for tag in window_tags {
                window_record.add_tag(tag.clone());
            }
            let search_options = SearchOptions {
                min_tags: options.min_window_tags,
                min_adjacencies: 0,
                max_targets: options.max_targets,
            };
            let hits = search_record(&window_record, index, search_options);
            let dominant = hits.iter().max_by(|a, b| {
                b.placement_score
                    .partial_cmp(&a.placement_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.best_genome.cmp(&b.best_genome))
            });
            let empty = Placement {
                query_id: String::new(),
                query_length: 0,
                query_tags: 0,
                query_adjacencies: 0,
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
                confidence: crate::search::scoring::PlacementConfidence::Exploratory,
                secondary_genome: None,
                secondary_score: None,
                chimera_warning: false,
                best_block: None,
                candidates: Vec::new(),
            };
            let hit = dominant.unwrap_or(&empty);
            let raw_transition = previous_cluster
                .as_ref()
                .map(|previous| previous != &hit.best_species_cluster)
                .unwrap_or(false);
            let inferred_confidence = if raw_transition {
                hit.placement_score.min(hit.secondary_score.unwrap_or(0.0))
            } else {
                0.0
            };
            let transition =
                raw_transition && inferred_confidence >= options.min_transition_confidence;
            segments.push(ChimeraSegment {
                query_id: format!("{}:{}-{}", query_stub_name(record, contig_id), start, end),
                segment,
                start,
                end,
                best_cluster: hit.best_species_cluster.clone(),
                best_genome: hit.best_genome.clone(),
                score: hit.placement_score,
                second_cluster: hit
                    .secondary_genome
                    .as_ref()
                    .and_then(|genome_id| index.genome_ids.get(genome_id))
                    .and_then(|genome_idx| index.genomes.get(*genome_idx as usize))
                    .and_then(|genome| genome.metadata.species_cluster.clone())
                    .unwrap_or_else(|| "*".to_string()),
                second_score: hit.secondary_score.unwrap_or(0.0),
                cluster_transition: transition,
                chimera_confidence: inferred_confidence,
            });
            previous_cluster = Some(hit.best_species_cluster.clone());
            start = end.saturating_add(1);
            segment += 1;
        }
    }
    segments
}

pub fn write_chimerism(
    path: impl AsRef<std::path::Path>,
    segments: &[ChimeraSegment],
) -> anyhow::Result<()> {
    use std::io::Write;
    let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(out, "query_id\tsegment\tstart\tend\tbest_cluster\tbest_genome\tscore\tsecond_cluster\tsecond_score\tcluster_transition\tchimera_confidence")?;
    for segment in segments {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{:.6}",
            segment.query_id,
            segment.segment,
            segment.start,
            segment.end,
            segment.best_cluster,
            segment.best_genome,
            segment.score,
            segment.second_cluster,
            segment.second_score,
            segment.cluster_transition,
            segment.chimera_confidence,
        )?;
    }
    Ok(())
}
