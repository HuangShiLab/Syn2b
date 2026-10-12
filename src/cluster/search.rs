//! Query-side gene-cluster search and stable report output.

use crate::cluster::features::parse_features;
use crate::cluster::features::token_weight;
use crate::cluster::index::ClusterIndex;
use crate::cluster::model::{build_clusters, GeneCluster};
use crate::cluster::scoring::{W_ADJACENCY, W_COPY_NUMBER, W_FUNCTIONAL, W_PROJECTED_ORDER};
use anyhow::Result;
use parking_lot::Mutex;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
pub struct ClusterSearchOptions {
    pub min_functional_jaccard: f64,
    pub min_synteny: f64,
    pub max_targets: usize,
}

impl Default for ClusterSearchOptions {
    fn default() -> Self {
        Self {
            min_functional_jaccard: 0.25,
            min_synteny: 0.10,
            max_targets: 25,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ClusterHit {
    pub query_cluster_id: String,
    pub target_cluster_id: String,
    pub target_genome: String,
    pub target_contig: String,
    pub functional_jaccard: f64,
    pub weighted_functional_similarity: f64,
    pub adjacency_score: f64,
    pub projected_order_score: f64,
    pub copy_number_score: f64,
    pub cluster_score: f64,
    pub missing_query_genes: usize,
    pub extra_target_genes: usize,
    pub inversion_flag: bool,
    pub split_cluster_flag: bool,
}

/// Search-time scoring avoids the general `score_pair` allocations. Results are
/// numerically identical; the first MVP used the general path.
fn fast_score(
    query: &GeneCluster,
    query_weighted: &BTreeMap<String, f64>,
    query_adjacencies: &BTreeSet<(String, String)>,
    query_representatives: &[String],
    query_total_count: u64,
    target: &GeneCluster,
    target_representatives: &[String],
) -> (f64, f64, f64, f64, f64, f64, bool) {
    let mut intersection = 0usize;
    let mut union = 0usize;
    let mut weighted_intersection = 0.0f64;
    let mut weighted_union = 0.0f64;
    for token in query.token_counts.keys().chain(target.token_counts.keys()) {
        let left_query = query.token_counts.contains_key(token);
        let right_target = target.token_counts.contains_key(token);
        intersection += usize::from(left_query && right_target);
        union += usize::from(left_query || right_target);
        let left = query_weighted.get(token).copied().unwrap_or_default();
        let right = token_weight(token)
            * target.token_counts.get(token).copied().unwrap_or_default() as f64;
        weighted_intersection += left.min(right);
        weighted_union += left.max(right);
    }
    let functional_jaccard = if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    };
    let weighted_functional_similarity = if weighted_union <= f64::EPSILON {
        0.0
    } else {
        weighted_intersection / weighted_union
    };

    let target_adjacencies: BTreeSet<(String, String)> = target
        .adjacencies
        .iter()
        .map(|(a, b)| {
            if a <= b {
                (a.clone(), b.clone())
            } else {
                (b.clone(), a.clone())
            }
        })
        .collect();
    let query_adjacency_count = query_adjacencies.len();
    let target_adjacency_count = target_adjacencies.len();
    let adjacency_intersection = query_adjacencies.intersection(&target_adjacencies).count();
    let adjacency_union = query_adjacencies.union(&target_adjacencies).count();
    let adjacency_score = if query_adjacency_count == 0 && target_adjacency_count == 0 {
        1.0
    } else if adjacency_union == 0 {
        0.0
    } else {
        adjacency_intersection as f64 / adjacency_union as f64
    };

    let target_tokens: BTreeSet<&String> = target.token_counts.keys().collect();
    let query_tokens: BTreeSet<&String> = query.token_counts.keys().collect();
    let projected_query: Vec<String> = query_representatives
        .iter()
        .filter(|token| target_tokens.contains(token))
        .cloned()
        .collect();
    let projected_target: Vec<String> = target_representatives
        .iter()
        .filter(|token| query_tokens.contains(token))
        .cloned()
        .collect();
    let denominator = 1
        .max(projected_query.len().saturating_sub(1))
        .max(1.max(projected_target.len().saturating_sub(1)));
    let forward = lcs(&projected_query, &projected_target);
    let reverse = lcs(
        &projected_query.iter().rev().cloned().collect::<Vec<_>>(),
        &projected_target,
    );
    let projected_order_score = (forward.max(reverse) as f64 / denominator as f64).clamp(0.0, 1.0);
    let inversion_flag = reverse > forward.saturating_mul(3 / 2) && reverse >= 2;

    let mut copy_distance = 0u64;
    for token in query.token_counts.keys().chain(target.token_counts.keys()) {
        let left = query.token_counts.get(token).copied().unwrap_or_default();
        let right = target.token_counts.get(token).copied().unwrap_or_default();
        copy_distance += left.abs_diff(right) as u64;
    }
    let copy_scale = 2 * query_total_count
        .max(
            target
                .token_counts
                .values()
                .map(|count| *count as u64)
                .sum(),
        )
        .max(1);
    let copy_number_score = (1.0 - copy_distance as f64 / copy_scale as f64).clamp(0.0, 1.0);
    let cluster_score = (W_FUNCTIONAL * weighted_functional_similarity
        + W_ADJACENCY * adjacency_score
        + W_PROJECTED_ORDER * projected_order_score
        + W_COPY_NUMBER * copy_number_score)
        .clamp(0.0, 1.0);
    (
        functional_jaccard,
        weighted_functional_similarity,
        adjacency_score,
        projected_order_score,
        copy_number_score,
        cluster_score,
        inversion_flag,
    )
}

fn lcs(a: &[String], b: &[String]) -> usize {
    if a.is_empty() || b.is_empty() {
        return 0;
    }
    let mut previous = vec![0usize; b.len() + 1];
    let mut current = vec![0usize; b.len() + 1];
    for item_a in a {
        for (j, item_b) in b.iter().enumerate() {
            current[j + 1] = if item_a == item_b {
                previous[j] + 1
            } else {
                previous[j + 1].max(current[j])
            };
        }
        std::mem::swap(&mut previous, &mut current);
        current.fill(0);
    }
    previous[b.len()]
}

#[derive(Debug)]
pub struct ClusterSearchCache {
    representatives: Mutex<BTreeMap<u32, Arc<[String]>>>,
}

impl ClusterSearchCache {
    pub fn new(_index: &ClusterIndex) -> Self {
        Self {
            representatives: Mutex::new(BTreeMap::new()),
        }
    }

    fn representatives(&self, index: &ClusterIndex, cluster_id: u32) -> Arc<[String]> {
        let mut cache = self.representatives.lock();
        cache
            .entry(cluster_id)
            .or_insert_with(|| {
                index.clusters[cluster_id as usize]
                    .gene_tokens
                    .iter()
                    .filter_map(|tokens| tokens.keys().next().cloned())
                    .collect::<Vec<_>>()
                    .into()
            })
            .clone()
    }
}

pub fn search_clusters(
    query: &GeneCluster,
    index: &ClusterIndex,
    options: ClusterSearchOptions,
) -> Vec<ClusterHit> {
    search_clusters_cached(query, index, &ClusterSearchCache::new(index), options)
}

pub fn search_clusters_cached(
    query: &GeneCluster,
    index: &ClusterIndex,
    cache: &ClusterSearchCache,
    options: ClusterSearchOptions,
) -> Vec<ClusterHit> {
    let mut candidates: BTreeSet<u32> = BTreeSet::new();
    // Seed with orthologue/domain evidence (weight >=2). Product keywords are
    // retained for scoring, but their high cross-cluster frequency makes them
    // poor primary seeds at 10k+ scale.
    for token in query.token_counts.keys() {
        if token_weight(token) < 2.0 {
            continue;
        }
        for id in index.token_postings.get(token).into_iter().flatten() {
            candidates.insert(*id);
        }
    }
    for (a, b) in &query.adjacencies {
        let key = if a <= b {
            (a.clone(), b.clone())
        } else {
            (b.clone(), a.clone())
        };
        for id in index.adjacency_postings.get(&key).into_iter().flatten() {
            candidates.insert(*id);
        }
    }

    let mut query_weighted = BTreeMap::new();
    for (token, count) in &query.token_counts {
        query_weighted.insert(token.clone(), token_weight(token) * *count as f64);
    }
    let query_adjacency_set: BTreeSet<(String, String)> = query
        .adjacencies
        .iter()
        .map(|(a, b)| {
            if a <= b {
                (a.clone(), b.clone())
            } else {
                (b.clone(), a.clone())
            }
        })
        .collect();
    let query_representatives: Vec<String> = query
        .gene_tokens
        .iter()
        .filter_map(|tokens| tokens.keys().next().cloned())
        .collect();
    let query_total_count = query.token_counts.values().map(|count| *count as u64).sum();
    let candidate_ids: Vec<u32> = candidates.into_iter().collect();
    let score_target = move |id: u32| -> Option<ClusterHit> {
        let target = &index.clusters[id as usize];
        let mut cheap_intersection = 0usize;
        let mut cheap_union = 0usize;
        for token in query.token_counts.keys().chain(target.token_counts.keys()) {
            let left = query.token_counts.contains_key(token);
            let right = target.token_counts.contains_key(token);
            cheap_intersection += usize::from(left && right);
            cheap_union += usize::from(left || right);
        }
        let cheap_jaccard = if cheap_union == 0 {
            0.0
        } else {
            cheap_intersection as f64 / cheap_union as f64
        };
        if cheap_jaccard < options.min_functional_jaccard {
            return None;
        }
        let (
            functional_jaccard,
            weighted_functional_similarity,
            adjacency_score,
            projected_order_score,
            copy_number_score,
            cluster_score,
            inversion_flag,
        ) = fast_score(
            query,
            &query_weighted,
            &query_adjacency_set,
            &query_representatives,
            query_total_count,
            target,
            cache.representatives(index, id).as_ref(),
        );
        if functional_jaccard < options.min_functional_jaccard
            || adjacency_score < options.min_synteny
        {
            return None;
        }
        let missing_query_genes = query
            .token_counts
            .keys()
            .filter(|token| !target.token_counts.contains_key(token.as_str()))
            .count();
        let extra_target_genes = target
            .token_counts
            .keys()
            .filter(|token| !query.token_counts.contains_key(token.as_str()))
            .count();
        let query_span = query.end.saturating_sub(query.start) + 1;
        let target_span = target.end.saturating_sub(target.start) + 1;
        let split_cluster_flag = query_span > 0 && target_span * 2 < query_span;
        Some(ClusterHit {
            query_cluster_id: query.cluster_id.clone(),
            target_cluster_id: target.cluster_id.clone(),
            target_genome: target.genome_id.clone(),
            target_contig: target.contig_id.clone(),
            functional_jaccard,
            weighted_functional_similarity,
            adjacency_score,
            projected_order_score,
            copy_number_score,
            cluster_score,
            missing_query_genes,
            extra_target_genes,
            inversion_flag,
            split_cluster_flag,
        })
    };
    // Rayon has per-call setup cost; exact-seed searches often have few targets.
    let mut hits: Vec<ClusterHit> = if candidate_ids.len() < 1024 {
        candidate_ids.into_iter().filter_map(score_target).collect()
    } else {
        candidate_ids
            .par_iter()
            .filter_map(|id| score_target(*id))
            .collect()
    };
    hits.sort_by(|a, b| {
        b.cluster_score
            .partial_cmp(&a.cluster_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.target_cluster_id.cmp(&b.target_cluster_id))
    });
    hits.truncate(options.max_targets.max(1));
    hits
}

pub fn search_feature_table(
    path: impl AsRef<std::path::Path>,
    index: &ClusterIndex,
    max_gap: u64,
    min_genes: usize,
    options: ClusterSearchOptions,
) -> Result<Vec<ClusterHit>> {
    let clusters = build_clusters(parse_features(path)?, max_gap, min_genes)?;
    let cache = ClusterSearchCache::new(index);
    let mut hits = Vec::new();
    for cluster in &clusters {
        hits.extend(search_clusters_cached(cluster, index, &cache, options));
    }
    Ok(hits)
}

pub fn write_cluster_hits(path: impl AsRef<std::path::Path>, hits: &[ClusterHit]) -> Result<()> {
    let mut out = std::io::BufWriter::new(std::fs::File::create(path.as_ref())?);
    writeln!(out, "query_cluster_id\ttarget_cluster_id\ttarget_genome\ttarget_contig\tfunctional_jaccard\tweighted_functional_similarity\tadjacency_score\tprojected_order_score\tcopy_number_score\tcluster_score\tmissing_query_genes\textra_target_genes\tinversion_flag\tsplit_cluster_flag")?;
    for hit in hits {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            hit.query_cluster_id,
            hit.target_cluster_id,
            hit.target_genome,
            hit.target_contig,
            hit.functional_jaccard,
            hit.weighted_functional_similarity,
            hit.adjacency_score,
            hit.projected_order_score,
            hit.copy_number_score,
            hit.cluster_score,
            hit.missing_query_genes,
            hit.extra_target_genes,
            hit.inversion_flag,
            hit.split_cluster_flag,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::features::Feature;
    use crate::cluster::index::ClusterIndex;
    use crate::cluster::model::build_clusters;

    fn feature(id: &str, start: u64, end: u64, ko: &str, cluster: &str) -> Feature {
        Feature {
            gene_id: id.into(),
            genome_id: "genome".into(),
            contig_id: "c".into(),
            start,
            end,
            strand: '+',
            product: "conserved_function".into(),
            ko: ko.into(),
            pfam: String::new(),
            cog: String::new(),
            eggnog: String::new(),
            smash_tool: String::new(),
            cluster_id: Some(cluster.into()),
        }
    }

    #[test]
    fn exact_query_recovers_indexed_cluster() {
        let target = build_clusters(
            vec![
                feature("t1", 1, 10, "K001", "target"),
                feature("t2", 20, 30, "K002", "target"),
                feature("t3", 40, 50, "K003", "target"),
            ],
            1_000,
            2,
        )
        .unwrap();
        let index = ClusterIndex::build(target).unwrap();
        let query = build_clusters(
            vec![
                feature("q1", 1, 10, "K001", "query"),
                feature("q2", 20, 30, "K002", "query"),
                feature("q3", 40, 50, "K003", "query"),
            ],
            1_000,
            2,
        )
        .unwrap();
        let hits = search_clusters(&query[0], &index, ClusterSearchOptions::default());
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].target_cluster_id, "genome:target");
        assert!(hits[0].cluster_score > 0.99);
    }
}
