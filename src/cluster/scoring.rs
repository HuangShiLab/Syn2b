//! Functional, adjacency, projected-order, and copy-number cluster scoring.

use crate::cluster::features::token_weight;
use crate::cluster::model::GeneCluster;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const W_FUNCTIONAL: f64 = 0.45;
pub const W_ADJACENCY: f64 = 0.30;
pub const W_PROJECTED_ORDER: f64 = 0.15;
pub const W_COPY_NUMBER: f64 = 0.10;

fn weighted_counts(cluster: &GeneCluster) -> BTreeMap<String, f64> {
    cluster
        .token_counts
        .iter()
        .map(|(token, count)| (token.clone(), token_weight(token) * *count as f64))
        .collect()
}

pub fn functional_scores(query: &GeneCluster, target: &GeneCluster) -> (f64, f64) {
    let q: BTreeSet<&String> = query.token_counts.keys().collect();
    let t: BTreeSet<&String> = target.token_counts.keys().collect();
    let intersection = q.intersection(&t).count();
    let union = q.union(&t).count();
    let unweighted = if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    };

    let qw = weighted_counts(query);
    let tw = weighted_counts(target);
    let mut intersection_weight = 0.0;
    let mut union_weight = 0.0;
    for token in qw.keys().chain(tw.keys()).collect::<BTreeSet<_>>() {
        let left = qw.get(token).copied().unwrap_or_default();
        let right = tw.get(token).copied().unwrap_or_default();
        intersection_weight += left.min(right);
        union_weight += left.max(right);
    }
    let weighted = if union_weight <= f64::EPSILON {
        0.0
    } else {
        intersection_weight / union_weight
    };
    (unweighted, weighted.clamp(0.0, 1.0))
}

fn canonical_adjacencies(cluster: &GeneCluster) -> BTreeSet<(String, String)> {
    cluster
        .adjacencies
        .iter()
        .map(|(a, b)| {
            if a <= b {
                (a.clone(), b.clone())
            } else {
                (b.clone(), a.clone())
            }
        })
        .collect()
}

pub fn adjacency_score(query: &GeneCluster, target: &GeneCluster) -> f64 {
    let q = canonical_adjacencies(query);
    let t = canonical_adjacencies(target);
    if q.is_empty() && t.is_empty() {
        return 1.0;
    }
    let intersection = q.intersection(&t).count();
    let union = q.union(&t).count();
    if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    }
}

fn representative_tokens(cluster: &GeneCluster) -> Vec<String> {
    cluster
        .gene_tokens
        .iter()
        .filter_map(|tokens| tokens.keys().next().cloned())
        .collect()
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

/// Gene-aware order score: remove genes absent from the other cluster before
/// computing LCS, then normalize by the possible projected adjacencies.
pub fn projected_order_score(query: &GeneCluster, target: &GeneCluster) -> (f64, bool) {
    let target_tokens: BTreeSet<&String> = target.token_counts.keys().collect();
    let query_tokens: BTreeSet<&String> = query.token_counts.keys().collect();
    let q: Vec<String> = representative_tokens(query)
        .into_iter()
        .filter(|token| target_tokens.contains(token))
        .collect();
    let t: Vec<String> = representative_tokens(target)
        .into_iter()
        .filter(|token| query_tokens.contains(token))
        .collect();
    let denominator = 1
        .max(q.len().saturating_sub(1))
        .max(1.max(t.len().saturating_sub(1)));
    let forward = lcs(&q, &t);
    let reverse = lcs(&q.into_iter().rev().collect::<Vec<_>>(), &t);
    let score = (forward.max(reverse) as f64 / denominator as f64).clamp(0.0, 1.0);
    let inversion = reverse > forward.saturating_mul(3 / 2) && reverse >= 2;
    (score, inversion)
}

pub fn copy_number_score(query: &GeneCluster, target: &GeneCluster) -> f64 {
    let mut tokens: BTreeSet<&String> = query.token_counts.keys().collect();
    tokens.extend(target.token_counts.keys());
    let mut distance = 0u64;
    for token in tokens {
        let left = query.token_counts.get(token).copied().unwrap_or_default();
        let right = target.token_counts.get(token).copied().unwrap_or_default();
        distance += left.abs_diff(right) as u64;
    }
    let scale = 2 * query
        .token_counts
        .values()
        .chain(target.token_counts.values())
        .map(|count| *count as u64)
        .sum::<u64>()
        .max(1);
    (1.0 - distance as f64 / scale as f64).clamp(0.0, 1.0)
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ClusterScores {
    pub functional_jaccard: f64,
    pub weighted_functional_similarity: f64,
    pub adjacency_score: f64,
    pub projected_order_score: f64,
    pub copy_number_score: f64,
    pub cluster_score: f64,
    pub inversion_flag: bool,
}

pub fn score_pair(query: &GeneCluster, target: &GeneCluster) -> ClusterScores {
    let (functional_jaccard, weighted_functional_similarity) = functional_scores(query, target);
    let adjacency_score = adjacency_score(query, target);
    let (projected_order_score, inversion_flag) = projected_order_score(query, target);
    let copy_number_score = copy_number_score(query, target);
    let cluster_score = W_FUNCTIONAL * weighted_functional_similarity
        + W_ADJACENCY * adjacency_score
        + W_PROJECTED_ORDER * projected_order_score
        + W_COPY_NUMBER * copy_number_score;
    ClusterScores {
        functional_jaccard,
        weighted_functional_similarity,
        adjacency_score,
        projected_order_score,
        copy_number_score,
        cluster_score: cluster_score.clamp(0.0, 1.0),
        inversion_flag,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::features::Feature;
    use crate::cluster::model::build_clusters;

    fn feature(id: &str, start: u64, end: u64, ko: &str) -> Feature {
        Feature {
            gene_id: id.into(),
            genome_id: "g".into(),
            contig_id: "c".into(),
            start,
            end,
            strand: '+',
            product: format!("p_{id}"),
            ko: ko.into(),
            pfam: String::new(),
            cog: String::new(),
            eggnog: String::new(),
            smash_tool: String::new(),
            cluster_id: Some("q".into()),
        }
    }

    fn cluster(features: Vec<Feature>) -> GeneCluster {
        build_clusters(features, 1_000, 2).unwrap().remove(0)
    }

    #[test]
    fn identical_clusters_score_one() {
        let genes = vec![
            feature("a", 1, 10, "K001"),
            feature("b", 20, 30, "K002"),
            feature("c", 40, 50, "K003"),
        ];
        let q = cluster(genes.clone());
        let t = cluster(genes);
        let scores = score_pair(&q, &t);
        assert!(scores.weighted_functional_similarity > 0.999);
        assert!(scores.adjacency_score > 0.999);
        assert!(scores.projected_order_score > 0.999);
        assert!(scores.cluster_score > 0.999);
    }

    #[test]
    fn middle_gene_loss_retains_projected_order() {
        let q = cluster(vec![
            feature("a", 1, 10, "K001"),
            feature("b", 20, 30, "K002"),
            feature("c", 40, 50, "K003"),
        ]);
        let t = cluster(vec![
            feature("a", 1, 10, "K001"),
            feature("c", 40, 50, "K003"),
        ]);
        let scores = score_pair(&q, &t);
        assert!(scores.projected_order_score > 0.999);
    }

    #[test]
    fn inversion_is_flagged() {
        let q = cluster(vec![
            feature("a", 1, 10, "K001"),
            feature("b", 20, 30, "K002"),
            feature("c", 40, 50, "K003"),
        ]);
        let t = cluster(vec![
            feature("c", 1, 10, "K003"),
            feature("b", 20, 30, "K002"),
            feature("a", 40, 50, "K001"),
        ]);
        let scores = score_pair(&q, &t);
        assert!(scores.inversion_flag);
    }
}
