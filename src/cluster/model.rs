//! Gene-cluster representation and deterministic `.s2cl` I/O.

use crate::cluster::features::{feature_tokens, Feature};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

pub const CLUSTER_FILE_MAGIC: &[u8; 6] = b"S2CL\x01\x00";

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct GeneCluster {
    pub cluster_id: String,
    pub genome_id: String,
    pub contig_id: String,
    pub start: u64,
    pub end: u64,
    pub genes: Vec<Feature>,
    /// Each entry is the token set for one ordered gene.
    pub gene_tokens: Vec<BTreeMap<String, u32>>,
    pub token_counts: BTreeMap<String, u32>,
    /// Canonical token adjacency. Direction is retained as an ordered pair;
    /// scoring currently treats pair presence as orientation-free.
    pub adjacencies: BTreeSet<(String, String)>,
}

use std::collections::BTreeSet;

impl GeneCluster {
    fn new(cluster_id: String, mut genes: Vec<Feature>) -> Self {
        genes.sort_by_key(|gene| (gene.start, gene.end, gene.gene_id.clone()));
        let start = genes.first().map(|gene| gene.start).unwrap_or_default();
        let end = genes.iter().map(|gene| gene.end).max().unwrap_or_default();
        let contig_id = genes
            .first()
            .map(|gene| gene.contig_id.clone())
            .unwrap_or_default();
        let genome_id = genes
            .first()
            .map(|gene| gene.genome_id.clone())
            .unwrap_or_default();
        let mut gene_tokens = Vec::new();
        let mut token_counts = BTreeMap::new();
        let mut ordered_tokens = Vec::new();
        for gene in &genes {
            let mut counts = BTreeMap::new();
            for token in feature_tokens(gene) {
                *counts.entry(token.clone()).or_default() += 1;
                *token_counts.entry(token.clone()).or_default() += 1;
            }
            // One representative per gene makes projected-order comparison stable.
            let representative = counts.keys().next().cloned().unwrap_or_default();
            if !representative.is_empty() {
                ordered_tokens.push(representative);
            }
            gene_tokens.push(counts);
        }
        let mut adjacencies = BTreeSet::new();
        for pair in ordered_tokens.windows(2) {
            adjacencies.insert((pair[0].clone(), pair[1].clone()));
        }
        Self {
            cluster_id,
            genome_id,
            contig_id,
            start,
            end,
            genes,
            gene_tokens,
            token_counts,
            adjacencies,
        }
    }
}

/// Build gene clusters. Explicit `cluster_id` always wins. Otherwise genes are
/// chained only when they share genome/contig and the gap is at most `max_gap`.
pub fn build_clusters(
    mut features: Vec<Feature>,
    max_gap: u64,
    min_genes: usize,
) -> Result<Vec<GeneCluster>> {
    if max_gap == 0 {
        bail!("--max-gap must be > 0");
    }
    if min_genes == 0 {
        bail!("--min-genes must be > 0");
    }
    features.sort_by(|a, b| {
        a.genome_id
            .cmp(&b.genome_id)
            .then(a.contig_id.cmp(&b.contig_id))
            .then(a.cluster_id.cmp(&b.cluster_id))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.gene_id.cmp(&b.gene_id))
    });

    let mut explicit: BTreeMap<(String, String, String, String), Vec<Feature>> = BTreeMap::new();
    let mut implicit: Vec<Feature> = Vec::new();
    for feature in features {
        match &feature.cluster_id {
            Some(id) => explicit
                .entry((
                    feature.genome_id.clone(),
                    feature.contig_id.clone(),
                    id.clone(),
                    String::new(),
                ))
                .or_default()
                .push(feature),
            None => implicit.push(feature),
        }
    }

    let mut groups = explicit;
    let mut current: Vec<Feature> = Vec::new();
    for feature in implicit {
        let contiguous = current.last().is_some_and(|previous| {
            previous.genome_id == feature.genome_id
                && previous.contig_id == feature.contig_id
                && feature.start.saturating_sub(previous.end) <= max_gap
        });
        if contiguous {
            current.push(feature);
        } else {
            if current.len() >= min_genes {
                groups
                    .entry((
                        current[0].genome_id.clone(),
                        current[0].contig_id.clone(),
                        String::new(),
                        String::new(),
                    ))
                    .or_default()
                    .append(&mut current);
            }
            current.clear();
            current.push(feature);
        }
    }
    if current.len() >= min_genes {
        groups
            .entry((
                current[0].genome_id.clone(),
                current[0].contig_id.clone(),
                String::new(),
                String::new(),
            ))
            .or_default()
            .append(&mut current);
    }

    let mut clusters = Vec::new();
    for ((genome_id, contig_id, explicit_id, _), mut genes) in groups {
        if genes.len() < min_genes {
            continue;
        }
        let cluster_id = if explicit_id.is_empty() {
            format!(
                "{genome_id}:{contig_id}:{:09}-{:09}",
                genes.first().unwrap().start,
                genes.last().unwrap().end
            )
        } else {
            format!("{genome_id}:{explicit_id}")
        };
        clusters.push(GeneCluster::new(cluster_id, std::mem::take(&mut genes)));
    }
    clusters.sort_by(|a, b| a.cluster_id.cmp(&b.cluster_id));
    Ok(clusters)
}

pub fn write_cluster_file(path: impl AsRef<Path>, clusters: &[GeneCluster]) -> Result<()> {
    let path = path.as_ref();
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut out = BufWriter::new(File::create(path)?);
    out.write_all(CLUSTER_FILE_MAGIC)?;
    bincode::serialize_into(&mut out, &clusters).context("failed to serialize .s2cl clusters")?;
    out.flush()?;
    Ok(())
}

pub fn read_cluster_file(path: impl AsRef<Path>) -> Result<Vec<GeneCluster>> {
    let path = path.as_ref();
    let mut input = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 6];
    input.read_exact(&mut magic)?;
    if &magic != CLUSTER_FILE_MAGIC {
        bail!("{} is not a v1 .s2cl file", path.display());
    }
    let clusters = bincode::deserialize_from(&mut input).context("invalid .s2cl payload")?;
    Ok(clusters)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::features::Feature;

    fn gene(id: &str, start: u64, end: u64, ko: &str, cluster: Option<&str>) -> Feature {
        Feature {
            gene_id: id.into(),
            genome_id: "g".into(),
            contig_id: "c".into(),
            start,
            end,
            strand: '+',
            product: "transporter".into(),
            ko: ko.into(),
            pfam: String::new(),
            cog: String::new(),
            eggnog: String::new(),
            smash_tool: String::new(),
            cluster_id: cluster.map(str::to_string),
        }
    }

    #[test]
    fn explicit_cluster_id_overrides_gap_grouping() {
        let clusters = build_clusters(
            vec![
                gene("a", 1, 100, "K001", Some("bgc")),
                gene("b", 900_000, 900_100, "K002", Some("bgc")),
            ],
            5_000,
            2,
        )
        .unwrap();
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].genes.len(), 2);
        assert_eq!(clusters[0].cluster_id, "g:bgc");
    }

    #[test]
    fn implicit_clusters_respect_gap_and_minimum() {
        let features = vec![
            gene("a", 1, 100, "K001", None),
            gene("b", 2_000, 2_100, "K002", None),
            gene("far", 100_000, 100_100, "K003", None),
        ];
        let clusters = build_clusters(features, 5_000, 2).unwrap();
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].genes.len(), 2);
    }

    #[test]
    fn cluster_file_roundtrip_is_stable() {
        let clusters = build_clusters(
            vec![
                gene("a", 1, 100, "K001", Some("x")),
                gene("b", 500, 600, "K002", Some("x")),
            ],
            1_000,
            2,
        )
        .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("genome.s2cl");
        write_cluster_file(&path, &clusters).unwrap();
        let back = read_cluster_file(&path).unwrap();
        assert_eq!(back, clusters);
    }
}
