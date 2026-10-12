//! `.s2cldb` v1: deterministic token/adjacency postings for gene clusters.

use crate::cluster::model::{read_cluster_file, GeneCluster};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub const CLUSTER_INDEX_MAGIC: &[u8; 8] = b"S2CLDB\x01\x00";
pub const CLUSTER_INDEX_VERSION: u32 = 1;

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct ClusterIndex {
    pub version: u32,
    pub clusters: Vec<GeneCluster>,
    pub cluster_ids: BTreeMap<String, u32>,
    pub token_postings: BTreeMap<String, BTreeSet<u32>>,
    pub adjacency_postings: BTreeMap<(String, String), BTreeSet<u32>>,
}

impl ClusterIndex {
    pub fn new() -> Self {
        Self {
            version: CLUSTER_INDEX_VERSION,
            clusters: Vec::new(),
            cluster_ids: BTreeMap::new(),
            token_postings: BTreeMap::new(),
            adjacency_postings: BTreeMap::new(),
        }
    }

    pub fn build(mut clusters: Vec<GeneCluster>) -> Result<Self> {
        clusters.sort_by(|a, b| a.cluster_id.cmp(&b.cluster_id));
        let mut index = Self::new();
        for cluster in clusters {
            if index.cluster_ids.contains_key(&cluster.cluster_id) {
                bail!("duplicate cluster id {}", cluster.cluster_id);
            }
            let id = index.clusters.len() as u32;
            for token in cluster.token_counts.keys() {
                index
                    .token_postings
                    .entry(token.clone())
                    .or_default()
                    .insert(id);
            }
            for (a, b) in &cluster.adjacencies {
                let key = if a <= b {
                    (a.clone(), b.clone())
                } else {
                    (b.clone(), a.clone())
                };
                index.adjacency_postings.entry(key).or_default().insert(id);
            }
            index.cluster_ids.insert(cluster.cluster_id.clone(), id);
            index.clusters.push(cluster);
        }
        index.validate()?;
        Ok(index)
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let mut out = BufWriter::new(File::create(path)?);
        out.write_all(CLUSTER_INDEX_MAGIC)?;
        bincode::serialize_into(&mut out, self).context("failed to serialize cluster index")?;
        out.flush()?;
        Ok(())
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let mut input = BufReader::new(File::open(path)?);
        let mut magic = [0u8; 8];
        input.read_exact(&mut magic)?;
        if &magic != CLUSTER_INDEX_MAGIC {
            bail!("{} is not .s2cldb v1", path.display());
        }
        let index: Self =
            bincode::deserialize_from(&mut input).context("invalid cluster index payload")?;
        index.validate()?;
        Ok(index)
    }

    pub fn validate(&self) -> Result<()> {
        if self.version != CLUSTER_INDEX_VERSION {
            bail!("unsupported .s2cldb version {}", self.version);
        }
        if self.clusters.len() != self.cluster_ids.len() {
            bail!("cluster index id map is inconsistent");
        }
        for (cluster_id, id) in &self.cluster_ids {
            if (*id as usize) >= self.clusters.len()
                || &self.clusters[*id as usize].cluster_id != cluster_id
            {
                bail!("cluster index id map is inconsistent");
            }
        }
        Ok(())
    }
}

impl Default for ClusterIndex {
    fn default() -> Self {
        Self::new()
    }
}

pub fn collect_cluster_files(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    for input in inputs {
        if input.is_dir() {
            for entry in std::fs::read_dir(input)? {
                let path = entry?.path();
                if path.extension().and_then(|value| value.to_str()) == Some("s2cl") {
                    paths.push(path);
                }
            }
        } else if input.is_file() {
            paths.push(input.clone());
        } else {
            bail!("cluster input does not exist: {}", input.display());
        }
    }
    paths.sort();
    paths.dedup();
    if paths.is_empty() {
        bail!("no .s2cl files found");
    }
    Ok(paths)
}

pub fn build_cluster_index(inputs: &[PathBuf]) -> Result<ClusterIndex> {
    let mut clusters = Vec::new();
    for path in collect_cluster_files(inputs)? {
        clusters.extend(
            read_cluster_file(&path)
                .with_context(|| format!("failed to read {}", path.display()))?,
        );
    }
    if clusters.is_empty() {
        bail!("no gene clusters found in input");
    }
    ClusterIndex::build(clusters)
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
            product: "p".into(),
            ko: ko.into(),
            pfam: String::new(),
            cog: String::new(),
            eggnog: String::new(),
            smash_tool: String::new(),
            cluster_id: Some("x".into()),
        }
    }

    #[test]
    fn cluster_index_roundtrip_and_postings() {
        let clusters = build_clusters(
            vec![
                feature("a", 1, 10, "K001"),
                feature("b", 20, 30, "K002"),
                feature("c", 40, 50, "K003"),
            ],
            1_000,
            2,
        )
        .unwrap();
        let index = ClusterIndex::build(clusters).unwrap();
        assert!(index.token_postings.contains_key("KO:K001"));
        assert_eq!(index.clusters.len(), 1);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("clusters.s2cldb");
        index.write(&path).unwrap();
        assert_eq!(ClusterIndex::read(&path).unwrap(), index);
    }

    #[test]
    fn duplicate_cluster_ids_are_rejected() {
        let clusters = build_clusters(
            vec![feature("a", 1, 10, "K001"), feature("b", 20, 30, "K002")],
            1_000,
            2,
        )
        .unwrap();
        let mut duplicate = clusters.clone();
        duplicate.push(clusters[0].clone());
        assert!(ClusterIndex::build(duplicate).is_err());
    }
}
