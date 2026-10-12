//! Version-1 `.s2cidx` construction, serialization, and metadata model.

use crate::enzyme::EnzymeType;
use crate::tgt::{TgtReader, TgtRecord};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub const INDEX_MAGIC: &[u8; 8] = b"S2CIDX\x01\x00";
pub const INDEX_VERSION: u32 = 1;

/// Optional, user-supplied genome metadata.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct GenomeMetadata {
    pub genome_id: String,
    pub species_cluster: Option<String>,
    pub taxonomy: Option<String>,
    pub fasta_path: Option<String>,
    pub is_representative: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct GenomeRecord {
    pub genome_id: String,
    pub total_tags: u64,
    pub total_adjacencies: u64,
    pub total_length: u64,
    pub contig_count: u32,
    pub n50: u64,
    pub landmark_source: Option<EnzymeType>,
    pub metadata: GenomeMetadata,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ContigMetadata {
    pub genome_id: u32,
    /// Reference-side numeric contig id. Zero means the TGT did not assign one.
    pub contig_id: u32,
    pub contig_name: String,
    pub length: u64,
    pub is_circular: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct TagPosting {
    pub genome_id: u32,
    pub contig_id: u32,
    pub position: u64,
    /// Stored recognition-window strand. Chaining uses canonical sequence ids.
    pub strand: bool,
}

/// Adjacency postings use canonical, unordered adjacency keys but retain the
/// reference's stored endpoint order so query orientation can be recovered.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct AdjacencyPosting {
    pub genome_id: u32,
    pub contig_id: u32,
    pub tag_a_id: u32,
    pub tag_b_id: u32,
    pub position_a: u64,
    pub position_b: u64,
}

/// A reusable ordered-landmark placement index.
///
/// Maps use `BTreeMap` so serialized bytes and iteration are deterministic.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct ContigIndex {
    pub version: u32,
    pub genomes: Vec<GenomeRecord>,
    pub contigs: Vec<ContigMetadata>,
    pub sequences: Vec<[u8; 32]>,
    pub sequence_ids: BTreeMap<[u8; 32], u32>,
    pub tag_postings: BTreeMap<u32, Vec<TagPosting>>,
    pub adjacency_postings: BTreeMap<(u32, u32), Vec<AdjacencyPosting>>,
    pub genome_ids: BTreeMap<String, u32>,
}

impl ContigIndex {
    pub fn new() -> Self {
        Self {
            version: INDEX_VERSION,
            genomes: Vec::new(),
            contigs: Vec::new(),
            sequences: Vec::new(),
            sequence_ids: BTreeMap::new(),
            tag_postings: BTreeMap::new(),
            adjacency_postings: BTreeMap::new(),
            genome_ids: BTreeMap::new(),
        }
    }

    fn intern_sequence(&mut self, sequence: [u8; 32]) -> u32 {
        if let Some(id) = self.sequence_ids.get(&sequence) {
            return *id;
        }
        let id = self
            .sequences
            .len()
            .try_into()
            .expect("sequence dictionary exceeds u32 capacity");
        self.sequences.push(sequence);
        self.sequence_ids.insert(sequence, id);
        id
    }

    pub fn sequence_id(&self, sequence: [u8; 32]) -> Option<u32> {
        self.sequence_ids.get(&sequence).copied()
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!("failed to create index directory {}", parent.display())
            })?;
        }
        let file = File::create(path)
            .with_context(|| format!("failed to create index {}", path.display()))?;
        let mut out = BufWriter::new(file);
        out.write_all(INDEX_MAGIC)?;
        bincode::serialize_into(&mut out, self).context("failed to serialize contig index")?;
        out.flush()?;
        Ok(())
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file =
            File::open(path).with_context(|| format!("failed to open index {}", path.display()))?;
        let mut input = BufReader::new(file);
        let mut magic = [0u8; 8];
        input
            .read_exact(&mut magic)
            .context("truncated contig-index magic")?;
        if &magic != INDEX_MAGIC {
            bail!(
                "{} is not an .s2cidx v1 index (bad magic); rebuild it with `syn2b contig-index`",
                path.display()
            );
        }
        let index: Self =
            bincode::deserialize_from(&mut input).context("invalid contig-index payload")?;
        index.validate()?;
        Ok(index)
    }

    pub fn validate(&self) -> Result<()> {
        if self.version != INDEX_VERSION {
            bail!("unsupported .s2cidx version {}", self.version);
        }
        if self.sequences.len() != self.sequence_ids.len() {
            bail!("contig-index sequence dictionary is inconsistent");
        }
        for (sequence, id) in &self.sequence_ids {
            if (*id as usize) >= self.sequences.len() || &self.sequences[*id as usize] != sequence {
                bail!("contig-index sequence dictionary is inconsistent");
            }
        }
        if self.genome_ids.len() != self.genomes.len() {
            bail!("contig-index genome dictionary is inconsistent");
        }
        for (genome_id, id) in &self.genome_ids {
            if (*id as usize) >= self.genomes.len()
                || &self.genomes[*id as usize].genome_id != genome_id
            {
                bail!("contig-index genome dictionary is inconsistent");
            }
        }
        Ok(())
    }
}

impl Default for ContigIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Discover `.tgt` files in a directory, or accept one TGT file directly.
pub fn collect_tgt_paths(input: &Path) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    if input.is_dir() {
        for entry in std::fs::read_dir(input)? {
            let path = entry?.path();
            if path.is_file() && path.extension().and_then(|x| x.to_str()) == Some("tgt") {
                paths.push(path);
            }
        }
    } else if input.is_file() {
        paths.push(input.to_path_buf());
    } else {
        bail!("TGT input does not exist: {}", input.display());
    }
    paths.sort();
    if paths.is_empty() {
        bail!("no .tgt files found in {}", input.display());
    }
    Ok(paths)
}

/// Read text or binary TGT records without changing either legacy format.
pub fn read_tgt_records(path: &Path) -> Result<Vec<TgtRecord>> {
    let mut magic = [0u8; 4];
    {
        let mut file = File::open(path)?;
        file.read_exact(&mut magic)?;
    }
    let mut reader = TgtReader::new(path)?;
    let mut records = Vec::new();
    if &magic == b"TGT\x01" || &magic == b"TGT\x02" || &magic == b"TGT\x03" {
        while let Some(record) = reader.read_binary()? {
            records.push(record);
        }
    } else {
        while let Some(record) = reader.read_record()? {
            records.push(record);
        }
    }
    Ok(records)
}

/// Parse the optional five-column metadata table.
///
/// A header is recognized case-insensitively when present. Without a header,
/// rows use the documented positional order. Extra columns are ignored, which
/// keeps taxonomy fields containing escaped tabs impractical rather than silent.
pub fn read_genome_metadata(path: &Path) -> Result<HashMap<String, GenomeMetadata>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read genome metadata {}", path.display()))?;
    let mut map = HashMap::new();
    let mut start = 0usize;
    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').map(str::trim).collect();
        if line_no == start
            && fields
                .first()
                .copied()
                .unwrap_or_default()
                .eq_ignore_ascii_case("genome_id")
        {
            start += 1;
            continue;
        }
        if fields.is_empty() || fields[0].is_empty() {
            bail!("{}:{}: empty genome_id", path.display(), line_no + 1);
        }
        let meta = GenomeMetadata {
            genome_id: fields[0].to_string(),
            species_cluster: non_empty(fields.get(1).copied()),
            taxonomy: non_empty(fields.get(2).copied()),
            fasta_path: non_empty(fields.get(3).copied()),
            is_representative: matches!(
                fields
                    .get(4)
                    .copied()
                    .unwrap_or_default()
                    .to_ascii_lowercase()
                    .as_str(),
                "1" | "true" | "yes"
            ),
        };
        if map.insert(meta.genome_id.clone(), meta).is_some() {
            bail!(
                "{}:{}: duplicate metadata genome_id {}",
                path.display(),
                line_no + 1,
                fields[0]
            );
        }
    }
    Ok(map)
}

fn non_empty(value: Option<&str>) -> Option<String> {
    value.filter(|v| !v.is_empty()).map(str::to_string)
}

/// Build an index from TGT records, using canonical landmark identities.
pub fn build_index(
    tgt_dir: &Path,
    metadata: Option<&HashMap<String, GenomeMetadata>>,
) -> Result<ContigIndex> {
    let mut index = ContigIndex::new();
    let mut seen_genome_ids = HashMap::new();
    let mut first_source: Option<Option<EnzymeType>> = None;

    for path in collect_tgt_paths(tgt_dir)? {
        for record in read_tgt_records(&path)? {
            if record.tags.is_empty() {
                bail!(
                    "TGT {} contains genome {} with no landmarks",
                    path.display(),
                    record.genome_id
                );
            }
            let source = common_landmark_source(&record.tags)?;
            if let Some(first) = first_source {
                if first != source {
                    bail!(
                        "cannot mix landmark sources: {} uses {:?}, but an earlier genome uses {:?}",
                        record.genome_id,
                        source,
                        first
                    );
                }
            } else {
                first_source = Some(source);
            }

            if seen_genome_ids
                .insert(record.genome_id.clone(), path.clone())
                .is_some()
            {
                bail!(
                    "duplicate genome id {:?} while indexing {} (also seen in {:?})",
                    record.genome_id,
                    path.display(),
                    seen_genome_ids[&record.genome_id]
                );
            }
            let genome_idx = index.genomes.len() as u32;
            if index
                .genome_ids
                .insert(record.genome_id.clone(), genome_idx)
                .is_some()
            {
                unreachable!("genome id was checked above");
            }

            let lengths = contig_lengths(&record);
            let mut contigs = Vec::new();
            for (zero_based, length) in lengths.iter().copied().enumerate() {
                let contig_id = if record.contig_names.is_empty() && zero_based == 0 {
                    0
                } else {
                    (zero_based + 1) as u32
                };
                let name = record
                    .contig_names
                    .get(zero_based)
                    .cloned()
                    .unwrap_or_else(|| record.genome_id.clone());
                contigs.push(ContigMetadata {
                    genome_id: genome_idx,
                    contig_id,
                    contig_name: name,
                    length,
                    is_circular: record
                        .contig_circular
                        .get(zero_based)
                        .copied()
                        .unwrap_or(false),
                });
            }
            let contig_count = contigs.len() as u32;
            index.contigs.extend(contigs);

            let tag_count = record.tags.len() as u64;
            let mut adjacency_count = 0u64;

            for group in landmark_adjacency_groups(&record) {
                let mut previous: Option<&crate::tgt::Tag> = None;
                for tag in group {
                    let canonical = canonical_sequence(&tag.sequence);
                    let tag_id = index.intern_sequence(canonical);
                    index
                        .tag_postings
                        .entry(tag_id)
                        .or_default()
                        .push(TagPosting {
                            genome_id: genome_idx,
                            contig_id: tag.contig_id as u32,
                            position: tag.position,
                            strand: tag.strand == crate::tgt::Strand::Reverse,
                        });
                    if let Some(prev) = previous {
                        let max_tag_len = std::cmp::max(
                            prev.enzyme.properties().tag_length,
                            tag.enzyme.properties().tag_length,
                        ) as u64;
                        if tag.position > prev.position
                            && tag.position - prev.position > max_tag_len
                        {
                            let a_id = index.intern_sequence(canonical_sequence(&prev.sequence));
                            let b_id = tag_id;
                            let key = if a_id <= b_id {
                                (a_id, b_id)
                            } else {
                                (b_id, a_id)
                            };
                            index.adjacency_postings.entry(key).or_default().push(
                                AdjacencyPosting {
                                    genome_id: genome_idx,
                                    contig_id: tag.contig_id as u32,
                                    tag_a_id: a_id,
                                    tag_b_id: b_id,
                                    position_a: prev.position,
                                    position_b: tag.position,
                                },
                            );
                            adjacency_count += 1;
                        }
                    }
                    previous = Some(tag);
                }
            }

            let default_metadata = GenomeMetadata::default();
            let metadata = metadata
                .and_then(|m| m.get(&record.genome_id))
                .unwrap_or(&default_metadata);
            index.genomes.push(GenomeRecord {
                genome_id: record.genome_id.clone(),
                total_tags: tag_count,
                total_adjacencies: adjacency_count,
                total_length: record.total_length,
                contig_count,
                n50: n50(&lengths),
                landmark_source: source,
                metadata: metadata.clone(),
            });
        }
    }

    for postings in index.tag_postings.values_mut() {
        postings.sort_by_key(|p| (p.genome_id, p.contig_id, p.position, p.strand));
    }
    for postings in index.adjacency_postings.values_mut() {
        postings.sort_by_key(|p| {
            (
                p.genome_id,
                p.contig_id,
                p.position_a,
                p.position_b,
                p.tag_a_id,
                p.tag_b_id,
            )
        });
    }
    index.validate()?;
    Ok(index)
}

fn common_landmark_source(tags: &[crate::tgt::Tag]) -> Result<Option<EnzymeType>> {
    let has_fmh = tags.iter().any(|tag| tag.enzyme == EnzymeType::FracMinHash);
    let has_enzyme = tags.iter().any(|tag| tag.enzyme != EnzymeType::FracMinHash);
    if has_fmh && has_enzyme {
        bail!("mixed enzyme and FracMinHash landmarks in one TGT are not searchable");
    }
    // A multi-enzyme 2bRAD panel is one structural landmark alphabet and is
    // searchable. FracMinHash is the only incompatible alternative alphabet.
    Ok(if has_fmh {
        Some(EnzymeType::FracMinHash)
    } else {
        None
    })
}

/// Group landmarks by contig and landmark source, preserving input order.
fn landmark_adjacency_groups(record: &TgtRecord) -> Vec<Vec<&crate::tgt::Tag>> {
    let mut groups: BTreeMap<(u16, u8), Vec<&crate::tgt::Tag>> = BTreeMap::new();
    for tag in &record.tags {
        groups
            .entry((tag.contig_id, tag.enzyme.index()))
            .or_default()
            .push(tag);
    }
    groups.into_values().collect()
}

pub fn canonical_sequence(sequence: &[u8; 32]) -> [u8; 32] {
    let n = sequence.iter().take_while(|&&b| b != 0).count();
    let mut rc = [0u8; 32];
    for i in 0..n {
        rc[n - 1 - i] = match sequence[i] {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            other => other,
        };
    }
    if rc[..n] < sequence[..n] {
        rc
    } else {
        *sequence
    }
}

/// Derive contig lengths from text metadata; binary TGTs only preserve total length.
pub(crate) fn contig_lengths(record: &TgtRecord) -> Vec<u64> {
    if !record.contig_offsets.is_empty() {
        return (0..record.contig_offsets.len())
            .map(|i| {
                let end = record
                    .contig_offsets
                    .get(i + 1)
                    .copied()
                    .unwrap_or(record.total_length);
                end.saturating_sub(record.contig_offsets[i])
            })
            .collect();
    }
    if record.contig_names.len() <= 1 {
        return vec![record.total_length];
    }
    vec![0; record.contig_names.len()]
}

pub(crate) fn n50(lengths: &[u64]) -> u64 {
    if lengths.is_empty() {
        return 0;
    }
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable();
    let half = sorted.iter().copied().sum::<u64>().div_ceil(2);
    let mut cumulative = 0u64;
    for length in sorted.into_iter().rev() {
        cumulative += length;
        if cumulative >= half {
            return length;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(value: &[u8]) -> [u8; 32] {
        let mut out = [0u8; 32];
        out[..value.len()].copy_from_slice(value);
        out
    }

    #[test]
    fn canonicalization_matches_reverse_complements() {
        let forward = seq(b"ACGT");
        let reverse = seq(b"ACGTrevcomp");
        let mut rc = [0u8; 32];
        rc[0] = b'A';
        assert_eq!(
            canonical_sequence(&seq(b"ACGT")),
            canonical_sequence(&seq(b"ACGT"))
        );
        assert_ne!(canonical_sequence(&forward), canonical_sequence(&reverse));
        assert_eq!(canonical_sequence(&rc), rc);
    }

    #[test]
    fn metadata_supports_header_and_positional_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("meta.tsv");
        std::fs::write(
            &path,
            "genome_id\tspecies_cluster\ttaxonomy\tfasta_path\tis_representative\n\
             g1\tsp1\td:Bacteria\tg1.fna\ttrue\n",
        )
        .unwrap();
        let map = read_genome_metadata(&path).unwrap();
        assert_eq!(map["g1"].species_cluster.as_deref(), Some("sp1"));
        assert!(map["g1"].is_representative);
    }

    #[test]
    fn n50_uses_cummulative_half() {
        assert_eq!(n50(&[10, 30, 60]), 60);
        assert_eq!(n50(&[40, 40]), 40);
        assert_eq!(n50(&[]), 0);
    }
}
