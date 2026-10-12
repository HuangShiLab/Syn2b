//! Feature-table and minimal GFF ingestion for Syn2b-Cluster.
//!
//! Syn2b does not call genes or predict function. It consumes existing
//! annotations and normalizes them into the stable 13-column feature table.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

pub const FEATURE_COLUMNS: [&str; 13] = [
    "gene_id",
    "genome_id",
    "contig_id",
    "start",
    "end",
    "strand",
    "product",
    "ko",
    "pfam",
    "cog",
    "eggnog",
    "smash_tool",
    "cluster_id",
];

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Feature {
    pub gene_id: String,
    pub genome_id: String,
    pub contig_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub product: String,
    pub ko: String,
    pub pfam: String,
    pub cog: String,
    pub eggnog: String,
    pub smash_tool: String,
    pub cluster_id: Option<String>,
}

fn clean_optional(value: &str) -> String {
    value.trim().trim_end_matches(',').to_string()
}

fn normalize_strand(value: &str) -> char {
    match value.trim() {
        "+" | "1" => '+',
        "-" | "-1" => '-',
        _ => '+',
    }
}

pub fn feature_tokens(feature: &Feature) -> BTreeSet<String> {
    let mut tokens = BTreeSet::new();
    for (prefix, value) in [
        ("KO", &feature.ko),
        ("PFAM", &feature.pfam),
        ("COG", &feature.cog),
        ("EGGNOG", &feature.eggnog),
        ("SMASH", &feature.smash_tool),
    ] {
        for item in value.split(|c| c == ';' || c == ',' || c == '|') {
            let item = item.trim();
            if !item.is_empty() {
                tokens.insert(format!("{}:{}", prefix, item.to_ascii_uppercase()));
            }
        }
    }
    let product: String = feature
        .product
        .to_ascii_lowercase()
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|word| !word.is_empty() && !matches!(*word, "protein" | "hypothetical"))
        .take(4)
        .collect::<Vec<_>>()
        .join("_");
    if !product.is_empty() {
        tokens.insert(format!("PRODUCT:{product}"));
    }
    tokens
}

pub fn token_weight(token: &str) -> f64 {
    match token.split(':').next().unwrap_or_default() {
        "KO" | "EGGNOG" | "SMASH" => 3.0,
        "PFAM" | "COG" | "AMRCLASS" => 2.0,
        _ => 1.0,
    }
}

/// Parse the documented feature table. A header is required; extra columns are
/// rejected so a producer cannot silently shift fields.
pub fn parse_features(path: impl AsRef<Path>) -> Result<Vec<Feature>> {
    let path = path.as_ref();
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to open feature table {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    reader
        .read_line(&mut header)
        .context("failed to read feature-table header")?;
    let columns: Vec<String> = header
        .trim()
        .split('\t')
        .map(str::trim)
        .map(String::from)
        .collect();
    if columns != FEATURE_COLUMNS {
        bail!(
            "feature table header must be exactly: {}",
            FEATURE_COLUMNS.join("\t")
        );
    }

    let mut features = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line
            .with_context(|| format!("failed to read {} line {}", path.display(), line_no + 2))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != columns.len() {
            bail!(
                "{}:{}: expected {} columns, got {}",
                path.display(),
                line_no + 2,
                columns.len(),
                fields.len()
            );
        }
        let value = |name: &str| -> Result<String> {
            let index = columns
                .iter()
                .position(|column| column == name)
                .expect("column exists because header was checked");
            Ok(clean_optional(fields[index]))
        };
        let start = value("start")?
            .parse::<u64>()
            .with_context(|| format!("{}:{}: invalid start", path.display(), line_no + 2))?;
        let end = value("end")?
            .parse::<u64>()
            .with_context(|| format!("{}:{}: invalid end", path.display(), line_no + 2))?;
        if end < start {
            bail!("{}:{}: end < start", path.display(), line_no + 2);
        }
        features.push(Feature {
            gene_id: value("gene_id")?,
            genome_id: value("genome_id")?,
            contig_id: value("contig_id")?,
            start,
            end,
            strand: normalize_strand(&value("strand")?),
            product: value("product")?,
            ko: value("ko")?,
            pfam: value("pfam")?,
            cog: value("cog")?,
            eggnog: value("eggnog")?,
            smash_tool: value("smash_tool")?,
            cluster_id: {
                let value = value("cluster_id")?;
                if value.is_empty() {
                    None
                } else {
                    Some(value)
                }
            },
        });
    }
    Ok(features)
}

pub fn write_features(path: impl AsRef<Path>, features: &[Feature]) -> Result<()> {
    let mut out = std::io::BufWriter::new(std::fs::File::create(path.as_ref())?);
    writeln!(out, "{}", FEATURE_COLUMNS.join("\t"))?;
    for feature in features {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            feature.gene_id,
            feature.genome_id,
            feature.contig_id,
            feature.start,
            feature.end,
            feature.strand,
            feature.product,
            feature.ko,
            feature.pfam,
            feature.cog,
            feature.eggnog,
            feature.smash_tool,
            feature.cluster_id.as_deref().unwrap_or(""),
        )?;
    }
    Ok(())
}

/// Minimal GFF3/GTF parser. Functional evidence is accepted only from explicit
/// attributes; Syn2b never infers KO/Pfam/product roles.
pub fn parse_gff(path: impl AsRef<Path>, genome_id: &str) -> Result<Vec<Feature>> {
    let file = std::fs::File::open(path.as_ref())
        .with_context(|| format!("failed to open GFF {}", path.as_ref().display()))?;
    let mut features = Vec::new();
    for (line_no, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            bail!(
                "{}:{}: expected 9 GFF columns",
                path.as_ref().display(),
                line_no + 1
            );
        }
        if fields[2].eq_ignore_ascii_case("gene") || fields[2].eq_ignore_ascii_case("CDS") {
            let attributes = parse_gff_attributes(fields[8]);
            let gene_id = attributes
                .get("gene_id")
                .or_else(|| attributes.get("ID"))
                .or_else(|| attributes.get("Name"))
                .cloned()
                .unwrap_or_else(|| format!("{}_{}", genome_id, line_no + 1));
            let product = attributes
                .get("product")
                .or_else(|| attributes.get("Product"))
                .cloned()
                .unwrap_or_default();
            let ko = attributes
                .get("KO")
                .or_else(|| attributes.get("ko"))
                .or_else(|| attributes.get("KEGG_KO"))
                .cloned()
                .unwrap_or_default();
            let pfam = attributes
                .get("Pfam")
                .or_else(|| attributes.get("pfam"))
                .cloned()
                .unwrap_or_default();
            let cog = attributes
                .get("COG")
                .or_else(|| attributes.get("cog"))
                .cloned()
                .unwrap_or_default();
            let eggnog = attributes
                .get("eggNOG")
                .or_else(|| attributes.get("eggnog"))
                .cloned()
                .unwrap_or_default();
            let smash_tool = attributes
                .get("antiSMASH")
                .or_else(|| attributes.get("smash_tool"))
                .cloned()
                .unwrap_or_default();
            features.push(Feature {
                gene_id,
                genome_id: genome_id.to_string(),
                contig_id: fields[0].to_string(),
                start: fields[3].parse()?,
                end: fields[4].parse()?,
                strand: normalize_strand(fields[6]),
                product,
                ko,
                pfam,
                cog,
                eggnog,
                smash_tool,
                cluster_id: attributes.get("cluster_id").cloned(),
            });
        }
    }
    Ok(features)
}

fn parse_gff_attributes(text: &str) -> std::collections::BTreeMap<String, String> {
    let mut map = std::collections::BTreeMap::new();
    for item in text.split(';') {
        let item = item.trim();
        if item.is_empty() {
            continue;
        }
        if let Some((key, value)) = item.split_once('=') {
            map.insert(key.trim().to_string(), value.trim().to_string());
        }
    }
    map
}
