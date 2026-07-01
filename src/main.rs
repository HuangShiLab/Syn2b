//! 2bSyn — 2bRAD-based Synteny Detection Engine
//!
//! CLI application for in silico digestion, TGT format conversion,
//! multi-enzyme coverage analysis, and synteny detection.

// The crate is intentionally named `Syn2b` (not snake_case); silence the lint.
#![allow(non_snake_case)]

use clap::{Parser, Subcommand};
use anyhow::Result;
use std::path::PathBuf;

use Syn2b::enzyme::EnzymeType;
use Syn2b::tgt::{Strand, Tag, TgtRecord};

#[derive(Parser)]
#[command(name = "Syn2b")]
#[command(about = "2bRAD-based Synteny Detection Engine")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// In silico digest genomes with 2bRAD enzymes
    Digest {
        /// Input FASTA file or directory
        #[arg(short, long)]
        input: PathBuf,
        /// Output file or directory
        #[arg(short, long)]
        output: PathBuf,
        /// Comma-separated enzyme list (default: BcgI). Use "all" for all 16 enzymes.
        #[arg(short, long, default_value = "BcgI")]
        enzymes: String,
        /// Output format: text or binary
        #[arg(short, long, default_value = "text")]
        format: String,
    },
    /// Compute synteny between genomes using TGT
    Synteny {
        /// Input TGT file or directory
        #[arg(short, long)]
        input: PathBuf,
        /// Output file or directory
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Analyze multi-enzyme coverage statistics
    Coverage {
        /// Input FASTA file or directory
        #[arg(short, long)]
        input: PathBuf,
        /// Comma-separated enzyme list (default: all)
        #[arg(short, long, default_value = "all")]
        enzymes: String,
    },
    /// Convert between TGT text and binary formats
    Convert {
        /// Input TGT file
        #[arg(short, long)]
        input: PathBuf,
        /// Output TGT file
        #[arg(short, long)]
        output: PathBuf,
        /// Output format: text or binary
        #[arg(short, long, default_value = "binary")]
        format: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Digest { input, output, enzymes, format } => {
            println!("Syn2b digest — In silico genome digestion");
            println!("  input:  {:?}", input);
            println!("  output: {:?}", output);
            println!("  enzymes: {}", enzymes);
            println!("  format: {}", format);

            // Parse enzyme list
            let enzyme_list = parse_enzyme_list(&enzymes);
            match enzyme_list {
                Ok(list) => {
                    println!("  Using {} enzyme(s): {:?}", list.len(), 
                             list.iter().map(|e| format!("{:?}", e)).collect::<Vec<_>>().join(", "));
                }
                Err(e) => {
                    eprintln!("Error parsing enzyme list: {}", e);
                    std::process::exit(1);
                }
            }

            // TODO: Full pipeline implementation:
            // 1. Read FASTA file(s) from input
            // 2. For each genome, call Syn2b::enzyme::digest_multi_enzyme
            // 3. Merge tags with Syn2b::enzyme::merge_multi_enzyme_tags
            // 4. Write TGT output (text or binary) using Syn2b::tgt::TgtWriter

            println!("Digest subcommand — full implementation in progress");
            println!("This will:");
            println!("  1. Read FASTA from {:?}", input);
            println!("  2. Digest with the specified enzyme(s)");
            println!("  3. Write TGT output to {:?} (format: {})", output, format);
            Ok(())
        }

        Commands::Synteny { input, output } => {
            println!("Syn2b synteny — Synteny detection");
            println!("  input:  {:?}", input);
            println!("  output: {:?}", output);

            // TODO: Full pipeline implementation:
            // 1. Read TGT file(s) from input
            // 2. Build TagAdjacencyGraph
            // 3. Find common tags and build edges
            // 4. Extract linear paths and synteny blocks
            // 5. Write synteny report

            println!("Synteny subcommand — full implementation in progress");
            Ok(())
        }

        Commands::Coverage { input, enzymes } => {
            println!("Syn2b coverage — Multi-enzyme coverage analysis");
            println!("  input:   {:?}", input);
            println!("  enzymes: {}", enzymes);

            // Parse enzyme list
            let enzyme_list = parse_enzyme_list(&enzymes);
            match enzyme_list {
                Ok(list) => {
                    println!("  Analyzing coverage for {} enzyme(s)", list.len());
                }
                Err(e) => {
                    eprintln!("Error parsing enzyme list: {}", e);
                    std::process::exit(1);
                }
            }

            // TODO: Full pipeline:
            // 1. Read FASTA from input
            // 2. For each enzyme, digest and count tags
            // 3. Compute combined coverage statistics
            // 4. Report per-enzyme and combined metrics

            println!("Coverage subcommand — full implementation in progress");
            Ok(())
        }

        Commands::Convert { input, output, format } => {
            println!("Syn2b convert — TGT format conversion");
            println!("  input:  {:?}", input);
            println!("  output: {:?}", output);
            println!("  format: {}", format);

            let is_binary = format == "binary";

            // TODO: Full pipeline:
            // 1. Open input TGT file (auto-detect text/binary)
            // 2. Read all records
            // 3. Write to output in requested format

            println!("Convert subcommand — full implementation in progress");
            println!("Will convert to {} format", if is_binary { "binary" } else { "text" });
            Ok(())
        }
    }
}

/// Parse a comma-separated enzyme list string into a Vec<EnzymeType>.
/// Special value "all" returns all 16 enzymes.
fn parse_enzyme_list(s: &str) -> Result<Vec<EnzymeType>> {
    if s == "all" {
        return Ok(EnzymeType::all().to_vec());
    }

    let mut enzymes = Vec::new();
    for name in s.split(',') {
        let name = name.trim();
        let enzyme = parse_enzyme_name(name)?;
        enzymes.push(enzyme);
    }

    if enzymes.is_empty() {
        anyhow::bail!("No valid enzymes specified");
    }

    Ok(enzymes)
}

/// Parse a single enzyme name into an EnzymeType.
fn parse_enzyme_name(name: &str) -> Result<EnzymeType> {
    match name {
        "BcgI" => Ok(EnzymeType::BcgI),
        "AlfI" => Ok(EnzymeType::AlfI),
        "AloI" => Ok(EnzymeType::AloI),
        "BaeI" => Ok(EnzymeType::BaeI),
        "BplI" => Ok(EnzymeType::BplI),
        "BsaXI" => Ok(EnzymeType::BsaXI),
        "BslFI" => Ok(EnzymeType::BslFI),
        "Bsp24I" => Ok(EnzymeType::Bsp24I),
        "CjeI" => Ok(EnzymeType::CjeI),
        "CjePI" => Ok(EnzymeType::CjePI),
        "CspCI" => Ok(EnzymeType::CspCI),
        "FalI" => Ok(EnzymeType::FalI),
        "HaeIV" => Ok(EnzymeType::HaeIV),
        "Hin4I" => Ok(EnzymeType::Hin4I),
        "PpiI" => Ok(EnzymeType::PpiI),
        "PsrI" => Ok(EnzymeType::PsrI),
        _ => anyhow::bail!("Unknown enzyme: {}", name),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_enzyme_list_all() {
        let list = parse_enzyme_list("all").unwrap();
        assert_eq!(list.len(), 16);
    }

    #[test]
    fn test_parse_enzyme_list_single() {
        let list = parse_enzyme_list("BcgI").unwrap();
        assert_eq!(list.len(), 1);
        assert_eq!(list[0], EnzymeType::BcgI);
    }

    #[test]
    fn test_parse_enzyme_list_multiple() {
        let list = parse_enzyme_list("BcgI,AlfI,AloI").unwrap();
        assert_eq!(list.len(), 3);
        assert_eq!(list[0], EnzymeType::BcgI);
        assert_eq!(list[1], EnzymeType::AlfI);
        assert_eq!(list[2], EnzymeType::AloI);
    }

    #[test]
    fn test_parse_enzyme_name_invalid() {
        assert!(parse_enzyme_name("NotAnEnzyme").is_err());
    }
}
