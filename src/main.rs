//! Syn2b — 2bRAD-based Synteny Detection Engine

use anyhow::Result;
use bsyn::enzyme::EnzymeType;
use bsyn::tgt::{Tag, TgtRecord};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "syn2b",
    about = "2bRAD-based Synteny Detection Engine",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract landmarks from a genome into a TGT.
    ///
    /// `--mode 2brad` (default) digests with Type IIB restriction enzymes.
    /// `--mode fracminhash` selects k-mers by `h(canonical) < u64::MAX / scale`
    /// instead; `--enzymes` is then unused and `--kmer` / `--scale` apply.
    Digest {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, default_value = "BcgI")]
        enzymes: String,
        #[arg(short, long, default_value = "text")]
        format: String,
        /// Landmark source: `2brad` or `fracminhash` (aliases: `enzyme`, `fmh`).
        #[arg(long, default_value = "2brad")]
        mode: String,
        /// FracMinHash k-mer length, 1..=32. Ignored in 2brad mode.
        #[arg(long, default_value_t = 31)]
        kmer: usize,
        /// FracMinHash compression: one k-mer in `scale` is kept, so expected
        /// landmark spacing is `scale` bp. Ignored in 2brad mode.
        #[arg(long, default_value_t = 1000)]
        scale: u64,
    },
    Synteny {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
    },
    Coverage {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, default_value = "all")]
        enzymes: String,
    },
    Convert {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, default_value = "binary")]
        format: String,
    },
    Scaffold {
        #[arg(short, long)]
        reference: PathBuf,
        #[arg(short, long)]
        draft: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, default_value = "3")]
        min_tags: usize,
    },
    /// Build a reusable ordered-landmark contig-placement index (.s2cidx).
    ContigIndex {
        #[arg(long)]
        genome_tgt_dir: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        genome_metadata: Option<PathBuf>,
    },
    /// Search query contigs against a contig-placement index.
    SearchContig {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "5")]
        min_tags: usize,
        #[arg(long, default_value = "3")]
        min_adjacencies: usize,
        #[arg(long, default_value_t = 25)]
        max_targets: usize,
    },
    /// Search query contigs and write chained-block placement details.
    PlaceContig {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        details: PathBuf,
        #[arg(long, default_value = "5")]
        min_tags: usize,
        #[arg(long, default_value = "3")]
        min_adjacencies: usize,
        #[arg(long, default_value_t = 25)]
        max_targets: usize,
    },
    /// Segment query contigs and report dominant genome/cluster transitions.
    ScanChimera {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value = "50kb")]
        window: String,
        #[arg(long, default_value_t = 10)]
        min_window_tags: usize,
        #[arg(long, default_value_t = 25)]
        max_targets: usize,
        #[arg(long, default_value_t = 0.20)]
        min_transition_confidence: f64,
    },

    /// Normalize an existing GFF3/GTF or feature annotation table.
    ClusterIngest {
        #[arg(long)]
        genome_id: Option<String>,
        #[arg(long, conflicts_with = "features")]
        gff: Option<PathBuf>,
        #[arg(long, conflicts_with = "gff")]
        features: Option<PathBuf>,
        #[arg(long)]
        output: PathBuf,
    },
    /// Build gene clusters from a normalized feature table (.s2cl).
    ClusterBuild {
        #[arg(long)]
        features: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 5_000)]
        max_gap: u64,
        #[arg(long, default_value_t = 3)]
        min_genes: usize,
    },
    /// Build a reusable gene-cluster database (.s2cldb).
    ClusterIndex {
        #[arg(long, num_args = 1..)]
        cluster_files: Vec<PathBuf>,
        #[arg(long)]
        output: PathBuf,
    },
    /// Search query gene clusters against a gene-cluster database.
    SearchCluster {
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 5_000)]
        max_gap: u64,
        #[arg(long, default_value_t = 3)]
        min_genes: usize,
        #[arg(long, default_value_t = 0.25)]
        min_functional_jaccard: f64,
        #[arg(long, default_value_t = 0.10)]
        min_synteny: f64,
        #[arg(long, default_value_t = 25)]
        max_targets: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Digest {
            input,
            output,
            enzymes,
            format,
            mode,
            kmer,
            scale,
        } => {
            let src = parse_landmark_mode(&mode, kmer, scale)?;
            // `--enzymes` is only parsed in 2brad mode, so an unknown name is
            // still an error there but a leftover default is not one here.
            let ez = match src {
                LandmarkMode::TwoBRad => parse_enzyme_list(&enzymes)?,
                LandmarkMode::FracMinHash(_) => Vec::new(),
            };
            run_digest(&input, &output, &ez, &format, src)?;
            Ok(())
        }
        Commands::Synteny { input, output } => {
            run_synteny(&input, &output)?;
            Ok(())
        }
        Commands::Coverage { input: _, enzymes } => {
            println!(
                "Coverage — WIP ({} enzymes)",
                parse_enzyme_list(&enzymes)?.len()
            );
            Ok(())
        }
        Commands::Convert { .. } => {
            println!("Convert — WIP");
            Ok(())
        }
        Commands::Scaffold {
            reference,
            draft,
            output,
            min_tags,
        } => bsyn::search::run_scaffold(&reference, &draft, &output, min_tags),
        Commands::ContigIndex {
            genome_tgt_dir,
            output,
            genome_metadata,
        } => {
            let metadata = genome_metadata
                .as_deref()
                .map(bsyn::search::read_genome_metadata)
                .transpose()?;
            let index = bsyn::search::build_index(&genome_tgt_dir, metadata.as_ref())?;
            index.write(&output)?;
            println!(
                "Indexed {} genome(s), {} contig(s), {} landmark(s) -> {:?}",
                index.genomes.len(),
                index.contigs.len(),
                index.sequences.len(),
                output
            );
            Ok(())
        }
        Commands::SearchContig {
            index,
            query,
            output,
            min_tags,
            min_adjacencies,
            max_targets,
        } => {
            let index = bsyn::search::ContigIndex::read(&index)?;
            let records = bsyn::search::read_tgt_records(&query)?;
            let options = bsyn::search::SearchOptions {
                min_tags,
                min_adjacencies,
                max_targets,
            };
            let placements = bsyn::search::search_records(&records, &index, options);
            bsyn::search::write_placements(&output, &placements)?;
            println!(
                "Placed {} query contig(s) -> {:?}",
                placements.len(),
                output
            );
            Ok(())
        }
        Commands::PlaceContig {
            index,
            query,
            output,
            details,
            min_tags,
            min_adjacencies,
            max_targets,
        } => {
            let index = bsyn::search::ContigIndex::read(&index)?;
            let records = bsyn::search::read_tgt_records(&query)?;
            let options = bsyn::search::SearchOptions {
                min_tags,
                min_adjacencies,
                max_targets,
            };
            let placements = bsyn::search::search_records(&records, &index, options);
            bsyn::search::write_placements(&output, &placements)?;
            bsyn::search::write_placement_details(&details, &placements)?;
            let hits = output.with_extension("hits.tsv");
            bsyn::search::write_candidate_hits(&hits, &placements)?;
            println!(
                "Placed {} query contig(s) -> {:?}",
                placements.len(),
                output
            );
            Ok(())
        }
        Commands::ScanChimera {
            index,
            query,
            output,
            window,
            min_window_tags,
            max_targets,
            min_transition_confidence,
        } => {
            let index = bsyn::search::ContigIndex::read(&index)?;
            let records = bsyn::search::read_tgt_records(&query)?;
            let options = bsyn::search::ChimeraOptions {
                window: parse_size(&window)?,
                min_window_tags,
                max_targets,
                min_transition_confidence,
            };
            let mut segments = Vec::new();
            for record in &records {
                segments.extend(bsyn::search::scan_record(record, &index, options));
            }
            bsyn::search::write_chimerism(&output, &segments)?;
            println!(
                "Wrote {} chimerism window(s) -> {:?}",
                segments.len(),
                output
            );
            Ok(())
        }
        Commands::ClusterIngest {
            genome_id,
            gff,
            features,
            output,
        } => {
            let mut parsed = if let Some(path) = gff {
                let genome_id = genome_id
                    .as_deref()
                    .ok_or_else(|| anyhow::anyhow!("--genome-id is required with --gff"))?;
                bsyn::cluster::parse_gff(&path, genome_id)?
            } else if let Some(path) = features {
                let mut parsed = bsyn::cluster::parse_features(&path)?;
                if let Some(genome_id) = genome_id {
                    for feature in &mut parsed {
                        feature.genome_id = genome_id.clone();
                    }
                }
                parsed
            } else {
                anyhow::bail!("provide --gff or --features");
            };
            parsed.sort_by(|a, b| {
                a.genome_id
                    .cmp(&b.genome_id)
                    .then(a.contig_id.cmp(&b.contig_id))
                    .then(a.start.cmp(&b.start))
                    .then(a.gene_id.cmp(&b.gene_id))
            });
            bsyn::cluster::write_features(&output, &parsed)?;
            println!("Ingested {} feature(s) -> {:?}", parsed.len(), output);
            Ok(())
        }
        Commands::ClusterBuild {
            features,
            output,
            max_gap,
            min_genes,
        } => {
            let clusters = bsyn::cluster::build_clusters(
                bsyn::cluster::parse_features(&features)?,
                max_gap,
                min_genes,
            )?;
            bsyn::cluster::write_cluster_file(&output, &clusters)?;
            println!("Built {} gene cluster(s) -> {:?}", clusters.len(), output);
            Ok(())
        }
        Commands::ClusterIndex {
            cluster_files,
            output,
        } => {
            let index = bsyn::cluster::build_cluster_index(&cluster_files)?;
            index.write(&output)?;
            println!(
                "Indexed {} gene cluster(s), {} token posting(s) -> {:?}",
                index.clusters.len(),
                index.token_postings.len(),
                output
            );
            Ok(())
        }
        Commands::SearchCluster {
            query,
            index,
            output,
            max_gap,
            min_genes,
            min_functional_jaccard,
            min_synteny,
            max_targets,
        } => {
            let index = bsyn::cluster::ClusterIndex::read(&index)?;
            let options = bsyn::cluster::ClusterSearchOptions {
                min_functional_jaccard,
                min_synteny,
                max_targets,
            };
            let hits =
                bsyn::cluster::search_feature_table(&query, &index, max_gap, min_genes, options)?;
            bsyn::cluster::write_cluster_hits(&output, &hits)?;
            println!("Wrote {} cluster hit(s) -> {:?}", hits.len(), output);
            Ok(())
        }
    }
}

// ── Synteny ─────────────────────────────────────────────────────────────────

fn run_synteny(input_dir: &PathBuf, output_path: &PathBuf) -> Result<()> {
    use bsyn::synteny::graph::TagAdjacencyGraph;
    use bsyn::synteny::scoring::{pairwise_synteny_matrix, structural_synteny};
    use bsyn::tgt::TgtReader;
    use std::io::Write;

    // Collect all .tgt files from input directory
    let mut tgt_files = Vec::new();
    if input_dir.is_dir() {
        for entry in std::fs::read_dir(input_dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.extension().map_or(false, |ext| ext == "tgt") {
                tgt_files.push(path);
            }
        }
    } else if input_dir.extension().map_or(false, |ext| ext == "tgt") {
        tgt_files.push(input_dir.clone());
    }

    if tgt_files.is_empty() {
        anyhow::bail!("No .tgt files found in {:?}", input_dir);
    }

    // Sort deterministically so callers can control which genome is genome_A
    // (the reference for raw_inverted_fraction) by naming the file first
    // alphabetically.
    tgt_files.sort();

    println!("Found {} TGT file(s)", tgt_files.len());

    // Build graph
    let mut graph = TagAdjacencyGraph::new();
    let mut genome_ids = Vec::new();
    // Records are retained so the structural metric can work on tag order
    // directly. The graph collapses genomes that share an id, which silently
    // yields a meaningless score, so duplicates are rejected instead.
    let mut records = Vec::new();
    // (genome_id, is_fracminhash) per record, to reject a mixed-source run.
    let mut sources: Vec<(String, bool)> = Vec::new();

    for path in &tgt_files {
        let mut reader = TgtReader::new(path)?;
        while let Some(record) = reader.read_record()? {
            let genome_id = record.genome_id.clone();
            if genome_ids.contains(&genome_id) {
                anyhow::bail!(
                    "duplicate genome id {:?}: two genomes sharing an id collapse into \
                     one node set and produce a meaningless score. Rename one.",
                    genome_id
                );
            }
            // Enzyme tags and FracMinHash k-mers are different alphabets of
            // landmark. Mixing them in one run is not a comparison that yields a
            // small number — it yields an empty intersection, which reads as
            // "unrelated genomes" and is the wrong diagnosis. Caught here, where
            // the file names are still in hand to name in the message.
            let fmh = !record.tags.is_empty()
                && record
                    .tags
                    .iter()
                    .all(|t| t.enzyme == EnzymeType::FracMinHash);
            if let Some((first_id, first_fmh)) = sources.first() {
                if *first_fmh != fmh {
                    let name = |f: bool| if f { "FracMinHash" } else { "2bRAD enzyme" };
                    anyhow::bail!(
                        "landmark sources differ: {:?} is {}, {:?} is {}. \
                         The two select different kinds of landmark and share none, \
                         so a comparison across them would report no similarity \
                         rather than an error. Re-run `digest` on both with the same \
                         --mode.",
                        first_id,
                        name(*first_fmh),
                        genome_id,
                        name(fmh)
                    );
                }
            }
            sources.push((genome_id.clone(), fmh));
            println!(
                "  Adding genome: {} ({} {})",
                genome_id,
                record.tags.len(),
                if fmh { "FracMinHash landmarks" } else { "tags" }
            );
            graph.add_genome(&genome_id, &record);
            genome_ids.push(genome_id);
            records.push(record);
        }
    }

    graph.build_edges();

    println!(
        "Graph: {} nodes, {} edges, {} genomes",
        graph.node_count(),
        graph.edge_count(),
        graph.num_genomes
    );

    // Compute pairwise matrix
    let matrix = pairwise_synteny_matrix(&graph);

    // Write output
    let mut f = std::fs::File::create(output_path)?;
    writeln!(f, "# Syn2b Pairwise Synteny Matrix")?;
    writeln!(
        f,
        "# breakpoints = junctions: adjacencies of A that B does not"
    )?;
    writeln!(
        f,
        "#   conserve. Exactly 2 per inversion, 3 per translocation, 0"
    )?;
    writeln!(
        f,
        "#   under substitution loads to at least 5%. This is the column"
    )?;
    writeln!(f, "#   to use; junction coordinates are written alongside.")?;
    writeln!(
        f,
        "# scj_distance = symmetric difference of the adjacency sets, the"
    )?;
    writeln!(
        f,
        "#   single-cut-or-join distance; twice breakpoints when circular."
    )?;
    writeln!(
        f,
        "# inverted_fraction = shared landmarks whose orientation differs,"
    )?;
    writeln!(
        f,
        "#   over the orientation-informative ones. Measures how much of"
    )?;
    writeln!(
        f,
        "#   the genome moved; breakpoints measures how often. A fraction,"
    )?;
    writeln!(
        f,
        "#   so it is comparable with alignment-based synteny measures."
    )?;
    writeln!(
        f,
        "# raw_inverted_fraction = same numerator, but relative to genome_A"
    )?;
    writeln!(
        f,
        "#   (no majority-frame flip). Ranges in [0,1] and matches fixed-"
    )?;
    writeln!(f, "#   reference alignment methods such as dnadiff.")?;
    writeln!(
        f,
        "# observable_fraction = share of A's adjacencies that B can judge"
    )?;
    writeln!(
        f,
        "#   at all; the rest are stranded at contig ends. THIS is where"
    )?;
    writeln!(
        f,
        "#   contig count belongs -- as a discount on detection power, not"
    )?;
    writeln!(
        f,
        "#   as a divisor of breakpoints. Expect to recover about this"
    )?;
    writeln!(f, "#   fraction of the true junctions.")?;
    writeln!(
        f,
        "# structural = conserved adjacencies over their union. Kept for"
    )?;
    writeln!(
        f,
        "#   continuity, but it divides a discrete event by the landmark"
    )?;
    writeln!(f, "#   count, so one inversion moves it by ~0.001.")?;
    writeln!(
        f,
        "# legacy_adjacency = the original metric: all consecutive tags, no"
    )?;
    writeln!(
        f,
        "#   canonicalisation. Dominated by substitution load, not structure."
    )?;
    writeln!(
        f,
        "genome_A,genome_B,breakpoints,scj_distance,breakpoint_density,\
inverted_fraction,raw_inverted_fraction,orientation_mismatches,orientation_mismatches_raw,orientation_uninformative,\
observable_fraction,observable_adjacencies,structural,\
shared_tags,repeats_dropped,landmarks_collapsed,circular,legacy_adjacency"
    )?;

    // Junction coordinates go to their own file: they are the deliverable, and
    // a variable-length list does not belong in a matrix cell.
    let junction_path = output_path.with_extension("junctions.tsv");
    let mut jf = std::fs::File::create(&junction_path)?;
    writeln!(jf, "genome_A\tgenome_B\tjunction_pos_in_A")?;

    for i in 0..genome_ids.len() {
        for j in (i + 1)..genome_ids.len() {
            let gi = &genome_ids[i];
            let gj = &genome_ids[j];
            let legacy = matrix.get(&(gi.clone(), gj.clone())).unwrap_or(&0.0);
            match structural_synteny(&records[i], &records[j]) {
                Some(r) => {
                    writeln!(
                        f,
                        "{},{},{},{},{:.5},{:.5},{:.5},{},{},{},{:.4},{},{:.4},{},{},{},{},{:.4}",
                        gi,
                        gj,
                        r.breakpoints,
                        r.scj_distance,
                        r.breakpoint_density,
                        r.inverted_fraction,
                        r.raw_inverted_fraction,
                        r.orientation_mismatches,
                        r.orientation_mismatches_raw,
                        r.orientation_uninformative,
                        r.observable_fraction,
                        r.observable_adjacencies,
                        r.score,
                        r.shared_tags,
                        r.repeats_dropped,
                        r.landmarks_collapsed,
                        r.circular,
                        legacy
                    )?;
                    for pos in &r.junctions {
                        writeln!(jf, "{}\t{}\t{}", gi, gj, pos)?;
                    }
                }
                None => writeln!(
                    f,
                    "{},{},NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,0,0,false,{:.4}",
                    gi, gj, legacy
                )?,
            }
        }
    }

    println!("Synteny matrix written to {:?}", output_path);
    println!("Junction coordinates written to {:?}", junction_path);
    Ok(())
}

// ── Digest ──────────────────────────────────────────────────────────────────

fn run_digest(
    inp: &PathBuf,
    out: &PathBuf,
    enzymes: &[EnzymeType],
    fmt: &str,
    src: LandmarkMode,
) -> Result<()> {
    use bsyn::io::fasta::FastaReader;
    use bsyn::tgt::TgtWriter;
    let bin = fmt == "binary";
    let mut tw = TgtWriter::new(out)?;
    let mut r = FastaReader::new(inp)?;
    let mut cid = 1u16;
    let mut off = 0u64;
    let mut names = Vec::new();
    let mut circular = Vec::new();
    let mut all_tags: Vec<Tag> = Vec::new();
    let mut gid = String::new();
    let mut tot = 0u64;
    while let Some(rec) = r.next_record()? {
        if gid.is_empty() {
            gid = rec.id.clone();
        }
        let sl = rec.sequence.len() as u64;
        tot += sl;
        names.push(rec.id.clone());
        circular.push(header_says_circular(&rec.id, rec.description.as_deref()));
        match src {
            LandmarkMode::TwoBRad => {
                use rayon::prelude::*;
                let enzyme_tags: Vec<Vec<Tag>> = enzymes
                    .par_iter()
                    .map(|&e| {
                        bsyn::enzyme::digest::digest_genome_contig(&rec.sequence, e, cid, off)
                    })
                    .collect();
                for t in enzyme_tags {
                    all_tags.extend(t);
                }
            }
            // Not parallelised over anything: the rolling encoder is a single pass
            // that must carry state across the contig, and contigs are already the
            // unit of the outer loop.
            LandmarkMode::FracMinHash(f) => {
                all_tags.extend(f.landmarks(&rec.sequence, cid, off));
            }
        }
        off += sl;
        cid += 1;
    }
    let mut rec = TgtRecord::new(&gid, tot);
    for t in all_tags {
        rec.add_tag(t);
    }
    rec.contig_names = names;
    // All-false means the headers said nothing, which must stay distinguishable
    // from "every contig is linear" — closure falls back to the single-contig
    // rule only when the topology is unknown.
    if circular.iter().any(|&c| c) {
        rec.contig_circular = circular;
    }
    let mut r2 = FastaReader::new(inp)?;
    let mut o = 0u64;
    while let Some(fr) = r2.next_record()? {
        rec.contig_offsets.push(o);
        o += fr.sequence.len() as u64;
    }
    if bin {
        tw.write_binary(&rec)?;
    } else {
        tw.write_record(&rec)?;
    }
    let n_circ = rec.contig_circular.iter().filter(|&&c| c).count();
    let what = match src {
        LandmarkMode::TwoBRad => format!("{} enzyme(s)", enzymes.len()),
        LandmarkMode::FracMinHash(f) => format!("FracMinHash k={} scale={}", f.k, f.scale),
    };
    println!(
        "Digested {} contigs ({} circular), {} bp, {} landmarks [{}] -> {:?}",
        rec.contig_names.len(),
        n_circ,
        tot,
        rec.tags.len(),
        what,
        out
    );
    Ok(())
}

/// Whether a FASTA header declares the sequence circular.
///
/// Long-read assemblers say so in the header, but each in its own words: Flye and
/// Unicycler write `circular=true`, Canu writes `suggestCircular=yes`, and the
/// INSDC/NCBI convention is `[topology=circular]`. Some pipelines only append
/// `_circular` to the name. A header that says nothing leaves the topology
/// unknown, which is not the same as linear.
fn header_says_circular(id: &str, desc: Option<&str>) -> bool {
    let hay = format!("{} {}", id, desc.unwrap_or("")).to_ascii_lowercase();
    for key in [
        "circular=true",
        "circular=yes",
        "suggestcircular=yes",
        "topology=circular",
        "completeness=complete circular",
    ] {
        if hay.contains(key) {
            return true;
        }
    }
    // bare word forms: `[circular]`, trailing `_circular`, ` circular` as a token
    hay.split(|c: char| !c.is_ascii_alphanumeric())
        .any(|t| t == "circular")
        && !hay.contains("circular=false")
        && !hay.contains("circular=no")
}

// ── Helpers ─────────────────────────────────────────────────────────────────

/// Which rule picks the landmarks.
///
/// The two are not interchangeable inside one TGT: `structural_synteny` refuses a
/// comparison across sources, because enzyme tags and FracMinHash k-mers are
/// different alphabets of landmark and an empty intersection between them would
/// otherwise read as "unrelated genomes".
#[derive(Clone, Copy, Debug)]
enum LandmarkMode {
    TwoBRad,
    FracMinHash(bsyn::landmark::FracMinHash),
}

fn parse_landmark_mode(mode: &str, k: usize, scale: u64) -> Result<LandmarkMode> {
    match mode.to_ascii_lowercase().as_str() {
        "2brad" | "enzyme" | "digest" => Ok(LandmarkMode::TwoBRad),
        "fracminhash" | "fmh" | "sketch" => Ok(LandmarkMode::FracMinHash(
            bsyn::landmark::FracMinHash::new(k, scale)?,
        )),
        other => anyhow::bail!(
            "Unknown --mode {:?} (expected `2brad` or `fracminhash`)",
            other
        ),
    }
}

fn parse_enzyme_list(s: &str) -> Result<Vec<EnzymeType>> {
    if s == "all" {
        return Ok(EnzymeType::all().to_vec());
    }
    let mut v = Vec::new();
    for n in s.split(',') {
        v.push(parse_enzyme_name(n.trim())?);
    }
    if v.is_empty() {
        anyhow::bail!("No enzymes");
    }
    Ok(v)
}
/// Parse sizes accepted by the chimera CLI: `50000`, `50kb`, `1mb`, `2kbp`.
fn parse_size(value: &str) -> Result<u64> {
    let value = value.trim().to_ascii_lowercase();
    let (digits, multiplier) = if let Some(x) = value.strip_suffix("kbp") {
        (x, 1_024u64.pow(1))
    } else if let Some(x) = value.strip_suffix("kb") {
        (x, 1_024)
    } else if let Some(x) = value.strip_suffix("mbp") {
        (x, 1_024u64.pow(2))
    } else if let Some(x) = value.strip_suffix("mb") {
        (x, 1_024u64.pow(2))
    } else {
        (value.as_str(), 1)
    };
    let parsed = digits
        .parse::<u64>()
        .map_err(|_| anyhow::anyhow!("invalid size {:?}; expected e.g. 50000 or 50kb", value))?;
    parsed
        .checked_mul(multiplier)
        .ok_or_else(|| anyhow::anyhow!("size {:?} overflows u64", value))
}

fn parse_enzyme_name(n: &str) -> Result<EnzymeType> {
    match n {
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
        "FracMinHash" | "FMH" => anyhow::bail!(
            "FracMinHash is a landmark source, not an enzyme — use `--mode fracminhash`"
        ),
        _ => anyhow::bail!("Unknown enzyme: {}", n),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn t1() {
        assert_eq!(parse_enzyme_list("all").unwrap().len(), 16);
    }
    #[test]
    fn t2() {
        assert_eq!(parse_enzyme_list("BcgI").unwrap()[0], EnzymeType::BcgI);
    }
    #[test]
    fn t3() {
        assert!(parse_enzyme_name("X").is_err());
    }
}
