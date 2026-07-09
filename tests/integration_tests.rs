//! Integration tests for Syn2b.
//!
//! These tests verify end-to-end functionality across all modules.

use std::process::Command;

// ---------------------------------------------------------------------------
// Enzyme module tests
// ---------------------------------------------------------------------------

#[test]
fn test_enzyme_count() {
    // All 16 enzymes should be defined
    use bsyn::enzyme::EnzymeType;
    let all = EnzymeType::all();
    assert_eq!(all.len(), 16, "Expected 16 enzyme types");
}

#[test]
fn test_enzyme_type_variants() {
    use bsyn::enzyme::EnzymeType;
    // Verify all 16 variants exist
    let _ = EnzymeType::BcgI;
    let _ = EnzymeType::AlfI;
    let _ = EnzymeType::AloI;
    let _ = EnzymeType::BaeI;
    let _ = EnzymeType::BplI;
    let _ = EnzymeType::BsaXI;
    let _ = EnzymeType::BslFI;
    let _ = EnzymeType::Bsp24I;
    let _ = EnzymeType::CjeI;
    let _ = EnzymeType::CjePI;
    let _ = EnzymeType::CspCI;
    let _ = EnzymeType::FalI;
    let _ = EnzymeType::HaeIV;
    let _ = EnzymeType::Hin4I;
    let _ = EnzymeType::PpiI;
    let _ = EnzymeType::PsrI;
}

#[test]
fn test_enzyme_properties() {
    use bsyn::enzyme::EnzymeType;
    let bcgi = EnzymeType::BcgI.properties();
    assert_eq!(bcgi.enzyme_type, EnzymeType::BcgI);
    assert_eq!(bcgi.tag_length, 32);
}

// ---------------------------------------------------------------------------
// TGT module tests
// ---------------------------------------------------------------------------

#[test]
fn test_tag_creation() {
    use bsyn::enzyme::EnzymeType;
    use bsyn::tgt::{Strand, Tag};

    let seq = [b'A'; 32];
    let tag = Tag::new(seq, 100, EnzymeType::BcgI, Strand::Forward, 0);

    assert_eq!(tag.position, 100);
    assert_eq!(tag.enzyme, EnzymeType::BcgI);
    assert_eq!(tag.strand, Strand::Forward);
    assert_eq!(tag.sequence_str(), "A".repeat(32));
}

#[test]
fn test_hamming_distance() {
    use bsyn::enzyme::EnzymeType;
    use bsyn::tgt::{Strand, Tag};

    let seq1 = [b'A'; 32];
    let mut seq2 = [b'A'; 32];
    seq2[0] = b'T';
    seq2[1] = b'G';

    let tag1 = Tag::new(seq1, 0, EnzymeType::BcgI, Strand::Forward, 0);
    let tag2 = Tag::new(seq2, 0, EnzymeType::BcgI, Strand::Forward, 0);

    assert_eq!(tag1.hamming_distance(&tag2), 2);
}

#[test]
fn test_tgt_record_creation() {
    use bsyn::enzyme::EnzymeType;
    use bsyn::tgt::{Strand, Tag, TgtRecord};

    let mut record = TgtRecord::new("E_coli_K12", 4_641_652);
    assert_eq!(record.genome_id, "E_coli_K12");
    assert_eq!(record.total_length, 4_641_652);
    assert_eq!(record.tag_count(), 0);

    let tag = Tag::new([b'A'; 32], 1000, EnzymeType::BcgI, Strand::Forward, 0);
    record.add_tag(tag);

    assert_eq!(record.tag_count(), 1);
    assert_eq!(record.gaps.len(), 0); // First tag has no preceding gap
}

#[test]
fn test_tgt_record_gaps() {
    use bsyn::enzyme::EnzymeType;
    use bsyn::tgt::{Strand, Tag, TgtRecord};

    let mut record = TgtRecord::new("test_genome", 10000);

    let tag1 = Tag::new([b'A'; 32], 100, EnzymeType::BcgI, Strand::Forward, 0);
    let tag2 = Tag::new([b'C'; 32], 1500, EnzymeType::BcgI, Strand::Forward, 0);
    let tag3 = Tag::new([b'G'; 32], 3000, EnzymeType::AlfI, Strand::Reverse, 0);

    record.add_tag(tag1);
    record.add_tag(tag2);
    record.add_tag(tag3);

    assert_eq!(record.tag_count(), 3);
    assert_eq!(record.gaps.len(), 2);
    assert_eq!(record.gaps[0].size, 1400); // 1500 - 100
    assert_eq!(record.gaps[1].size, 1500); // 3000 - 1500
}

#[test]
fn test_gap_statistics() {
    use bsyn::enzyme::EnzymeType;
    use bsyn::tgt::{Strand, Tag, TgtRecord};

    let mut record = TgtRecord::new("test", 100000);

    // Add tags with known gaps: 1000, 2000, 3000
    record.add_tag(Tag::new([b'A'; 32], 0, EnzymeType::BcgI, Strand::Forward, 0));
    record.add_tag(Tag::new([b'C'; 32], 1000, EnzymeType::BcgI, Strand::Forward, 0));
    record.add_tag(Tag::new([b'G'; 32], 3000, EnzymeType::BcgI, Strand::Forward, 0));
    record.add_tag(Tag::new([b'T'; 32], 6000, EnzymeType::BcgI, Strand::Forward, 0));

    assert_eq!(record.mean_gap(), 2000.0); // (1000 + 2000 + 3000) / 3
    assert_eq!(record.max_gap(), 3000);
}

// ---------------------------------------------------------------------------
// Synteny module tests
// ---------------------------------------------------------------------------

#[test]
fn test_tag_adjacency_graph_creation() {
    use bsyn::synteny::TagAdjacencyGraph;

    let graph = TagAdjacencyGraph::new();
    assert_eq!(graph.node_count(), 0);
    assert_eq!(graph.edge_count(), 0);
    assert_eq!(graph.num_genomes, 0);
}

// ---------------------------------------------------------------------------
// IO module tests
// ---------------------------------------------------------------------------

#[test]
fn test_fasta_record_fields() {
    use bsyn::io::FastaRecord;

    let record = FastaRecord {
        id: "test_id".to_string(),
        sequence: b"ATGC".to_vec(),
        description: Some("description".to_string()),
    };
    assert_eq!(record.id, "test_id");
    assert_eq!(record.sequence, b"ATGC");
    assert_eq!(record.description, Some("description".to_string()));
    assert_eq!(record.sequence.len(), 4);
}

// ---------------------------------------------------------------------------
// Utility function tests
// ---------------------------------------------------------------------------

#[test]
fn test_reverse_complement() {
    use bsyn::utils::reverse_complement;

    assert_eq!(reverse_complement(b"ATGC"), vec![b'G', b'C', b'A', b'T']);
    assert_eq!(reverse_complement(b"AAAA"), vec![b'T', b'T', b'T', b'T']);
}

#[test]
fn test_gc_content() {
    use bsyn::utils::gc_content;

    assert_eq!(gc_content(b"ATGC"), 0.5);
    assert_eq!(gc_content(b"ATAT"), 0.0);
    assert_eq!(gc_content(b"GCGC"), 1.0);
}

// ---------------------------------------------------------------------------
// CLI tests
// ---------------------------------------------------------------------------

#[test]
fn test_cli_help() {
    // Verify the CLI binary can print help
    let output = Command::new("cargo")
        .args(["run", "--", "--help"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("Failed to execute cargo run -- --help");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    // Help output might be in stdout or stderr depending on how clap outputs
    let combined = format!("{}{}", stdout, stderr);
    
    assert!(
        combined.contains("2bRAD-based Synteny Detection Engine") || 
        combined.contains("Commands:") ||
        combined.contains("digest") ||
        combined.contains("synteny"),
        "CLI help should contain expected content. Got stdout: {}, stderr: {}",
        stdout, stderr
    );
}

// ---------------------------------------------------------------------------
// Round-trip test (binary format)
// ---------------------------------------------------------------------------

#[test]
fn test_binary_roundtrip() {
    use bsyn::enzyme::EnzymeType;
    use bsyn::tgt::{Strand, Tag, TgtRecord};
    use bsyn::tgt::{TgtReader, TgtWriter};
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    // Create a test record
    let mut record = TgtRecord::new("test_genome", 10000);
    record.add_tag(Tag::new([b'A'; 32], 100, EnzymeType::BcgI, Strand::Forward, 0));
    record.add_tag(Tag::new([b'C'; 32], 1500, EnzymeType::AlfI, Strand::Reverse, 0));
    record.add_tag(Tag::new([b'G'; 32], 3000, EnzymeType::BaeI, Strand::Forward, 0));

    // Write binary
    let tmpfile = NamedTempFile::new().unwrap();
    let path = tmpfile.path().to_path_buf();
    {
        let mut writer = TgtWriter::new(&path).unwrap();
        writer.write_binary(&record).unwrap();
        writer.flush().unwrap();
    }

    // Read binary back
    let mut reader = TgtReader::new(&path).unwrap();
    let read_record = reader.read_binary().unwrap();

    assert!(read_record.is_some());
    let read = read_record.unwrap();
    assert_eq!(read.genome_id, "test_genome");
    assert_eq!(read.total_length, 10000);
    assert_eq!(read.tag_count(), 3);
}

// ---------------------------------------------------------------------------
// Digestion tests
// ---------------------------------------------------------------------------

#[test]
fn test_bcg_i_digest_empty_sequence() {
    use bsyn::enzyme::digest::digest_genome;
    use bsyn::enzyme::EnzymeType;

    let seq = b"";
    let tags = digest_genome(seq, EnzymeType::BcgI);
    assert!(tags.is_empty());
}

#[test]
fn test_bcg_i_digest_no_sites() {
    use bsyn::enzyme::digest::digest_genome;
    use bsyn::enzyme::EnzymeType;

    // Sequence with no BcgI recognition site
    let seq = b"ATATATATATATATATATATATATATATATATATATATATATATATAT";
    let tags = digest_genome(seq, EnzymeType::BcgI);
    assert!(tags.is_empty());
}

#[test]
fn test_bcg_i_digest_real_site() {
    use bsyn::enzyme::digest::digest_genome;
    use bsyn::enzyme::EnzymeType;

    // BcgI: tag_length=32, forward anchors at offset 10=CGA, offset 19=TGC
    let seq = b"AAAAAAAAAACGAAAAAAATGCAAAAAAAAAA";
    assert_eq!(seq.len(), 32);
    let tags = digest_genome(seq, EnzymeType::BcgI);
    assert_eq!(tags.len(), 1);
    assert_eq!(tags[0].position, 0);
}
