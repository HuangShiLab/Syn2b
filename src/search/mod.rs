//! Syn2b-Contig: reusable ordered-landmark search and placement.

pub mod chaining;
pub mod chimera;
pub mod index;
pub mod io;
pub mod query;
pub mod scaffold;
pub mod scoring;

pub use chimera::scan_record;
pub use chimera::{write_chimerism, ChimeraOptions, ChimeraSegment};
pub use index::{build_index, canonical_sequence};
pub use index::{
    AdjacencyPosting, ContigIndex, ContigMetadata, GenomeMetadata, TagPosting, INDEX_MAGIC,
    INDEX_VERSION,
};
pub use io::{
    collect_tgt_paths, read_genome_metadata, read_tgt_records, write_candidate_hits,
    write_placement_details, write_placements,
};
pub use query::{search_record, Placement, SearchOptions};
pub use query::{search_records, CandidateSummary};
pub use scaffold::run_scaffold;
pub use scoring::PlacementConfidence;
