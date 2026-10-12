//! Stable search-format helpers re-exported for CLI and integration tests.

pub use crate::search::index::{collect_tgt_paths, read_genome_metadata, read_tgt_records};
pub use crate::search::query::{write_candidate_hits, write_placement_details, write_placements};
