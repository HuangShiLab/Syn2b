//! 2bSyn — 2bRAD-based Synteny Detection Engine
//!
//! A Rust-based alignment-free synteny detection engine that leverages
//! Type IIB restriction enzyme tags (2bRAD tags) from up to 16 enzymes
//! to detect genomic structural variations and infer synteny relationships
//! between microbial genomes.

pub mod enzyme;
pub mod io;
pub mod synteny;
pub mod tgt;
pub mod utils;
