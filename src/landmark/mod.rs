//! Landmark sources — where the ordered, identifiable points along a genome come from.
//!
//! Syn2b's structural mathematics never depended on restriction digestion. Every
//! metric in [`crate::synteny::scoring`] consumes only a list of
//! `(canonical identity, position, contig, orientation)`, so any rule that picks
//! reproducible loci can drive it. Two are available:
//!
//! - **2bRAD** ([`crate::enzyme::digest`]) — Type IIB restriction sites. Fixed
//!   anchors, biologically realisable, and the only mode that corresponds to an
//!   actual wet-lab protocol.
//! - **FracMinHash** ([`fracminhash`]) — `h(canonical(kmer)) < u64::MAX / scale`.
//!   Context-free and genome-independent, with a continuous density knob.
//!
//! The two differ in one way that the downstream code must know about, and it is
//! the reason this is a mode rather than a drop-in: **run collapse applies only to
//! the enzyme path.** Type IIB enzymes produce overlapping cut sites within
//! `MIN_TAG_SEPARATION`, so one physical locus can yield several tags and they must
//! be collapsed to a single representative. FracMinHash selects each position
//! independently, so two landmarks 20 bp apart are two genuine loci; collapsing
//! them would destroy real signal. See [`crate::synteny::scoring::collapse_runs`].

pub mod fracminhash;

pub use fracminhash::FracMinHash;
