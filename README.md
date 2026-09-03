# Syn2b — 2bRAD-based Synteny Detection Engine

**Syn2b** is a Rust, alignment-free engine for detecting genome synteny and
structural variation from **2bRAD tags** — the short, fixed-position sequences
produced by *Type IIB* restriction enzymes. Instead of aligning whole genomes,
Syn2b represents each genome as an ordered series of sparse anchor tags (the
**Tag–Gap–Tag / TGT** model) and infers synteny from how the *order and
adjacency* of those tags are conserved across genomes.

> **Status: active research tool.** The core pipeline (digest → TGT → graph →
> synteny/scaffold) is implemented, unit-tested, and benchmarked. The engine now
> supports both 2bRAD enzyme landmarks and FracMinHash landmarks, and has been
> validated on synthetic rearrangements, GTDB-scale pairs, and published isolate
> cohorts. See [Project status & roadmap](#project-status--roadmap).

---

## Table of contents

- [Motivation](#motivation)
- [How it works](#how-it-works)
- [Features](#features)
- [Repository layout](#repository-layout)
- [The TGT format](#the-tgt-format)
- [Landmark sources — 2bRAD or FracMinHash](#landmark-sources--2brad-or-fracminhash)
  - [Full comparison: `docs/LANDMARK_COMPARISON.md`](docs/LANDMARK_COMPARISON.md)
- [Type IIB restriction enzymes](#type-iib-restriction-enzymes)
- [Building](#building)
- [Command-line usage](#command-line-usage)
- [Library usage](#library-usage)
- [Synteny algorithm](#synteny-algorithm)
- [Scoring metrics](#scoring-metrics)
- [Project status & roadmap](#project-status--roadmap)
- [Validation & limitations](#validation--limitations)
- [Testing](#testing)
- [References](#references)
- [License](#license)

---

## Motivation

2bRAD-M and related reduced-representation methods profile a microbial genome
with roughly **1,000–3,000 species-specific tags** (27–33 bp), each anchored at a
fixed, restriction-map-determined coordinate. These tags are essentially a very
sparse, deterministic *minimizer set*: their **relative order is fixed by the
genome**, so a rearrangement (inversion, translocation, large indel) shows up as
a change in which tags are adjacent.

Existing synteny tools cannot consume this data directly. SynTracker, for
example, needs ~5 kb homologous regions with flanking context for BLAST/DECIPHER
alignment — a 27 bp isolated tag carries neither the length nor the context it
requires. Syn2b is a **dedicated method** that works natively on sparse tag
series, extending 2bRAD-M from purely *taxonomic* profiling toward *strain-level
structural* comparison.

Design ideas borrow from **ntSynt** (minimizer-graph synteny) and
**KmerAperture** (ordered k-mer series for structural variation), treating each
tag as a graph node and each observed adjacency as an edge.

---

## How it works

```
   FASTA genome(s)
        │
        │  in silico digestion  (enzyme/digest.rs)
        ▼
   2bRAD tags  ──────────────►  TGT record  (one per genome per enzyme)
        │                        Tag–Gap–Tag: ordered tags + inter-tag gaps
        │                        + contig_names / contig_offsets / contig_id
        │  text or binary I/O   (tgt/reader.rs, tgt/writer.rs)
        ▼
   ┌─────────────────┐
   │  Synteny mode   │          (synteny/graph.rs, synteny/scoring.rs)
   │  TagAdjacencyGraph           • nodes = unique tag sequences
   │    • build_edges()           • edges = adjacencies, weighted
   │    • simplify()    ─────►  Pairwise matrix (Jaccard similarity)
 │    • linear_paths() ─────►  Synteny backbones
   └─────────────────┘
        │
   ┌─────────────────┐
   │ Scaffold mode   │          (main.rs scaffold subcommand)
   │  Map draft contigs         • Evaluate fwd/rev orientation
   │  onto reference            • Sort by median ref position
   │  via shared tags   ────►  AGP v2.1 output
   └─────────────────┘
```

---

## Features

- **In silico digestion** with 16 Type IIB restriction enzymes (IUPAC-degenerate
  recognition sites supported: BaeI, HaeIV, Hin4I).
- **TGT (Tag–Gap–Tag) representation** — a compact, structured genome model of
  ordered tags plus the base-pair gaps between them. TGT v2 includes contig
  metadata (`contig_id`, `contig_names`, `contig_offsets`).
- **Two on-disk formats** — a human-readable text format and a compact,
  fixed-layout binary format, with round-trip conversion.
- **Streaming FASTA reader** — genomes are parsed record-by-record without
  loading the whole file into memory.
- **Tag adjacency graph** — ntSynt-inspired graph construction, low-support
  simplification, and linear-path (backbone) extraction.
- **Synteny quantification** — path synteny score, all-vs-all pairwise matrix,
  adjacency Jaccard, Kendall's τ on tag order, and breakpoint counts.
- **Synteny blocks & indels** — block extraction with per-genome coordinates,
  size filtering, and indel detection from inter-tag distance differences.
- **Scaffold subcommand** — map draft genome contigs onto a reference genome
  using shared 2bRAD tags, with orientation detection and AGP v2.1 output.

---

## Repository layout

```
syn2b/
├── Cargo.toml
├── src/
│   ├── main.rs                # CLI entry point (clap): digest (--mode 2brad|fracminhash) / synteny / scaffold / coverage / convert
│   ├── lib.rs                 # public API re-exports (crate name: `bsyn`)
│   ├── tgt/                   # TGT core data structures
│   │   ├── tag.rs             # Tag  (32 bp sequence + position + enzyme + strand + contig_id)
│   │   ├── gap.rs             # Gap  (inter-tag distance in bp)
│   │   ├── record.rs          # TgtRecord (one genome: ordered tags + gaps + contig metadata)
│   │   ├── writer.rs          # TgtWriter (text + binary output, v2 format)
│   │   └── reader.rs          # TgtReader (text + binary input, v2 format)
│   ├── enzyme/                # Restriction-enzyme definitions
│   │   ├── enzyme.rs          # EnzymeType (16 enzymes + the FracMinHash marker) + IUPAC matching
│   │   └── digest.rs          # in silico digestion (pattern matching + tag extraction)
│   ├── landmark/              # Landmark sources — what the structural metrics consume
│   │   └── fracminhash.rs     # FracMinHash selection: h(canonical(kmer)) < u64::MAX / scale
│   ├── synteny/               # Synteny-detection engine
│   │   ├── graph.rs           # TagAdjacencyGraph, TagNode, AdjacencyEdge
│   │   ├── scoring.rs         # synteny_score, pairwise matrix, Jaccard, Kendall τ, breakpoints
│   │   └── blocks.rs          # SyntenyBlock extraction, indel detection, size filtering
│   ├── io/
│   │   └── fasta.rs           # streaming FASTA reader (FastaReader / FastaRecord)
│   └── utils/
│       └── mod.rs             # reverse_complement, is_valid_dna, gc_content
└── tests/
    └── integration_tests.rs   # end-to-end tests (enzyme, TGT, graph, I/O, binary round-trip)
```

The library crate is named **`bsyn`** and the binary is **`syn2b`**.

---

## The TGT format

A **TGT record** describes one genome as an ordered list of tags with the gaps
between them. `gaps[i]` is the distance (bp) between `tags[i]` and `tags[i+1]`,
so a record with *N* tags has *N − 1* gaps.

TGT v2 adds **multi-contig support**: each tag carries a `contig_id` (u16),
and the record stores `contig_names` and `contig_offsets` for mapping back to
FASTA headers.

### Text format

```
#contigs=NC_000913:4641652;NC_000914:5000000
>genome_id|length=9641652
BcgI:ATCG… -1313- GCTA… -1298- TTAA…
AlfI:CGAT… -892-  AATT… -1567- GCAT…
```

- A header comment `#contigs=name:length;...` lists contig names and lengths.
- A header line `>genome_id|length=<bp>` opens each record.
- Tag lines are grouped by enzyme; each tag is written as `Enzyme:SEQUENCE`.
- Gaps appear between consecutive tags as `-<size>-`.

The reader reconstructs tag positions from the cumulative gaps and validates that
each parsed gap matches the value recomputed from tag positions.

### Binary format (v2)

A fixed-layout, little-endian binary encoding:

```
Header (48 bytes)
  [0..4]    Magic           b"TGT\x02"
  [4..8]    Version         u32 (= 2)
  [8..16]   Genome length   u64
  [16..20]  Tag count       u32
  [20..22]  Enzyme count    u16
  [22..24]  Contig count    u16
  [24..48]  Reserved

Tag table (tag_count × 48 bytes)
  [0..32]   Tag sequence    raw bytes (zero-padded)
  [32..40]  Position        u64
  [40..41]  Enzyme index    u8 (0–15)
  [41..42]  Strand          u8 (0 = forward, 1 = reverse)
  [42..44]  Contig ID       u16
  [44..48]  Reserved

Gap table ((tag_count − 1) × 4 bytes)
  [0..4]    Gap size        u32

Contig name table (variable)
  For each contig: u16 name_len + name bytes
```

On read, the magic bytes and version are verified and the stored gap table is
cross-checked against gaps recomputed from tag positions.

---

## Landmark sources — 2bRAD or FracMinHash

> Measured head-to-head comparison, with the mechanism behind each difference:
> **[docs/LANDMARK_COMPARISON.md](docs/LANDMARK_COMPARISON.md)**.

Syn2b's structural mathematics never depended on restriction digestion. Every metric
in `synteny::scoring` consumes only a list of `(canonical identity, position, contig,
orientation)`, so any rule that picks reproducible loci can drive it. Two are
available, selected with `digest --mode`:

| | `--mode 2brad` (default) | `--mode fracminhash` |
|---|---|---|
| selection | Type IIB restriction sites | `h(canonical(kmer)) < u64::MAX / scale` |
| density | discrete: 1, 2, 4 or 16 enzymes | continuous, via `--scale` |
| wet-lab realisable | yes | no — in silico only |
| run collapse | required | must not be applied |

`--scale` is expected landmark spacing in bp, so `--scale 750` on *E. coli* K-12
gives 6,034 landmarks against the four-enzyme panel's 6,079 — a matched-density
comparison. Measured across `--scale` 250 to 2000, observed density tracks
`4.54 Mb / scale` to within 1.7%.

**Why FracMinHash and not minimizers.** Minimizer selection is window-relative, so a
single substitution re-selects a whole neighbourhood and landmark identity stops
being stable across genomes — which an adjacency-based structural metric cannot
tolerate. FracMinHash is *context-free*: whether a k-mer is kept depends on the k-mer
alone. It is also genome-independent, unlike bottom-*s* MinHash, whose cutoff is the
*s*-th smallest hash of that particular genome.

**The one place the modes diverge, and why it is a mode rather than a drop-in.** Type
IIB enzymes cut on both sides of a site they recognise once, so a single physical
locus can yield several tags within `MIN_TAG_SEPARATION` and they must be collapsed
to one representative. FracMinHash evaluates each position independently, so two
selected k-mers 20 bp apart are two genuine loci; collapsing them would delete real
landmarks — symmetrically, and therefore silently. `collapse_runs` is gated on the
source, which is read from the tags themselves rather than from a record-level flag
that could disagree with them. Comparing a 2bRAD TGT against a FracMinHash TGT is
refused with a message, not scored as zero similarity.

### Verified equivalence, E. coli K-12

Both sources at ~6,000 landmarks. Every structural control gives the same answer:

| control | 4-enzyme panel | FracMinHash k=31 s=750 | truth |
|---|---|---|---|
| self-comparison | 0 bp, 0 SCJ, obs 1.0000 | 0 bp, 0 SCJ, obs 1.0000 | collinear |
| origin rotation 1.2 Mb | 0 bp, 0 SCJ | 0 bp, 0 SCJ | no change |
| 500 kb inversion | 2 bp, 4 SCJ, f 0.10919 | 2 bp, 4 SCJ, f 0.11224 | 2, 4, f 0.1101 |
| 1 / 2 / 3 / 5 inversions | 2/4, 4/8, 6/12, 10/20 | 2/4, 4/8, 6/12, 10/20 | 2R, 4R |
| 40-contig shatter | 0 bp, obs 0.9930 | 0 bp, obs 0.9934 | 0 bp, K = 40 |

### Where they differ: near-duplicate landmarks

Substitution ladder, no rearrangement, so every junction is false:

| substitution | enzyme bp | enzyme SCJ | enzyme kept | FMH bp | FMH SCJ | FMH kept |
|---|---|---|---|---|---|---|
| 0.1% | 0 | 0 | 89.5% | 0 | 0 | 94.6% |
| 1% | 0 | 0 | 67.7% | 0 | 0 | 71.1% |
| 3% | 0 | **6** | 36.4% | 0 | **0** | 37.6% |
| 5% | 0 | **18** | 19.1% | 0 | **0** | 19.0% |

`breakpoints` is 0 for both, because the >=2-landmark relocation rule
(`docs/MATH_REVIEW.md`) rejects the paralog-convergence artifact. But `scj_distance`
is the *unfiltered* symmetric difference, and there the enzyme path still carries 6
and 18 while FracMinHash carries none. The mechanism needs a unique landmark sitting one substitution away from a
**multi-copy family**: the family is dropped by the per-genome uniqueness filter,
but in a diverged genome, once enough of its copies are destroyed the survivor
becomes unique and collides with the other locus. Measured on E. coli K-12, counting
unique landmarks within one substitution of a multi-copy family in either
orientation:

| source | unique landmarks | multi-copy families | at risk | share |
|---|---|---|---|---|
| BcgI | 2,809 | 13 | 7 | 0.249% |
| four-enzyme panel | 5,889 | 38 | 20 | 0.340% |
| FracMinHash k=31 s=1582 | 2,776 | 17 | **0** | **0.000%** |
| FracMinHash k=31 s=750 | 5,880 | 39 | **0** | **0.000%** |

FracMinHash carries just as many genuine multi-copy families — repeats are a property
of the genome, not of the selection rule — but none of its unique landmarks sits one
substitution from one. Enzyme landmarks must contain a recognition motif, so they are
crammed into a far smaller region of sequence space and near-collisions are
correspondingly likelier; FracMinHash k-mers are drawn from the whole space with no
shared constraint.

So the relocation rule exists to protect the enzyme path from a failure mode the
sketch path does not have.

## Type IIB restriction enzymes

Type IIB enzymes cut on **both** sides of their recognition site, excising a
short, defined fragment — ideal for producing fixed-length 2bRAD tags. Syn2b
ships definitions for 16 of them (`src/enzyme/enzyme.rs`). Recognition sites use
anchor motifs plus IUPAC constraints (`Y` = C/T; `R` = A/G; `[GAC]` = G/A/C).

| Enzyme  | Tag length | IUPAC degenerate | Notes |
|---------|:----------:|:----------------:|-------|
| BcgI    | 32         | —                | Two anchors: CGA @10, TGC @19 |
| AlfI    | 32         | —                | Palindrome: GCA/TGC |
| AloI    | 27         | —                | |
| BaeI    | 28         | Y (fwd), R (rev) | `[CT]` @19 fwd, `[AG]` @8 rev |
| BplI    | 27         | —                | Palindrome: GAG/CTC |
| BsaXI   | 27         | —                | |
| BslFI   | 25         | —                | |
| Bsp24I  | 27         | —                | |
| CjeI    | 28         | —                | |
| CjePI   | 27         | —                | |
| CspCI   | 33         | —                | |
| FalI    | 27         | —                | Palindrome: AAG/CTT |
| HaeIV   | 27         | Y + R            | Y@9+R@15 fwd; Y@11+R@17 rev |
| Hin4I   | 27         | Y + [GAC]        | Y@10+[GAC]@16 fwd; [CTG]@10+R@16 rev |
| PpiI    | 27         | —                | |
| PsrI    | 27         | —                | |

Tag lengths are derived from the 2bRAD-M reference implementation and cross-
validated against Fast2bRAD-M. The pattern-matching engine supports both exact
anchor motifs and IUPAC bitmask constraints (bit0=A, bit1=T, bit2=C, bit3=G).

---

## Building

Requires a stable Rust toolchain (Rust 2021 edition; install via
[rustup](https://rustup.rs)).

```bash
cd syn2b
cargo build --release
```

Dependencies (from `Cargo.toml`): `clap` (CLI), `anyhow` (errors), `serde` +
`bincode` (serialization), `parking_lot`, and `rayon` (parallelism);
`tempfile` for tests.

---

## Command-line usage

```
syn2b <COMMAND> [OPTIONS]

Commands:
  digest     In silico digest genomes with 2bRAD enzymes
  synteny    Compute synteny between genomes using TGT
  scaffold   Map draft contigs onto a reference using shared tags
  coverage   Analyze multi-enzyme coverage statistics
  convert    Convert between TGT text and binary formats
```

Examples:

```bash
# Digest a genome with the default enzyme (BcgI) into a text TGT file
syn2b digest -i genome.fasta -o genome.tgt

# Digest with all 16 enzymes, writing the compact binary format
syn2b digest -i genome.fasta -o genome.btgt --enzymes all --format binary

# Select landmarks by FracMinHash instead of by digestion. Same TGT, same
# downstream analysis; --enzymes is unused and --kmer/--scale apply.
syn2b digest -i genome.fasta -o genome.tgt --mode fracminhash --kmer 31 --scale 1000

# Compute synteny across a set of TGT records
syn2b synteny -i tgts/ -o synteny_matrix.csv

# Scaffold: map draft contigs onto a reference genome
syn2b scaffold -r reference.tgt -d draft.tgt -o scaffolds.agp --min-tags 3

# Multi-enzyme coverage statistics
syn2b coverage -i genome.fasta --enzymes all

# Convert a text TGT to binary (or vice versa)
syn2b convert -i genome.tgt -o genome.btgt --format binary
```

Common options: `-i/--input`, `-o/--output`, `-e/--enzymes`
(comma-separated names or `all`), `-f/--format` (`text` | `binary`).

### Scaffold subcommand

The `scaffold` command maps draft contigs onto a reference genome:

1. Load reference and draft TGT records.
2. For each draft contig, evaluate forward and reverse-complement orientations
   by matching tags against the reference.
3. Use a count-ratio heuristic (>2× difference) to pick the dominant
   orientation.
4. Sort contigs by median reference position.
5. Output AGP v2.1 with real contig lengths and estimated gap sizes.

This has been validated on E. coli K-12 self-scaffold (4 reversed contigs
correctly identified) and ABHQ draft (135 contigs → 45 anchored at
`min_tags=3`).

---

## Library usage

The crate exposes its functionality as the `bsyn` library. Representative API:

```rust
use bsyn::enzyme::{Enzyme, EnzymeType};
use bsyn::enzyme::digest::{digest_genome, digest_genome_contig};
use bsyn::tgt::{Tag, TgtRecord, TgtReader, TgtWriter, Strand};
use bsyn::synteny::{TagAdjacencyGraph, extract_synteny_blocks, synteny_score};
use std::path::Path;

// 1. Enzyme properties
let bcgi = Enzyme::properties(EnzymeType::BcgI);
assert_eq!(bcgi.tag_length, 32);
assert_eq!(EnzymeType::all().len(), 16);

// 2. In silico digestion → tags (multi-contig aware)
// BcgI: 32-bp window, anchors at offset 10 (CGA) and 19 (TGC)
let seq = b"AAAAAAAAAACGAAAAAAATGCAAAAAAAAAA";
let tags = digest_genome_contig(seq, EnzymeType::BcgI, 1, 0);

// 3. Assemble a TGT record (gaps auto-computed from positions)
let mut rec = TgtRecord::new("genome_a", 4_641_652);
for t in tags { rec.add_tag(t); }
rec.contig_names = vec!["NC_000913".to_string()];
rec.contig_offsets = vec![0];
println!("tags={}, mean_gap={:.1}", rec.tag_count(), rec.mean_gap());

// 4. Persist / load (text or binary)
{
    let mut writer = TgtWriter::new(Path::new("genome_a.tgt"))?;
    writer.write_record(&rec)?;
}
let mut reader = TgtReader::new(Path::new("genome_a.tgt"))?;
while let Some(r) = reader.read_record()? { /* ... */ }

// 5. Build the adjacency graph across genomes and extract synteny
let mut g = TagAdjacencyGraph::new();
g.add_genome("genome_a", &rec);
// g.add_genome("genome_b", &rec_b);  // load another record
g.build_edges();
g.simplify(2);                       // keep adjacencies supported by ≥2 genomes
for path in g.linear_paths() {
    println!("backbone score = {:.3}", synteny_score(&path, &g));
}
let blocks = extract_synteny_blocks(&g);
```

`TgtRecord` also offers `median_gap()`, `max_gap()`, and `coverage_fraction()`
(estimated fraction of the genome covered by tag bases). `Tag` provides
`sequence_str()` and `hamming_distance()`. Utilities in `bsyn::utils` include
`reverse_complement`, `is_valid_dna`, and `gc_content`.

The crate exposes its functionality as the `bsyn` library. Representative API:

```rust
use bsyn::enzyme::EnzymeType;
use bsyn::enzyme::enzyme::Enzyme;
use bsyn::enzyme::digest::digest_genome_contig;
use bsyn::tgt::{Tag, TgtRecord, TgtReader, TgtWriter, Strand};
use bsyn::synteny::{TagAdjacencyGraph, extract_synteny_blocks, synteny_score};
use std::path::Path;

// 1. Enzyme properties
let bcgi = Enzyme::properties(EnzymeType::BcgI);
assert_eq!(bcgi.tag_length, 32);
assert_eq!(EnzymeType::all().len(), 16);

// 2. In silico digestion → tags (multi-contig aware)
let tags = digest_genome_contig(b"CGA......TGC...", EnzymeType::BcgI, 1, 0);

// 3. Assemble a TGT record (gaps auto-computed from positions)
let mut rec = TgtRecord::new("genome_a", 4_641_652);
for t in tags { rec.add_tag(t); }
rec.contig_names = vec!["NC_000913".to_string()];
rec.contig_offsets = vec![0];
println!("tags={}, mean_gap={:.1}", rec.tag_count(), rec.mean_gap());

// 4. Persist / load (text or binary)
TgtWriter::new(Path::new("genome_a.tgt"))?.write_record(&rec)?;
let mut reader = TgtReader::new(Path::new("genome_a.tgt"))?;
while let Some(r) = reader.read_record()? { /* ... */ }

// 5. Build the adjacency graph across genomes and extract synteny
let mut g = TagAdjacencyGraph::new();
g.add_genome("genome_a", &rec);
g.add_genome("genome_b", &rec_b);
g.build_edges();
g.simplify(2);                       // keep adjacencies supported by ≥2 genomes
for path in g.linear_paths() {
    println!("backbone score = {:.3}", synteny_score(&path, &g));
}
let blocks = extract_synteny_blocks(&g);
```

`TgtRecord` also offers `median_gap()`, `max_gap()`, and `coverage_fraction()`
(estimated fraction of the genome covered by tag bases). `Tag` provides
`sequence_str()` and `hamming_distance()`. Utilities in `bsyn::utils` include
`reverse_complement`, `is_valid_dna`, and `gc_content`.

---

## Synteny algorithm

Implemented in `src/synteny/graph.rs` and `src/synteny/blocks.rs`:

1. **Graph construction** — each genome's TGT record is ingested with
   `add_genome`. Identical tag sequences are de-duplicated into a single
   **`TagNode`**, which records `(position, strand, contig_id)` per genome.
   `build_edges` then walks each genome's tag order and creates a directed
   **`AdjacencyEdge`** for every consecutive pair; the edge **weight** is the
   number of distinct genomes exhibiting that adjacency, and it stores the
   inter-tag distance per genome.

2. **Simplification** — `simplify(min_weight)` drops edges supported by fewer
   than `min_weight` genomes (noise / unanchored adjacencies), then removes nodes
   left with degree 0. This is the ntSynt-inspired anchoring step.

3. **Backbone extraction** — `linear_paths()` returns maximal chains of degree-2
   nodes. Each chain is a **synteny backbone**: a stretch of tags whose adjacency
   is consistently conserved. Cyclic components are handled separately.

4. **Block extraction** — `extract_synteny_blocks` turns each backbone into a
   **`SyntenyBlock`** carrying, for every genome that contains the whole chain,
   the `(start, end, strand)` coordinates. Blocks require ≥2 genomes.
   `filter_blocks_by_size` keeps blocks above a minimum genomic span.

5. **Indel detection** — `detect_indels` compares the inter-tag distances of a
   block across genome pairs; differences beyond a threshold (≥10% of the mean
   distance and ≥100 bp) are reported as insertions/deletions.

Common-tag queries are also available: `find_common_tags()` (exact sequence
match across all genomes) and `find_common_tags_tolerance(hamming)`
(union-find clustering within a Hamming radius).

---

## Scoring metrics

From `src/synteny/scoring.rs`:

| Function | Input | Range | Meaning |
|---|---|:---:|---|
| `synteny_score(path, graph)` | a backbone path | 0–1 | mean edge weight ÷ #genomes, times a √-length bonus that saturates at 10 tags |
| `pairwise_synteny_matrix(graph)` | whole graph | 0–1 | all-vs-all Jaccard of per-genome adjacency sets |
| `adjacency_jaccard(rec_a, rec_b)` | two records | 0–1 | Jaccard of adjacent tag-sequence pairs (works without a graph) |
| `kendall_tag_order(rec_a, rec_b)` | two records | −1–1 | Kendall's τ rank correlation on the order of shared tags |
| `breakpoint_count(rec_a, rec_b)` | two records | ≥0 | number of adjacencies present in one genome but not the other (symmetric difference) |
| `windowed_synteny_score(...)` | two records | 0–1 | sliding-window Jaccard + position correlation for local synteny |

`adjacency_jaccard`, `kendall_tag_order`, and `breakpoint_count` operate directly
on tag *sequences*, so they can compare two `TgtRecord`s without first building a
graph.

---

## Project status & roadmap

**What is implemented and unit-tested**

- TGT data model v2: `Tag` (with `contig_id`), `Gap`, `TgtRecord` (with
  `contig_names`, `contig_offsets`, gap statistics, coverage fraction).
- TGT text I/O (`TgtReader::read_record` / `TgtWriter::write_record`) and
  binary I/O (`read_binary` / `write_binary`) with round-trip parsing and gap
  validation.
- Streaming FASTA reader (`FastaReader`).
- Enzyme catalog: 16 `EnzymeType` variants with corrected tag lengths
  (27–33 bp), index round-trips, and IUPAC degenerate base support.
- In silico digestion: `digest_genome` (single contig) and
  `digest_genome_contig` (multi-contig with cumulative offset).
- Tag adjacency graph: ingestion, edge building, simplification, linear-path and
  common-tag extraction (including inversion-break behavior).
- Synteny blocks, indel detection, size filtering.
- All scoring metrics (Jaccard, Kendall τ, breakpoint count, windowed scores).
- **Scaffold subcommand**: draft-to-reference contig mapping with orientation
  detection and AGP v2.1 output.
- DNA utilities (`reverse_complement`, `is_valid_dna`, `gc_content`).
- **Pluggable landmark sources**: 2bRAD enzyme digestion or FracMinHash
  sketches (`--mode fracminhash`), producing the same TGT downstream.
- **Structural-synteny metrics**: `breakpoint_count`, `inverted_fraction`
  (majority-frame and fixed-reference), `observable_fraction`, and junction
  coordinates.
- Integration tests: 18 tests covering enzyme catalog, TGT round-trip, binary
  I/O, digestion, graph creation, FASTA parsing, CLI help, and FracMinHash
  round-trip.

**Completed milestones**

1. **Performance optimization.** Digestion now uses `memchr` skip-based search
   and per-enzyme `rayon` parallelism; the four-enzyme panel is substantially
   faster than the byte-by-byte implementation.
2. **Detection-power characterization.** Measured L50 resolution on *E. coli*
   K-12: ~2.6 kb for inversions (BcgI) and ~1.2 kb for the four-enzyme panel;
   translocations require roughly half the length. See
   `docs/PHASE2_DETECTION_POWER.md`.
3. **Mathematical analysis of fragmentation.** `breakpoint_count` now requires
   positive contradiction, making it invariant to contig count; the residual
   information loss is captured by `observable_fraction`. See
   `docs/MATH_REVIEW.md`.
4. **GTDB-scale validation.** Inverted fraction agrees with dnadiff at r = 0.936
   across 43,334 held-out pairs (r = 0.996 for ≥97% ANIm), with a closed-form
   standard-error model. See `docs/MATH_REVIEW.md` §7.
5. **SynTracker cohort replication.** Four published isolate collections (E. coli
   hypermutator, H. pylori, N. gonorrhoeae, S. rimosus) reproduced the expected
   SNP-vs-SV evolutionary signatures.
6. **Syn2bANI application benchmarks.** The same synteny engine underpins
   Syn2bANI, which adds an MLE ANI estimator and a pairwise `struct` wrapper.
   Application case studies in the Syn2bANI-paper repository demonstrate
   phenotype-relevant architecture differences that ANI alone misses:
   *H. pylori* cagPAI rearrangement states across Correa's cascade,
   *B. longum* abfA gene-cluster deletion in a constipation-treatment cohort,
   *E. coli* O157:H7 host-source differentiation, and FDA-ARGOS *S. aureus*
   lineage structure. See https://github.com/HuangShiLab/Syn2bANI-paper.

**Current priorities**

1. **Closed-genome phenotype cohorts.** Select near-closed genomes from GTDB for
   species with known large-inversion biology (*S. pneumoniae*, *S. enterica*,
   *B. pertussis*, *P. aeruginosa*) and validate junction coordinates against
   dnadiff.
2. **FracMinHash design-rule experiments.** Test whether the error-model
   constants transfer across landmark sources, or are panel-specific.
3. **Syn2bANI integration.** Keep the Syn2b synteny engine aligned with the
   Syn2bANI pairwise `struct` wrapper and shared TGT format; maintain
   cross-references between `docs/MATH_REVIEW.md` / `docs/PHASE2_DETECTION_POWER.md`
   and the Syn2bANI manuscript.

---

## Validation & limitations

A benchmark against the **real SynTracker** (Enav, Paz & Ley, *Nat Biotechnol*
2024) on **10 complete, single-contig *C. acnes* genomes** (45 pairs) produced a
sobering, important result:

- The Syn2b metric that *appears* to validate the approach — **tag-Jaccard**
  (presence/absence), raw *r* = 0.98 vs SynTracker — is **sequence identity in
  disguise**: it is ~indistinguishable from Mash distance (*r* = 0.997), and its
  **partial correlation with SynTracker controlling for Mash collapses to
  ≈ +0.06**. It adds essentially nothing beyond sequence divergence. This is what
  produced earlier, over-optimistic "r ≈ 0.85" claims.
- The **genuine** synteny metric — **tag adjacency** — carries a *real but weak*
  structural signal: raw *r* = 0.67 (Spearman 0.79), **partial *r* ≈ 0.37**
  controlling for Mash (*p* ≈ 0.015). But its dynamic range on co-linear genomes
  is tiny (0.000–0.0025), so discrimination power is low.
- On these genomes SynTracker's own APSS is ~98% explained by sequence
  divergence, i.e. synteny was not an independent strain-discriminating axis for
  this species.

A subsequent four-species replication on published isolate cohorts (Enav et al.
2024) produced the expected evolutionary signatures: *E. coli* hypermutator
(SNP-driven, anchor adjacency flat despite wide ANI range), *S. rimosus*
(clonal ANI but highly variable anchor adjacency, SV-driven), *N. gonorrhoeae*
and *H. pylori* (mixed SNP + SV modes). These results are held in the
Syn2bANI-paper repository under `results/syntracker_validation/`.

**Practical limitations** of the tag-order approach:

- **Resolution ≈ 1–2 kb.** With ~1,660 bp average tag spacing, rearrangements
  smaller than the gap are invisible.
- **32 bp tags** carry limited sequence information and may map ambiguously to
  paralogs/repeats.
- **Reference-dependent ordering** — tag coordinates require a reference; de novo
  comparison is not supported.
- **Sparse representation** — only ~1.5% of the genome is sampled.

**Takeaway:** "correlates with SynTracker" is *not* evidence of capturing
synteny, because sequence divergence itself correlates with APSS across a
phylogeny. Any synteny claim must use **partial correlation** (control for
ANI/Mash) or otherwise **decouple structure from SNP load** — e.g. simulate known
inversions/translocations on a single complete genome at fixed divergence and
show that tag-adjacency tracks the structural change while Mash does not. That is
the recommended next experiment.

See the top-level repository for the full analyses:
`SYNTRACKER_vs_Syn2b_COMPLETE_REPORT.md`, `executive_summary.md`, and
`syntracker_vs_syn2b_COMPLETE.png`.

---

## Testing

```bash
cd syn2b
cargo test            # unit + integration tests (22 passed, 0 failed)
cargo test --release  # optimized tests
cargo build --release # release build (0 errors, 0 warnings)
```

Every library module carries unit tests (enzyme catalog, digestion, TGT text/binary
I/O and gap validation, FASTA parsing, graph construction/simplification/
paths, block extraction/indels, and all scoring metrics).
`tests/integration_tests.rs` covers end-to-end scenarios including binary
round-trip and CLI validation.

---

## References

- **SPEC:** [`../SPEC.md`](../SPEC.md) — full design specification for Syn2b.
- **ntSynt** — minimizer-graph synteny detection (inspiration for the tag
  adjacency graph).
- **KmerAperture** (2024) — ordered k-mer series for structural variation.
- **SynTracker** — Enav, Paz & Ley, *Nature Biotechnology*, 2024 (validation
  reference / gold standard).
- **2bRAD-M** — reduced-representation microbial profiling (source of the tags).
- **Fast2bRAD-M** — high-performance 2bRAD tag extraction (reference for tag
  lengths and recognition-site definitions).

---

## License

MIT (see `Cargo.toml`). © Syn2b Contributors.
