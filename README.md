# Syn2b — 2bRAD-based Synteny Detection Engine

**2bSyn** is a Rust, alignment-free engine for detecting genome synteny and
structural variation from **2bRAD tags** — the short, fixed-position sequences
produced by *Type IIB* restriction enzymes. Instead of aligning whole genomes,
2bSyn represents each genome as an ordered series of sparse anchor tags (the
**Tag–Gap–Tag / TGT** model) and infers synteny from how the *order and
adjacency* of those tags are conserved across genomes.

> **Status: research prototype.** The crate builds (`cargo build --release`) and
> its core data structures, TGT I/O, adjacency graph, and scoring metrics are
> implemented and unit-tested. The command-line pipeline, however, is still a
> scaffold and the in-silico digestion is not yet complete. See
> [Project status & roadmap](#project-status--roadmap) before relying on it.

---

## Table of contents

- [Motivation](#motivation)
- [How it works](#how-it-works)
- [Features](#features)
- [Repository layout](#repository-layout)
- [The TGT format](#the-tgt-format)
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
with roughly **1,000–3,000 species-specific 32 bp tags**, each anchored at a
fixed, restriction-map-determined coordinate. These tags are essentially a very
sparse, deterministic *minimizer set*: their **relative order is fixed by the
genome**, so a rearrangement (inversion, translocation, large indel) shows up as
a change in which tags are adjacent.

Existing synteny tools cannot consume this data directly. SynTracker, for
example, needs ~5 kb homologous regions with flanking context for BLAST/DECIPHER
alignment — a 32 bp isolated tag carries neither the length nor the context it
requires. 2bSyn is a **dedicated method** that works natively on sparse tag
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
   2bRAD tags  ──────────────►  TGT record  (one per genome)
        │                        Tag–Gap–Tag: ordered tags + inter-tag gaps
        │  text or binary I/O   (tgt/reader.rs, tgt/writer.rs)
        ▼
   TagAdjacencyGraph            (synteny/graph.rs)
     • nodes  = unique tag sequences (shared across genomes)
     • edges  = adjacencies, weighted by #supporting genomes
        │
        │  simplify()  → drop low-support edges & isolated nodes
        │  linear_paths()  → maximal degree-2 chains = synteny backbones
        ▼
   Synteny blocks + scores      (synteny/blocks.rs, synteny/scoring.rs)
     • genomic coordinates per genome
     • indel detection from gap-size differences
     • Jaccard / Kendall τ / breakpoint / path scores
```

---

## Features

- **In silico digestion** with 16 Type IIB restriction enzymes (IUPAC-degenerate
  recognition sites supported).
- **TGT (Tag–Gap–Tag) representation** — a compact, structured genome model of
  ordered tags plus the base-pair gaps between them.
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

---

## Repository layout

```
Syn2b/
├── Cargo.toml                  # package / lib / bin all named `Syn2b`
├── src/
│   ├── main.rs                # CLI entry point (clap): digest / synteny / coverage / convert
│   ├── lib.rs                 # public API re-exports (crate name: `Syn2b`)
│   ├── tgt/                   # TGT core data structures
│   │   ├── tag.rs             # Tag  (32 bp sequence + position + enzyme + strand)
│   │   ├── gap.rs             # Gap  (inter-tag distance in bp)
│   │   ├── record.rs          # TgtRecord (one genome: ordered tags + gaps + stats)
│   │   ├── writer.rs          # TgtWriter (text + binary output)
│   │   └── reader.rs          # TgtReader (text + binary input)
│   ├── enzyme/                # Restriction-enzyme definitions
│   │   ├── enzyme.rs          # EnzymeType (16 variants) + Enzyme properties
│   │   └── digest.rs          # in silico digestion (recognition-site search + tag extraction)
│   ├── synteny/               # Synteny-detection engine
│   │   ├── graph.rs           # TagAdjacencyGraph, TagNode, AdjacencyEdge
│   │   ├── scoring.rs         # synteny_score, pairwise matrix, Jaccard, Kendall τ, breakpoints
│   │   └── blocks.rs          # SyntenyBlock extraction, indel detection, size filtering
│   ├── io/
│   │   └── fasta.rs           # streaming FASTA reader (FastaReader / FastaRecord)
│   └── utils/
│       └── mod.rs             # reverse_complement, is_valid_dna, gc_content
└── tests/
    └── integration_tests.rs   # end-to-end tests (see status note below)
```

The Cargo package, the library crate, and the binary are all named **`Syn2b`**
(a Cargo *package* name may not start with a digit, which is why the project is
`Syn2b` and not `2bsyn`; the library and binary targets share the same name).

---

## The TGT format

A **TGT record** describes one genome as an ordered list of tags with the gaps
between them. `gaps[i]` is the distance (bp) between `tags[i]` and `tags[i+1]`,
so a record with *N* tags has *N − 1* gaps.

### Text format

```
>NC_000913|length=4641652
BcgI:ATCG… -1313- GCTA… -1298- TTAA…
AlfI:CGAT… -892-  AATT… -1567- GCAT…
```

- A header line `>genome_id|length=<bp>` opens each record.
- Tag lines are grouped by enzyme; each tag is written as `Enzyme:SEQUENCE`.
- Gaps appear between consecutive tags as `-<size>-`.

The reader reconstructs tag positions from the cumulative gaps and validates that
each parsed gap matches the value recomputed from tag positions.

### Binary format

A fixed-layout, little-endian binary encoding for efficient storage:

```
Header (32 bytes)
  [0..4]    Magic          b"TGT\x01"
  [4..8]    Version        u32 (= 1)
  [8..16]   Genome length  u64
  [16..20]  Tag count      u32
  [20..22]  Enzyme count   u16
  [22..32]  Reserved

Tag table (tag_count × 48 bytes)
  [0..32]   Tag sequence   raw bytes (zero-padded)
  [32..40]  Position       u64
  [40..41]  Enzyme index   u8 (0–15)
  [41..42]  Strand         u8 (0 = forward, 1 = reverse)
  [42..48]  Reserved

Gap table ((tag_count − 1) × 4 bytes)
  [0..4]    Gap size       u32
```

On read, the magic bytes and version are verified and the stored gap table is
cross-checked against gaps recomputed from tag positions.

---

## Type IIB restriction enzymes

Type IIB enzymes cut on **both** sides of their recognition site, excising a
short, defined fragment — ideal for producing fixed-length 2bRAD tags. 2bSyn
ships definitions for 16 of them (`src/enzyme/enzyme.rs`). Recognition sites use
IUPAC codes (`N` = any; `R` = A/G; `Y` = C/T; etc.), which the digester matches
via `matches_site`.

| Enzyme  | Recognition site   | Tag length | 5′/3′ cut offset |
|---------|--------------------|:----------:|:----------------:|
| BcgI    | `CGANNNNNNTGC`     | 32         | −10 / +10        |
| AlfI    | `GCANNNNNNTGC`     | 32         | −9 / +9          |
| AloI    | `GAACNNNNNNTCC`    | 32         | −11 / +11        |
| BaeI    | `ACNNNNGTAYC`      | 32         | −12 / +12        |
| BplI    | `GAGNNNNNCTC`      | 32         | −8 / +8          |
| BsaXI   | `ACNNNNNCTCC`      | 32         | −10 / +10        |
| BslFI   | `GGGAC`            | 28         | −7 / +7          |
| Bsp24I  | `GACNNNNNNTGG`     | 32         | −11 / +11        |
| CjeI    | `RYANNNNNNCTC`     | 32         | −10 / +10        |
| CjePI   | `GCANNNNNNGTG`     | 32         | −10 / +10        |
| CspCI   | `CAANNNNNGTGG`     | 32         | −11 / +11        |
| FalI    | `AAGNNNNNCTT`      | 32         | −9 / +9          |
| HaeIV   | `GAYNNNNNRTC`      | 32         | −9 / +9          |
| Hin4I   | `GAYNNNNNVTC`      | 32         | −10 / +10        |
| PpiI    | `GAACNNNNNCTC`     | 32         | −10 / +10        |
| PsrI    | `GAACNNNNNNTAC`    | 32         | −11 / +11        |

The recognition sites and cut offsets above are the values encoded in the
source; they are a simplified in-silico model rather than a full biochemical
specification (see [Project status](#project-status--roadmap)).

---

## Building

Requires a stable Rust toolchain (Rust 2021 edition; install via
[rustup](https://rustup.rs)).

```bash
cd Syn2b
cargo build --release
```

This compiles cleanly (only a couple of dead-code/unused-import warnings) and
produces the `Syn2b` binary at `target/release/Syn2b`. The Cargo package, the
library crate, and the binary target are all named `Syn2b` — a package name may
not start with a digit, which is why the project is not called `2bsyn`.

Dependencies (from `Cargo.toml`): `clap` (CLI), `anyhow` (errors), `serde` +
`bincode` (serialization), `parking_lot`, and `rayon` (parallelism);
`tempfile` for tests.

---

## Command-line usage

```
Syn2b <COMMAND> [OPTIONS]

Commands:
  digest    In silico digest genomes with 2bRAD enzymes
  synteny   Compute synteny between genomes using TGT
  coverage  Analyze multi-enzyme coverage statistics
  convert   Convert between TGT text and binary formats
```

Examples:

```bash
# Digest a genome with the default enzyme (BcgI) into a text TGT file
Syn2b digest -i genome.fasta -o genome.tgt

# Digest with all 16 enzymes, writing the compact binary format
Syn2b digest -i genomes/ -o out/ --enzymes all --format binary

# Multi-enzyme coverage statistics
Syn2b coverage -i genome.fasta --enzymes all

# Compute synteny across a set of TGT records
Syn2b synteny -i tgts/ -o synteny_report.txt

# Convert a text TGT to binary (or vice versa)
Syn2b convert -i genome.tgt -o genome.btgt --format binary
```

Common options: `-i/--input`, `-o/--output`, `-e/--enzymes`
(comma-separated names or `all`), `-f/--format` (`text` | `binary`).

> The subcommands currently parse and validate their arguments (including the
> enzyme list) and print what they *would* do — the end-to-end
> digest→TGT→graph→report orchestration is not yet wired up. The building blocks
> exist as library functions; see below.

---

## Library usage

The crate exposes its functionality as the `Syn2b` library. Representative API:

```rust
use Syn2b::enzyme::EnzymeType;
use Syn2b::enzyme::enzyme::Enzyme;
use Syn2b::enzyme::digest::digest_genome;
use Syn2b::tgt::{Tag, Gap, TgtRecord, TgtReader, TgtWriter, Strand};
use Syn2b::synteny::{TagAdjacencyGraph, extract_synteny_blocks, synteny_score};
use std::path::Path;

// 1. Enzyme properties
let bcgi = Enzyme::properties(EnzymeType::BcgI);   // recognition_site, tag_length, cut offsets
assert_eq!(EnzymeType::all().len(), 16);

// 2. In silico digestion → tags
let tags = digest_genome(b"CGA......TGC...", EnzymeType::BcgI);

// 3. Assemble a TGT record (gaps are auto-computed from positions)
let mut rec = TgtRecord::new("genome_a", 4_641_652);
for t in tags { rec.add_tag(t); }
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

`TgtRecord` also offers `median_gap()`, `max_gap()`, `coverage_fraction()`
(estimated fraction of the genome covered by tag bases), and `enzyme_count()`.
`Tag` provides `sequence_str()` and `hamming_distance()`; `Strand` provides
`to_u8()` / `from_u8()` for the binary format. Utilities in `Syn2b::utils` include
`reverse_complement`, `is_valid_dna`, and `gc_content`.

> Note: the `Enzyme::properties(...)` / `TagAdjacencyGraph::new()` signatures
> above reflect the **current source**. The [SPEC](../SPEC.md) and
> `tests/integration_tests.rs` describe a slightly different target API
> (`EnzymeType::properties(&self)`, `TagAdjacencyGraph::new(num_genomes)`,
> `digest_multi_enzyme`, …) that the library has not fully converged on yet.

---

## Synteny algorithm

Implemented in `src/synteny/graph.rs` and `src/synteny/blocks.rs`:

1. **Graph construction** — each genome's TGT record is ingested with
   `add_genome`. Identical 32 bp tag sequences are de-duplicated into a single
   **`TagNode`**, which records `(position, strand)` per genome. `build_edges`
   then walks each genome's tag order and creates a directed **`AdjacencyEdge`**
   for every consecutive pair; the edge **weight** is the number of distinct
   genomes exhibiting that adjacency, and it stores the inter-tag distance per
   genome.

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

`adjacency_jaccard`, `kendall_tag_order`, and `breakpoint_count` operate directly
on tag *sequences*, so they can compare two `TgtRecord`s without first building a
graph.

---

## Project status & roadmap

**What is implemented and unit-tested**

- TGT data model: `Tag`, `Gap`, `TgtRecord` (gap statistics, coverage fraction,
  enzyme count, text `Display`).
- TGT text I/O (`TgtReader::read_record` / `TgtWriter::write_record`) and the
  fixed-layout binary I/O (`read_binary` / `write_binary`, incl. `Strand`
  encode/decode), with gap validation.
- Streaming FASTA reader (`FastaReader`).
- Enzyme catalog: 16 `EnzymeType` variants, index round-trips, and per-enzyme
  `Enzyme::properties`.
- Tag adjacency graph: ingestion, edge building, simplification, linear-path and
  common-tag extraction (including inversion-break behavior).
- Synteny blocks, indel detection, size filtering.
- All scoring metrics.
- DNA utilities (`reverse_complement`, `is_valid_dna`, `gc_content`).

**Known gaps**

- **In-silico digestion is incomplete.** The tag window derived from
  `cut_offset_5`/`cut_offset_3` does not currently match each enzyme's declared
  `tag_length` (e.g. BcgI's ±10 offsets yield a 20 bp window vs. `tag_length =
  32`), so `digest_genome` does not yet emit correct-length tags on real
  sequence. The recognition-site search and IUPAC matching themselves work.
- **CLI subcommands are scaffolds** — they validate arguments and print intended
  actions but do not yet run the digest→TGT→graph→report pipeline.
- **Multi-enzyme digestion** (`digest_multi_enzyme`, `merge_multi_enzyme_tags`)
  and some type signatures described in [`SPEC.md`](../SPEC.md) / exercised by
  `tests/integration_tests.rs` are not yet implemented, so the integration-test
  target does not compile against the current library.
- **Seven library unit tests fail** for pre-existing reasons unrelated to the
  build: five test *helpers* overflow `u8` (`i * 100` with `i: u8`) in
  `scoring.rs` / `graph.rs`, and two text-format assertions compare against
  `sequence_str()`, which returns the full zero-padded 32-byte array. **58 of 65
  pass.**

**Roadmap**

1. Fix the digestion window so tag length matches each enzyme definition; validate
   against the SPEC targets (~1,169 BcgI tags on *E. coli* K-12; ~8.4% combined
   coverage for 16 enzymes).
2. Add multi-enzyme digestion and wire the CLI subcommands to the library.
3. Reconcile the library API with `SPEC.md` / the integration tests, and fix the
   pre-existing unit-test failures.
4. Add the structural-decoupling benchmark described below.

---

## Validation & limitations

A benchmark against the **real SynTracker** (Enav, Paz & Ley, *Nat Biotechnol*
2024) on **10 complete, single-contig *C. acnes* genomes** (45 pairs) produced a
sobering, important result:

- The 2bSyn metric that *appears* to validate the approach — **tag-Jaccard**
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
`SYNTRACKER_vs_2bSyn_COMPLETE_REPORT.md`, `executive_summary.md`, and
`syntracker_vs_2bsyn_COMPLETE.png`.

---

## Testing

```bash
cd Syn2b
cargo test --lib      # library unit tests
```

Every library module carries unit tests (enzyme catalog, digestion, TGT text and
binary I/O with gap validation, FASTA parsing, graph construction / simplification
/ paths, block extraction / indels, and the scoring metrics). At present **58 of
65 unit tests pass**; the 7 failures are pre-existing issues unrelated to the
build fix — mostly `u8` overflow inside test helpers (see
[Known gaps](#project-status--roadmap)). The higher-level
`tests/integration_tests.rs` targets the API from `SPEC.md` and does not yet
compile against the current library, so use `--lib` (or `cargo build`) for now
rather than a bare `cargo test`.

---

## References

- **SPEC:** [`../SPEC.md`](../SPEC.md) — full design specification for 2bSyn.
- **ntSynt** — minimizer-graph synteny detection (inspiration for the tag
  adjacency graph).
- **KmerAperture** (2024) — ordered k-mer series for structural variation.
- **SynTracker** — Enav, Paz & Ley, *Nature Biotechnology*, 2024 (validation
  reference / gold standard).
- **2bRAD-M** — reduced-representation microbial profiling (source of the tags).

---

## License

MIT (see `Cargo.toml`). © Syn2b Contributors.
