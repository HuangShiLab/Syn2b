# Syn2b Project Context — 2026-01-07

> **For**: Agent Swarm continuation (Syn2b development + paper writing)
> **Status**: Core pipeline complete, performance optimization & paper writing next
> **Repo**: https://github.com/HuangShiLab/Syn2b

---

## 1. Project Overview

**Syn2b** is a Rust-based, alignment-free synteny detection engine that leverages 2bRAD tags (Type IIB restriction enzyme fragments) to detect genomic structural variations between microbial genomes. Instead of aligning whole genomes, Syn2b represents each genome as an ordered series of sparse anchor tags (the **Tag-Gap-Tag / TGT** model) and infers synteny from tag adjacency conservation.

- **Crate name**: `bsyn` (library), `syn2b` (binary)
- **Language**: Rust 2021 edition
- **Repository**: `https://github.com/HuangShiLab/Syn2b`
- **Local path**: `/Users/shihuang/Downloads/Syn2b` (also symlinked as `syn2b`)

---

## 2. Architecture & Repository Layout

```
syn2b/
├── Cargo.toml
├── README.md
├── SPEC.md
├── src/
│   ├── main.rs          # CLI entry (clap): digest / synteny / scaffold / coverage / convert
│   ├── lib.rs           # Public API re-exports (crate = bsyn)
│   ├── tgt/             # TGT data structures
│   │   ├── tag.rs       # Tag (32bp seq + position + enzyme + strand + contig_id)
│   │   ├── gap.rs       # Gap (inter-tag distance in bp)
│   │   ├── record.rs    # TgtRecord (ordered tags + gaps + contig metadata)
│   │   ├── writer.rs    # TgtWriter (text + binary output, v2 format)
│   │   └── reader.rs    # TgtReader (text + binary input, v2 format)
│   ├── enzyme/          # Restriction enzyme definitions
│   │   ├── enzyme.rs    # 16 EnzymeType variants + Enzyme properties + IUPAC matching
│   │   └── digest.rs    # in silico digestion (pattern matching + tag extraction)
│   ├── synteny/         # Synteny detection engine
│   │   ├── graph.rs     # TagAdjacencyGraph, TagNode, AdjacencyEdge
│   │   ├── scoring.rs   # Jaccard, Kendall τ, breakpoints, windowed scores
│   │   └── blocks.rs    # SyntenyBlock extraction, indel detection, size filtering
│   ├── io/
│   │   └── fasta.rs     # Streaming FASTA reader (FastaReader / FastaRecord)
│   └── utils/
│       └── mod.rs       # reverse_complement, is_valid_dna, gc_content
└── tests/
    └── integration_tests.rs  # 17 end-to-end tests
```

---

## 3. Implemented Features

### 3.1 TGT Data Model (v2)
- `Tag` struct: 32-byte sequence, position (u64), enzyme index, strand (FWD/REV), contig_id (u16)
- `Gap` struct: inter-tag distance (u32)
- `TgtRecord`: ordered tags + auto-computed gaps + contig metadata (`contig_names`, `contig_offsets`)
- Text format: `>genome_id|length=NNN` header, `#contigs=name:length;...` comment, `Enzyme:SEQUENCE` tags with `-<gap>-` separators
- Binary format: fixed-layout little-endian, magic `TGT\x01`, version 1, gap validation on read

### 3.2 Enzyme Catalog (16 Type IIB Enzymes)
| Enzyme | Tag Len | IUPAC | Notes |
|--------|---------|-------|-------|
| BcgI | 32 | — | Two anchors: CGA@10, TGC@19 |
| AlfI | 32 | — | Palindrome: GCA/TGC |
| AloI | 27 | — | |
| BaeI | 28 | Y (fwd), R (rev) | `[CT]`@19 fwd, `[AG]`@8 rev |
| BplI | 27 | — | Palindrome: GAG/CTC |
| BsaXI | 27 | — | |
| BslFI | 25 | — | |
| Bsp24I | 27 | — | |
| CjeI | 28 | — | |
| CjePI | 27 | — | |
| CspCI | 33 | — | |
| FalI | 27 | — | Palindrome: AAG/CTT |
| HaeIV | 27 | Y + R | Y@9+R@15 fwd; Y@11+R@17 rev |
| Hin4I | 27 | Y + [GAC] | Y@10+[GAC]@16 fwd; [CTG]@10+R@16 rev |
| PpiI | 27 | — | |
| PsrI | 27 | — | |

- IUPAC bitmask: bit0=A, bit1=T, bit2=C, bit3=G
- Pattern matching: `Anchor` (exact motif at offset) + `IupacConstraint` (bitmask at offset)
- All tag lengths **cross-validated against Fast2bRAD-M** reference implementation

### 3.3 In Silico Digestion
- `digest_genome(sequence, enzyme)` → single contig
- `digest_genome_contig(sequence, enzyme, contig_id, offset)` → multi-contig aware with cumulative offset
- Algorithm: O(N × L × P) byte-by-byte sliding window (one window per enzyme pattern per position)
- Filters: anchor match + IUPAC constraint + `is_pure_atcg` (no N/degenerate bases in tag)
- Tags sorted by position and deduplicated

### 3.4 Synteny Detection Engine
- **TagAdjacencyGraph**:
  - `add_genome()` → deduplicate tags into `TagNode` (records per-genome position/strand/contig_id)
  - `build_edges()` → create `AdjacencyEdge` for consecutive tag pairs per genome; weight = #genomes sharing adjacency
  - `simplify(min_weight)` → drop low-support edges, remove isolated nodes (ntSynt-inspired)
  - `linear_paths()` → extract maximal degree-2 chains (synteny backbones)
- **SyntenyBlock extraction**: `extract_synteny_blocks()` → blocks with per-genome (start, end, strand) coordinates
- **Indel detection**: `detect_indels()` compares inter-tag distances across genome pairs; threshold: ≥10% of mean AND ≥100 bp
- **Size filtering**: `filter_blocks_by_size()`
- **Common tag queries**: `find_common_tags()` (exact), `find_common_tags_tolerance(hamming)` (union-find clustering)

### 3.5 Scoring Metrics (src/synteny/scoring.rs)
| Function | Input | Range | Meaning |
|---|---|:---:|---|
| `synteny_score(path, graph)` | backbone | 0–1 | mean edge weight ÷ #genomes × √-length bonus (saturates at 10 tags) |
| `pairwise_synteny_matrix(graph)` | whole graph | 0–1 | all-vs-all Jaccard of adjacency sets |
| `adjacency_jaccard(rec_a, rec_b)` | two records | 0–1 | Jaccard of adjacent tag-sequence pairs |
| `kendall_tag_order(rec_a, rec_b)` | two records | −1–1 | Kendall's τ on shared tag order |
| `breakpoint_count(rec_a, rec_b)` | two records | ≥0 | symmetric adjacency difference |
| `windowed_synteny_score(...)` | two records | 0–1 | sliding-window Jaccard + position correlation |

### 3.6 Scaffold Subcommand (main.rs)
1. Load reference and draft TGT records
2. For each draft contig, evaluate FWD and REV orientation by matching tags against reference
3. Count-ratio heuristic (>2× difference) picks dominant orientation
4. Sort contigs by median reference position
5. Output AGP v2.1 with real contig lengths and estimated gap sizes

- **Validated**: E. coli K-12 self-scaffold (4 reversed contigs correctly identified)
- **ABHQ draft**: 135 contigs → 45 anchored at `min_tags=3`

### 3.7 CLI Commands
```bash
syn2b digest -i genome.fasta -o genome.tgt
syn2b digest -i genome.fasta -o genome.btgt --enzymes all --format binary
syn2b synteny -i tgts/ -o synteny_matrix.csv
syn2b scaffold -r reference.tgt -d draft.tgt -o scaffolds.agp --min-tags 3
syn2b coverage -i genome.fasta --enzymes all
syn2b convert -i genome.tgt -o genome.btgt --format binary
```

---

## 4. Current Test Status

```bash
cargo test            # 99 tests passed, 0 failed
cargo build --release # 0 errors, 0 warnings
```

**Integration tests** (tests/integration_tests.rs):
- Enzyme catalog: 16 variants, correct tag lengths, index round-trip, all() count
- TGT round-trip: text write → read → verify tags/gaps
- Binary I/O: write_binary → read_binary → gap validation, magic/version verification
- Digestion: empty sequence, no sites, real BcgI/CjeI/CjePI sites, N-rejection, degenerate enzyme (BaeI/HaeIV/Hin4I)
- Graph creation: empty graph, single genome, two genome shared tag
- FASTA parsing: empty, single record, multi-record
- CLI: help output validation

---

## 5. Known Issues & Limitations

### 5.1 Performance: Critical Priority
- **Current**: `digest_genome_contig()` is O(N × L × P) — byte-by-byte scanning across entire sequence for each enzyme pattern
- **Fast2bRAD-M achieves ~15× speedup** via regex `find_at` + rewind for degenerate enzymes
- **Target optimizations**:
  1. **Skip-based search**: Use `memchr` for first anchor byte, then verify full pattern
  2. **Regex/aho-corasick**: For multi-pattern matching (especially degenerate enzymes)
  3. **Rayon parallel**: Per-enzyme digestion in parallel + batch processing with reusable buffers
  4. **Streaming**: For GTDB-scale (732K genomes, ~3.7 Tbp), need to avoid loading entire genomes
- **Scale estimate**: Current naive ~2–4 months for GTDB; target with optimization: ~1–2 weeks on workstation

### 5.2 Validation Results (SynTracker vs Syn2b — Important)
From benchmarking against 10 complete *C. acnes* genomes (45 pairs):
- **Tag-Jaccard** (presence/absence) correlates with SynTracker APSS (*r*=0.98) but is **sequence identity in disguise** — ~indistinguishable from Mash distance (*r*=0.997), partial correlation collapses to ~+0.06
- **Genuine synteny signal**: Tag adjacency carries real but weak structural signal: raw *r*=0.67 (Spearman 0.79), **partial *r*≈0.37** controlling for Mash (*p*≈0.015). Dynamic range on co-linear genomes is tiny (0.000–0.0025)
- **Takeaway**: "Correlates with SynTracker" is NOT evidence of capturing synteny. Must use **partial correlation** (control for ANI/Mash) or simulate known rearrangements on a single genome at fixed divergence to show tag-adjacency tracks structural change while Mash does not
- **Recommended next experiment**: Simulate inversions/translocations on a single complete genome at fixed divergence, show tag-adjacency tracks structural change while Mash does not

### 5.3 Practical Limitations
- Resolution ≈ 1–2 kb (average tag spacing ~1,660 bp; rearrangements smaller than gap are invisible)
- 32 bp tags carry limited sequence information; may map ambiguously to paralogs/repeats
- Reference-dependent ordering; de novo comparison not supported
- Sparse representation: only ~1.5% of genome sampled
- On *C. acnes*, SynTracker's own APSS is ~98% explained by sequence divergence (synteny was not an independent strain-discriminating axis for this species)

### 5.4 Code Debt
- Binary format is **v1** (magic `TGT\x01`) but README documents v2 layout with contig name table — **discrepancy needs fixing**
- `TgtWriter::write_binary()` writes v1 format (no contig name table); `TgtReader::read_binary()` reads v1 format
- Need to upgrade binary format to actually include contig name table (documented as v2 but not implemented)
- `digest_genome_contig()` does not extract reverse-complement tags (only forward strand) — need to verify if reverse tags should be extracted for degenerate enzymes
- Some integration tests use `#[allow(dead_code)]` for `contig_id` in `AnchoredContig`

---

## 6. Paper Writing Status

**Not yet started.** Need to:
1. Draft methodology section (Syn2b algorithm: TGT representation, adjacency graph, backbone extraction)
2. Benchmark results section (SynTracker comparison, partial correlation analysis, simulated rearrangement experiment)
3. Results: Scaffold validation (E. coli, ABHQ)
4. Discussion: Limitations, resolution, when Syn2b works vs doesn't
5. Introduction: 2bRAD-M background, synteny detection problem, alignment-free approaches
6. Related work: ntSynt, KmerAperture, SynTracker, Mash, MinHash

---

## 7. Next Steps (Priority Order)

### Immediate (Performance)
1. **Implement skip-based search** in `digest_genome_contig()`:
   - Use `memchr` for first anchor byte
   - For multi-anchor patterns, verify second anchor after first match
   - For degenerate enzymes, use regex `find_at` with rewind
2. **Add rayon parallelization** for multi-enzyme digestion
3. **Batch processing** with reusable buffers for GTDB-scale

### Short-term (Validation)
4. **Simulated rearrangement experiment**:
   - Take single complete genome (E. coli K-12)
   - Introduce known inversions/translocations at fixed divergence (e.g., 99% ANI)
   - Show tag-adjacency tracks structural change while Mash/ANI does not
   - This is the critical validation needed for paper claims
5. **Benchmark against larger dataset**: multiple species, not just *C. acnes*
6. **Fix binary format v2**: implement contig name table in `write_binary` and `read_binary`

### Medium-term (Features)
7. **Multi-enzyme synteny**: combine tags from multiple enzymes for denser coverage
8. **De novo comparison**: without reference dependency
9. **Visualization**: dotplot or graph visualization of synteny blocks

### Paper Writing
10. Draft full manuscript with results from simulated rearrangement experiment
11. Include all validation metrics with partial correlation controls

---

## 8. Key Code Patterns

### Adding a new enzyme
Edit `src/enzyme/enzyme.rs`:
```rust
pub enum EnzymeType { ... NewEnzyme }

impl Enzyme {
    pub fn properties(enzyme_type: EnzymeType) -> Enzyme {
        match enzyme_type {
            ...
            EnzymeType::NewEnzyme => Enzyme {
                enzyme_type: EnzymeType::NewEnzyme,
                tag_length: 27,
                patterns: vec![
                    Pattern {
                        anchors: vec![
                            Anchor { offset: 5, motif: b"GATC" },
                        ],
                        iupac: vec![
                            IupacConstraint { offset: 10, allowed: 0b0101 }, // Y = C/T
                        ],
                    },
                ],
            },
        }
    }
}
```

### IUPAC bitmask mapping
```rust
// bit0=A, bit1=T, bit2=C, bit3=G
const BASE_MASK: [u8; 256] = {
    let mut m = [0u8; 256];
    m[b'A' as usize] = 0b0001; m[b'a' as usize] = 0b0001;
    m[b'T' as usize] = 0b0010; m[b't' as usize] = 0b0010;
    m[b'C' as usize] = 0b0100; m[b'c' as usize] = 0b0100;
    m[b'G' as usize] = 0b1000; m[b'g' as usize] = 0b1000;
    m
};
// Y = C/T = 0b0110, R = A/G = 0b1001, [GAC] = 0b1101
```

### Using the library API
```rust
use bsyn::enzyme::{Enzyme, EnzymeType};
use bsyn::enzyme::digest::{digest_genome, digest_genome_contig};
use bsyn::tgt::{Tag, TgtRecord, TgtReader, TgtWriter, Strand};
use bsyn::synteny::{TagAdjacencyGraph, extract_synteny_blocks, synteny_score};

// Digest
let tags = digest_genome_contig(seq, EnzymeType::BcgI, 1, 0);

// Build record
let mut rec = TgtRecord::new("genome", 4_641_652);
for t in tags { rec.add_tag(t); }

// Write
let mut writer = TgtWriter::new(Path::new("out.tgt"))?;
writer.write_record(&rec)?;

// Graph
let mut g = TagAdjacencyGraph::new();
g.add_genome("g1", &rec);
g.build_edges();
g.simplify(2);
for path in g.linear_paths() {
    println!("score = {:.3}", synteny_score(&path, &g));
}
```

---

## 9. Dependencies

```toml
[dependencies]
clap = { version = "4", features = ["derive"] }
anyhow = "1"
serde = { version = "1", features = ["derive"] }
bincode = "1"
parking_lot = "0.12"
rayon = "1.7"

[dev-dependencies]
tempfile = "3"
```

`rayon` is listed but **not yet used** — need to add parallelization.

---

## 10. Git Repository State

- **Remote**: `git@github.com:HuangShiLab/Syn2b.git`
- **Branch**: `main`
- **History**: 4 commits from previous machine + 1 local commit (`dd010ef`) with all current work
- **All changes pushed**: Yes, `git status` is clean on main
- **No uncommitted changes**: `cargo test` passes, `cargo build --release` clean

---

## 11. Environment Notes

- **Rust toolchain**: `cargo` via `~/.cargo/bin`
- **Local workspace**: `/Users/shihuang/Downloads/Syn2b` (also `syn2b` symlink)
- **E. coli test data**: `/Users/shihuang/Downloads/Syn2b/data/` (E. coli K-12 and ABHQ draft genomes)
- **Benchmark scripts**: Python-based comparison scripts in `data/` directory (not yet organized into formal benchmark suite)

---

## 12. Critical Design Decisions

1. **Tag sequence is always 32 bytes** (zero-padded), even for shorter tag lengths (25, 27, 28 bp). This simplifies binary format and Hamming distance computation.
2. **Gaps are auto-computed** in `TgtRecord::add_tag()` — when reading, stored gaps are cross-validated against recomputed gaps.
3. **Strand is recorded per tag** but current digestion only extracts forward tags. Reverse tags may need to be added for completeness.
4. **Contig_id is u16** (max 65,535 contigs per genome), sufficient for microbial draft genomes.
5. **Binary format magic is `TGT\x01`** — need to upgrade to `TGT\x02` when adding contig name table.

---

**End of context file. Ready for Agent Swarm continuation.**
