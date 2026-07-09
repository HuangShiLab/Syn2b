# Syn2b Optimization Plan

## Overview
Four high-priority tasks for Syn2b optimization, to be executed in parallel by specialized sub-agents.

## Stage 1 — Parallel Implementation (3 Workers)

### Worker A: Performance_Engineer
**Scope**: `src/enzyme/digest.rs`, `src/main.rs`, `Cargo.toml`
**Goal**: Accelerate `digest_genome_contig()` from O(N×L×P) byte-by-byte scanning to skip-based search, and add rayon parallelization for multi-enzyme digestion.

**Tasks**:
1. Add `memchr` dependency to `Cargo.toml`.
2. Rewrite `digest_genome_contig()` to use `memchr::memchr` for the **first anchor byte** of each pattern. When a candidate position is found, verify the full anchor(s) and IUPAC constraints within the tag-length window. This skips large stretches of non-matching sequence.
3. For degenerate enzymes (BaeI, HaeIV, Hin4I), consider a regex-based or bitmask-accelerated approach. The Fast2bRAD-M approach uses `regex find_at` with rewind. Since our patterns are simpler (fixed offsets within a small window), memchr + manual verification is likely sufficient and avoids regex overhead for the common case.
4. Add rayon parallelization in `run_digest()`: parallelize the per-enzyme loop over `rec.sequence` using `rayon::join` or `par_iter`. When `--enzymes all` is used, each enzyme can be digested independently and results merged.
5. Ensure all existing digestion tests still pass. Do NOT change public API signatures unless necessary.

**Validation**:
- `cargo test` passes
- `cargo test --release` passes
- `cargo build --release` compiles with 0 warnings

### Worker B: Format_Engineer
**Scope**: `src/tgt/writer.rs`, `src/tgt/reader.rs`, `tests/integration_tests.rs`, `README.md`
**Goal**: Implement the documented TGT binary v2 format.

**Current state**: README documents v2 (`TGT\x02`, 48-byte header with `contig_count`, contig name table), but code writes v1 (`TGT\x01`, 32-byte header, no contig names).

**Tasks**:
1. Update `TgtWriter::write_binary()` to write v2 format:
   - Magic: `TGT\x02`
   - Version: `2`
   - Header: 48 bytes (add `contig_count: u16` at bytes 22..24)
   - After gap table, write contig name table: for each contig, `u16 name_len + name bytes`
   - Ensure backward compatibility is NOT required (this is an internal format; we bump version).
2. Update `TgtReader::read_binary()` to read v2 format:
   - Accept magic `TGT\x02` and version `2`
   - Read `contig_count` from header
   - After gap table, read contig name table and populate `record.contig_names`
   - Remove support for reading v1 (or keep it if trivial, but v1 was never used in production).
3. Update integration tests in `tests/integration_tests.rs` and unit tests in `writer.rs`/`reader.rs` to expect v2 magic and verify contig name round-trip.
4. Update README.md binary format section to match the actual implementation (fix any discrepancies).

**Validation**:
- `cargo test` passes
- `cargo test --release` passes
- Binary round-trip test verifies contig names are preserved

### Worker C: Validation_Engineer
**Scope**: New file `scripts/simulate_rearrangement.py`
**Goal**: Create a simulated rearrangement experiment to support the paper's core claim.

**Context**: The current validation against SynTracker on *C. acnes* showed that tag-Jaccard is "sequence identity in disguise" (r=0.997 with Mash). The genuine synteny signal is weak. The paper needs a stronger claim: "tag adjacency tracks structural change while Mash does not."

**Tasks**:
1. Load E. coli K-12 complete genome (NC_000913) from a FASTA file (assume path `data/e_coli_k12.fasta` or accept as CLI arg).
2. Generate a set of **derived genomes** at a fixed divergence level (e.g., ~99% ANI = ~1% substitutions) using a simple neutral mutator:
   - Randomly substitute bases at a given rate (e.g., 0.01).
   - Introduce **known structural variations** on top of the mutated genome:
     - Inversions of varying sizes (e.g., 10kb, 50kb, 100kb, 500kb, 1Mb)
     - Translocations (swap two segments)
   - Create a control set with **only substitutions** (no rearrangements).
3. Run Syn2b `digest` on each derived genome to produce `.tgt` files.
4. Compute two metrics between original and each derived genome:
   - **Mash distance** (or simple ANI proxy: fraction of matching k-mers or Jaccard of 21-mers)
   - **Syn2b tag-adjacency metrics**: `adjacency_jaccard`, `kendall_tag_order`, `breakpoint_count`
5. Generate a figure (matplotlib) showing:
   - X-axis: structural variation type/size
   - Y-axis: metric value (Mash distance vs Syn2b adjacency score)
   - The key result: Mash distance should be nearly identical for "substitutions only" vs "substitutions + rearrangements" (because Mash is insensitive to rearrangements at fixed SNP rate), while Syn2b adjacency metrics should clearly separate them.
6. Also output a CSV with raw results for statistical analysis.

**Validation**:
- Script runs end-to-end on E. coli K-12
- Figure clearly shows separation between rearranged and non-rearranged genomes for Syn2b metrics but not for Mash
- CSV contains all pairwise comparisons

## Merge Order
1. Worker B (format v2) — self-contained, low risk
2. Worker A (performance) — touches main digestion path
3. Worker C (validation script) — independent, just needs final binary

## Final Verification
- `cargo test` passes (all 99+ tests)
- `cargo build --release` clean
- `cargo clippy` (if available) with no new warnings
- Simulated rearrangement script produces expected figure
