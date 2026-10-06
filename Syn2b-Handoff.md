# Syn2b-Handoff: contig search, functional-cluster search, and structural-similarity extensions

**Date:** 2026-10-05  
**Repository:** `/Users/macstudio/Downloads/Syn2b`  
**Primary goal:** expand Syn2b from a pairwise synteny/structural-variation engine into a scalable search engine for contigs, scaffolds, and ordered functional gene clusters.  
**Non-goal:** do not fold these extensions into Syn2bANI. Syn2bANI remains focused on genome-level similarity search and ANI–structural evolutionary-mode analysis.  
**Language/runtime:** Rust crate `syn2b` / library `bsyn`.  
**Current code baseline:** `main`; core pipeline is `digest → TGT → synteny/scoring`; a basic single-reference `scaffold` command already exists.

---

## 1. Why these extensions matter

Syn2b currently represents genomes as ordered Type IIB restriction tags, or optionally FracMinHash landmarks. Its strength is not merely sequence similarity; it is the ability to preserve and search **ordered biological structure**.

The next natural abstraction is:

```text
genome
 └── scaffold / contig
      └── landmark adjacency / synteny block
           └── gene cluster / functional module
                └── encoded proteins
```

This creates three extensions:

| Priority | Extension | Working name | Main use case |
|---|---|---|---|
| **P0** | Contig/scaffold-level search and placement | **Syn2b-Contig** | MAG placement, binning refinement, chimerism detection |
| **P1** | Functional gene-cluster search | **Syn2b-Cluster** | BGC, AMR, secretion-system, metabolic-module search |
| **P2** | Protein 3D structural-similarity correlation | **Syn2b-Fold analysis** | Test whether genomic synteny predicts protein structural similarity |

P0 and P1 should become the primary Syn2b application directions. P2 is scientifically interesting but should initially be an exploratory/benchmark module, not the main software claim.

---

## 2. Current implementation baseline

### 2.1 Existing commands

Current CLI already supports:

```bash
syn2b digest   --input genome.fna --output genome.tgt
syn2b synteny  --input tgt_dir --output synteny.tsv
syn2b scaffold --reference ref.tgt --draft draft.tgt --output scaffold.agp
syn2b coverage
syn2b convert
```

The existing `scaffold` command anchors draft contigs onto a single reference using shared tags. It is useful but limited: it assumes one reference, one draft, and outputs AGP-like placement. It should become the seed for the generalized contig-search engine.

### 2.2 Relevant modules

| Path | Current role |
|---|---|
| `src/main.rs` | CLI dispatch |
| `src/enzyme/digest.rs` | Type IIB in-silico digestion |
| `src/landmark/fracminhash.rs` | FracMinHash landmark mode |
| `src/tgt/tag.rs` | Tag representation |
| `src/tgt/record.rs` | TGT genome record |
| `src/tgt/reader.rs`, `src/tgt/writer.rs` | TGT I/O |
| `src/synteny/graph.rs` | Tag adjacency graph |
| `src/synteny/scoring.rs` | Structural synteny scoring |
| `src/main.rs`, scaffold section | Basic single-reference contig anchoring |

### 2.3 Existing strengths

- TGT model preserves ordered landmarks.
- Supports 2bRAD enzyme tags and FracMinHash landmarks.
- Synteny scoring already handles canonical orientation, shared-tag restriction, contig boundaries, overlapping cut-site collapse, and circular-origin normalization.
- Basic AGP scaffold output exists.
- Rust implementation is suitable for large-scale indexing/search.

---

## 3. Syn2b-Contig: contig and scaffold search

### 3.1 Goal

Given one or more query contigs/scaffolds, rapidly search a reference genome or MAG database and report:

```text
best matching genome
best matching scaffold/contig
placement orientation
supported landmark adjacencies
genome coverage
species-cluster assignment if available
chimerism / multi-source warning
```

This is not a full aligner. It is a fast placement and search layer.

### 3.2 Primary use cases

1. **MAG placement**  
   Place a newly binned MAG or unbinned contig into an existing species/genome cluster.

2. **Contig-level search**  
   Find reference genomes containing conserved ordered structure.

3. **Binning refinement**  
   Add synteny evidence to coverage/composition-based binning.

4. **Chimerism detection**  
   Detect a contig whose different regions match unrelated genome clusters.

5. **MAG gap filling**  
   Identify whether a contig contributes missing syntenic regions to an existing MAG.

---

## 3.3 New CLI design

### Build database index

```bash
syn2b contig-index \
  --genome-tgt-dir reference_tgt_dir \
  --output reference.s2cidx
```

Optional genome metadata:

```bash
syn2b contig-index \
  --genome-tgt-dir reference_tgt_dir \
  --genome-metadata genomes.tsv \
  --output reference.s2cidx
```

`genomes.tsv` columns:

```text
genome_id    species_cluster    taxonomy    fasta_path    is_representative
```

### Search query contigs

```bash
syn2b search-contig \
  --index reference.s2cidx \
  --query query_contigs.tgt \
  --output placements.tsv \
  --min-tags 5 \
  --min-adjacencies 3 \
  --max-targets 25
```

### Detailed placement report

```bash
syn2b place-contig \
  --index reference.s2cidx \
  --query query_contigs.tgt \
  --output placements.tsv \
  --details placement_details.tsv
```

### Chimera scan

```bash
syn2b scan-chimera \
  --index reference.s2cidx \
  --query query_contigs.tgt \
  --output chimerism.tsv \
  --window 50kb \
  --min-window-tags 10
```

---

## 3.4 Index design

Create a new module:

```text
src/search/mod.rs
src/search/index.rs
src/search/query.rs
src/search/chaining.rs
src/search/scoring.rs
src/search/chimera.rs
src/search/io.rs
```

Recommended index name: `.s2cidx`.

The index should store:

1. Genome metadata.
2. Contig metadata.
3. Tag dictionary.
4. Tag-to-genome/contig postings.
5. Adjacency dictionary.
6. Adjacency-to-genome postings.
7. Genome-level summary statistics.
8. Optional species-cluster mapping.

### 3.4.1 Tag postings

For every landmark:

```text
tag_id → [(genome_id, contig_id, position, strand)]
```

Tags should be canonicalized so inverted regions retain matching landmarks.

### 3.4.2 Adjacency postings

For every adjacent landmark pair on the same contig:

```text
(tag_id_a, tag_id_b) → [(genome_id, contig_id, position_a, position_b, strand)]
```

Use canonical unordered adjacency:

```text
key = min(tag_id_a, tag_id_b), max(tag_id_a, tag_id_b)
```

For placement, orientation can be recovered separately.

### 3.4.3 Genome-level summaries

Store per genome:

```text
genome_id
total_tags
total_adjacencies
total_length
contig_count
n50
species_cluster
taxonomy
```

These summaries are needed for normalization and confidence scoring.

---

## 3.5 Search algorithm

### Step 1: seed lookup

For a query contig:

```text
shared tags = query tags ∩ reference tag postings
```

Generate candidate genomes using shared tags and shared adjacencies.

Recommended candidate score:

```text
seed_score =
  weighted_shared_tags
+ adjacency_support
+ unique_genome_support
```

### Step 2: adjacency chaining

For each candidate genome, collect matched tag pairs and chain them by reference coordinate.

Output chained blocks:

```text
query_start    query_end    ref_start    ref_end    orientation    n_tags    n_adjacencies
```

This is analogous to synteny-block chaining but generalized to query contigs.

### Step 3: placement scoring

For each query contig against candidate genome:

```text
PlacementScore =
  w1 × shared_tag_fraction
+ w2 × conserved_adjacency_fraction
+ w3 × chained_block_coverage
+ w4 × orientation_consistency
+ w5 × uniqueness_penalty
+ w6 × contig_interior_support
```

Recommended initial weights:

```text
w1 = 0.15
w2 = 0.30
w3 = 0.30
w4 = 0.10
w5 = 0.10
w6 = 0.05
```

These should later be calibrated on simulated and real datasets.

### Step 4: confidence calibration

Use simulated contigs to calibrate score to posterior probability:

```text
P(correct placement | PlacementScore, query_length, reference_quality)
```

Initial bins:

| Query length | Minimum tags | Minimum adjacencies |
|---:|---:|---:|
| ≥500 kb | ≥25 | ≥15 |
| 100–500 kb | ≥15 | ≥8 |
| 50–100 kb | ≥10 | ≥5 |
| 10–50 kb | ≥5 | ≥2 |
| <10 kb | exploratory only |

### Step 5: chimera segmentation

For long contigs, split into sliding windows.

For each window:

```text
window_id
start
end
best_genome
best_species_cluster
second_best_cluster
support_score
```

Flag chimera if:

```text
left_windows dominant cluster A
right_windows dominant cluster B
A ≠ B
both supports exceed threshold
junction near window transition
```

Output:

```text
query_id
breakpoint_start
breakpoint_end
left_cluster
right_cluster
left_support
right_support
confidence
```

---

## 3.6 Output formats

### `placements.tsv`

```text
query_id
query_length
query_tags
query_adjacencies
best_genome
best_species_cluster
best_contig
orientation
shared_tags
shared_adjacencies
chained_blocks
query_coverage
reference_support_fraction
placement_score
confidence
secondary_genome
secondary_score
chimera_warning
```

### `placement_details.tsv`

One row per chained block:

```text
query_id
query_block_id
query_start
query_end
reference_genome
reference_contig
reference_start
reference_end
orientation
n_tags
n_adjacencies
block_score
```

### `chimerism.tsv`

```text
query_id
segment
start
end
best_cluster
best_genome
score
second_cluster
second_score
cluster_transition
chimera_confidence
```

---

## 3.7 Benchmark plan

### Benchmark A: simulated fragmentation

Use complete isolate genomes from GTDB, HROM, and OAPGC.

Fragment each genome into contigs of:

```text
500 kb
250 kb
100 kb
50 kb
25 kb
10 kb
5 kb
```

For each genome:

```text
all fragments → known placement truth
```

Measure:

```text
placement precision
placement recall
F1
top-1 accuracy
top-5 accuracy
species-cluster accuracy
orientation accuracy
coverage error
runtime
peak memory
```

### Benchmark B: real MAG contigs

Use HROM and OAPGC representative MAGs.

For each MAG:

1. split into contigs if multi-contig;
2. leave-one-contig-out;
3. search against the full reference database excluding the source MAG.

Measure whether the missing contig maps back to the correct species cluster.

### Benchmark C: chimera simulation

Create artificial chimeric contigs:

```text
left 50–500 kb from genome A
right 50–500 kb from genome B
```

Use same-genus pairs, same-family pairs, and unrelated pairs.

Measure:

```text
chimera sensitivity
chimera specificity
breakpoint localization error
false-positive chimera rate
```

### Benchmark D: runtime and scale

Build indexes at multiple scales:

```text
1,000 genomes
10,000 genomes
50,000 genomes
100,000 genomes
```

Report:

```text
index build time
index size
queries/sec
peak RSS
disk usage
```

### Baselines

Compare against:

```text
minimap2 -x asm20
MUMmer/nucmer
skani search
FastANI
BLASTn
mash sketch/search
sourmash gather
```

Syn2b does not need to beat minimap2 in nucleotide-level sensitivity. It should win or remain competitive in:

```text
search speed
index reusability
low-memory query
structural context
chimerism detection
```

---

## 3.8 Acceptance criteria for Syn2b-Contig

Minimum acceptable:

```text
≥95% top-1 species-cluster accuracy for query contigs ≥100 kb
≥90% top-1 genome accuracy for query contigs ≥250 kb
runtime faster than minimap2 asm20 at ≥10,000-genome database scale
index reusable without rebuilding
no unmapped genome IDs
deterministic output
```

Strong target:

```text
≥98% species-cluster accuracy for contigs ≥100 kb
chimera sensitivity ≥90% for ≥250 kb artificial chimeras
chimera false-positive rate ≤5%
index search ≥10× faster than minimap2 asm20 for large databases
```

---

# 4. Syn2b-Cluster: functional gene-cluster search

## 4.1 Goal

Given a query gene cluster, search a genome/database for conserved or evolutionarily related gene clusters.

Examples:

```text
BGCs
AMR clusters
virulence loci
type III/IV/VII secretion systems
CRISPR-associated modules
conjugation modules
metabolic operons
stress-response clusters
```

The search should use both:

```text
functional gene composition
ordered gene adjacency / synteny
```

---

## 4.2 New CLI design

### Annotate or ingest annotations

Syn2b should not become a full gene annotator. It should ingest existing annotations.

```bash
syn2b cluster-ingest \
  --genome-id genome_id \
  --gff genes.gff \
  --eggnog proteins.emapper.annotations \
  --output genome.gclusters.tsv
```

Alternative:

```bash
syn2b cluster-ingest \
  --genome-id genome_id \
  --features features.tsv \
  --output genome.gclusters.tsv
```

`features.tsv` format:

```text
gene_id    contig    start    end    strand    product    ko    pfam    cog    eggnog    smash_tool
```

### Define gene clusters

```bash
syn2b cluster-build \
  --features genome.gclusters.tsv \
  --output genome.s2cl \
  --max-gap 5000 \
  --min-genes 3
```

### Build cluster database

```bash
syn2b cluster-index \
  --cluster-files '*.s2cl' \
  --output clusters.s2cldb
```

### Search

```bash
syn2b search-cluster \
  --query query_cluster.tsv \
  --index clusters.s2cldb \
  --output cluster_hits.tsv \
  --min-functional-jaccard 0.25 \
  --min-synteny 0.10
```

---

## 4.3 Gene-cluster model

Each gene is represented by several feature tokens:

```text
KO tokens
Pfam tokens
eggNOG OG tokens
COG tokens
antiSMASH domain tokens
product keyword tokens
```

Example:

```text
gene_1:
  KO:K03515
  PFAM:MutS_I
  COG:COG1193
  PRODUCT:mismatch_repair
```

Each cluster is then:

```text
gene set
ordered gene adjacency set
copy-number vector
contig coordinates
```

---

## 4.4 Scoring

### Functional similarity

```text
FunctionalJaccard =
  |query_genes ∩ target_genes| /
  |query_genes ∪ target_genes|
```

Better weighted version:

```text
WeightedFunctionalSimilarity =
  sum(weight_i for shared gene_i) /
  sum(weight_i for union gene_i)
```

Suggested weights:

```text
KO / eggNOG orthologue: 3
Pfam domain: 2
product keyword: 1
antiSMASH domain: 3
transporter class: 2
regulator: 1
```

### Synteny similarity

```text
AdjacencyScore =
  |ordered query adjacencies ∩ ordered target adjacencies| /
  |ordered query adjacencies ∪ ordered target adjacencies|
```

Handle missing genes with gene-aware adjacency projection.

For example:

```text
query: A → B → C → D
target: A → C → D
```

After deleting missing B:

```text
projected query: A → C → D
```

This allows conserved order despite gene loss.

### Copy-number score

```text
CopyScore =
  1 - normalized_copy_number_distance
```

Useful for:

```text
duplicated KS domains
repeated transporter genes
expanded regulator families
```

### Combined cluster score

```text
ClusterScore =
  0.45 × FunctionalJaccard
+ 0.30 × AdjacencyScore
+ 0.15 × ProjectedOrderScore
+ 0.10 × CopyScore
```

Output both individual components and combined score.

---

## 4.5 Gene-cluster index

New modules:

```text
src/cluster/mod.rs
src/cluster/features.rs
src/cluster/ingest.rs
src/cluster/model.rs
src/cluster/signature.rs
src/cluster/index.rs
src/cluster/search.rs
src/cluster/scoring.rs
src/io/gff.rs
src/io/features.rs
```

Recommended index name: `.s2cldb`.

Store:

```text
gene token dictionary
cluster signatures
gene-to-cluster postings
adjacency-to-cluster postings
genome-to-cluster mapping
optional taxonomy
optional BGC type
```

Use approximate indexing for large databases:

```text
MinHash for functional gene sets
LSH or blocked postings for adjacency signatures
```

---

## 4.6 Benchmark plan

### Benchmark A: MIBiG retrieval

Use MIBiG clusters as truth.

For each query BGC:

```text
leave out its MIBiG entry
search against genome-derived clusters
check whether top hits recover the correct BGC family
```

Metrics:

```text
top-1 family accuracy
top-5 family accuracy
precision@10
recall@10
mean reciprocal rank
```

### Benchmark B: synthetic rearrangement

Take real gene clusters and simulate:

```text
gene loss
gene duplication
inversion
translocation
operon splitting
contig truncation
```

For each event level, measure:

```text
cluster score retention
gene-order recovery
adjacency recovery
false cluster splitting
```

### Benchmark C: real AMR/virulence clusters

Use known clusters:

```text
beta-lactamase loci
vancomycin resistance clusters
macrolide resistance clusters
type IV secretion systems
type VI secretion systems
```

Validate that Syn2b-Cluster recovers known homologues and distinguishes structurally different variants.

### Benchmark D: runtime

Build indexes at:

```text
1,000 genomes
10,000 genomes
50,000 genomes
100,000 genomes
```

Report:

```text
queries/sec
index size
peak memory
build time
```

---

## 4.7 Baselines

Compare against:

```text
antiSMASH KnownClusterBlast
MultiGeneBlast
BLASTP-based gene order search
EggNOG profile search
DRAM/AMRFinder profile search
simple KO Jaccard
simple Pfam Jaccard
```

Syn2b-Cluster should show advantage when:

```text
gene order is conserved
sequence identity is low
genes are rearranged
clusters are fragmented
functional domains are conserved but gene calling differs
```

---

## 4.8 Acceptance criteria

Minimum:

```text
top-10 BGC family accuracy ≥80%
runtime ≥10× faster than naive BLASTP gene-order search
supports gene loss and copy-number variation
deterministic output
```

Strong target:

```text
top-5 family accuracy ≥90%
synthetic gene-order recovery AUC ≥0.95
AMR/BGC validation on independent datasets
index query <1 s for 10,000-genome database
```

---

# 5. Protein 3D structural-similarity correlation

This is scientifically interesting but should be P2 exploratory.

## 5.1 Biological question

Does Syn2b-detected gene-cluster or genome synteny predict protein 3D structural similarity?

This must be separated from sequence identity.

---

## 5.2 Three possible comparisons

### A. Orthologous protein comparison

For each orthologous protein pair:

```text
protein A in genome 1
orthologous protein A in genome 2
```

Compare:

```text
Syn2b genome/cluster similarity
Foldseek/TM-align structural similarity
```

### B. Gene-cluster protein-set comparison

Compare the whole protein set encoded by two clusters.

Use:

```text
Foldseek all-vs-all
best-hit structural score
domain architecture similarity
structural Jaccard
```

Compare against:

```text
Syn2b cluster synteny score
functional gene similarity
gene-order conservation
```

### C. Structural consequence of SV

For SV breakpoints intersecting genes, test whether events produce:

```text
gene truncation
gene fusion
domain loss
domain insertion
copy-number change
regulatory inversion
```

This is the most functionally interesting but hardest layer.

---

## 5.3 Recommended tools

Large-scale structural search:

```bash
foldseek easy-search query.fa target.fa out.m8 out_tmp
```

High-confidence validation:

```bash
TM-align
US-align
```

Protein structure prediction:

```text
ESMFold for large-scale screening
ColabFold / AlphaFold2 for selected proteins
OmegaFold for metagenomic sequences
```

---

## 5.4 Analysis design

### Step 1: select clusters

Use clusters with:

```text
≥5 genomes
complete or high-quality MAGs
known or predicted gene annotation
diverse ANI ranges
```

### Step 2: identify orthologous proteins

Use:

```text
eggNOG best orthologue
KO groups
Pfam domains
BLASTP reciprocal best hits
```

### Step 3: compute structural similarity

For each orthologous protein pair:

```text
Foldseek alignment
TM-align for validation
```

Metrics:

```text
TM-score
aligned length
RMSD
sequence identity
evalue
qcov / tcov
```

### Step 4: compute genomic similarity

For the corresponding genome pair:

```text
ANI
Syn2b breakpoints
Syn2b synteny score
gene-cluster synteny score
breakpoint density
shared adjacencies
```

### Step 5: correlation

Primary test:

```text
Spearman correlation between
Syn2b cluster synteny score
and Foldseek structural similarity
```

Critical control:

```text
partial out protein sequence identity
```

Without this control, the result may simply reflect sequence conservation.

### Step 6: stratify by divergence

Analyze separately:

```text
ANI ≥99.9%
99–99.9%
97–99%
95–97%
<95%
```

This shows whether structural and genomic similarity decouple at increasing divergence.

---

## 5.5 Expected outcomes

Possible findings:

### A. Strong local correlation

```text
gene-cluster synteny predicts protein structural similarity
```

This would support Syn2b-Cluster as a functional module search tool.

### B. Weak global correlation

```text
genome rearrangement does not strongly predict protein fold similarity
```

This is biologically plausible and still valuable.

### C. Sequence identity dominates

```text
after controlling sequence identity,
Syn2b similarity adds little structural signal
```

Then report honestly: Syn2b is a genomic structural search tool, not a protein structure predictor.

### D. Structural convergence

Some protein pairs show:

```text
low sequence identity
low Syn2b similarity
high Foldseek similarity
```

These may represent fold reuse, domain convergence, or ancient horizontal acquisition.

---

## 5.6 Acceptance criteria

Exploratory analysis:

```text
≥100 gene clusters
≥1,000 orthologous protein pairs
Foldseek and TM-align validation on subset
sequence-identity-controlled analysis
```

Strong paper-level claim:

```text
≥1,000 clusters
≥10,000 orthologous protein pairs
partial correlation significant after sequence-identity control
validated subset using TM-align/US-align
clear functional examples
```

---

# 6. Development roadmap

## Phase 0: clean baseline

### Tasks

1. Refactor existing scaffold logic into `src/search/`.
2. Keep old `syn2b scaffold` CLI as a compatibility wrapper.
3. Add integration tests for:
   - TGT loading;
   - contig metadata handling;
   - shared-tag lookup;
   - adjacency canonicalization;
   - orientation handling;
   - AGP output.
4. Add benchmark harness.

### Deliverables

```text
src/search/ skeleton
regression tests
benchmark runner
baseline performance report
```

### Exit criteria

Existing `syn2b scaffold` results unchanged after refactor.

---

## Phase 1: Syn2b-Contig MVP

### Tasks

1. Implement `.s2cidx` index writer/reader.
2. Implement shared-tag candidate search.
3. Implement adjacency postings.
4. Implement greedy chaining.
5. Implement placement scoring.
6. Implement `search-contig`.
7. Implement `place-contig`.
8. Add simulated fragmentation benchmark.
9. Add minimap2 baseline comparison.

### Deliverables

```text
syn2b contig-index
syn2b search-contig
syn2b place-contig
benchmark_contig_fragmentation.tsv
benchmark_contig_runtime.tsv
```

### Exit criteria

```text
≥95% species-cluster accuracy for ≥100 kb contigs
index reusable
runtime benchmark complete
```

---

## Phase 2: Syn2b-Contig chimera detection

### Tasks

1. Implement sliding-window search.
2. Implement dominant-cluster transition detection.
3. Calibrate chimera confidence.
4. Simulate chimeric contigs.
5. Benchmark against real MAG contamination cases.

### Deliverables

```text
syn2b scan-chimera
chimera benchmark report
```

### Exit criteria

```text
chimera sensitivity ≥90% for ≥250 kb chimeras
false-positive rate ≤5%
```

---

## Phase 3: Syn2b-Cluster MVP

### Tasks

1. Implement annotation ingestion.
2. Implement gene-token signatures.
3. Implement cluster model.
4. Implement functional Jaccard.
5. Implement adjacency-aware score.
6. Implement `.s2cldb` index.
7. Implement `search-cluster`.
8. Benchmark against MIBiG and real AMR/BGC clusters.

### Deliverables

```text
syn2b cluster-ingest
syn2b cluster-build
syn2b cluster-index
syn2b search-cluster
benchmark_cluster_retrieval.tsv
```

### Exit criteria

```text
top-10 MIBiG family accuracy ≥80%
runtime benchmark complete
```

---

## Phase 4: Syn2b-Cluster maturation

### Tasks

1. Add gene-aware adjacency projection.
2. Add copy-number score.
3. Add cluster completeness score.
4. Add inversion/translocation classification.
5. Add taxonomy-aware filtering.
6. Benchmark independent AMR/virulence datasets.

### Exit criteria

```text
top-5 family accuracy ≥90%
robust under gene loss/duplication/rearrangement
```

---

## Phase 5: protein structure correlation

### Tasks

1. Select 100–1,000 gene clusters.
2. Map proteins to structures using Foldseek/ESMFold.
3. Validate subset with TM-align/US-align.
4. Compare structural similarity with Syn2b cluster similarity.
5. Control for protein sequence identity.
6. Stratify by ANI/genome divergence.

### Deliverables

```text
protein_structure_correlation.tsv
partial_correlation_report.md
structural_convergence_examples.tsv
```

### Exit criteria

Clear reporting even if correlation is weak.

---

# 7. Proposed CLI architecture

## Current

```bash
syn2b digest
syn2b synteny
syn2b scaffold
syn2b coverage
syn2b convert
```

## Proposed

```bash
# Existing
syn2b digest
syn2b synteny
syn2b scaffold

# Contig extension
syn2b contig-index
syn2b search-contig
syn2b place-contig
syn2b scan-chimera

# Gene-cluster extension
syn2b cluster-ingest
syn2b cluster-build
syn2b cluster-index
syn2b search-cluster

# Optional analysis
syn2b fold-ingest
syn2b fold-correlate
```

Keep old `scaffold` command for compatibility, but internally route it through the new search/scoring layer after refactor.

---

# 8. File formats

## 8.1 Genome metadata

```text
genome_id
species_cluster
taxonomy
fasta_path
total_length
contig_count
n50
is_representative
```

## 8.2 Contig metadata

```text
genome_id
contig_id
contig_name
length
n50
is_circular
```

## 8.3 Feature table

```text
gene_id
genome_id
contig_id
start
end
strand
product
ko
pfam
cog
eggnog
smash_tool
cluster_id
```

## 8.4 Placement report

See section 3.6.

## 8.5 Cluster report

```text
query_cluster_id
target_cluster_id
target_genome
target_contig
functional_jaccard
weighted_functional_similarity
adjacency_score
projected_order_score
copy_number_score
cluster_score
missing_query_genes
extra_target_genes
inversion_flag
split_cluster_flag
```

## 8.6 Protein structure comparison

```text
protein_pair_id
genome_A
protein_A
genome_B
protein_B
orthology_method
sequence_identity
foldseek_align_length
foldseek_rmsd
foldseek_prob
tm_score
us_align_rmsd
cluster_syn2b_score
ANI
```

---

# 9. Statistical plan

## 9.1 Contig placement

Use:

```text
precision
recall
F1
top-k accuracy
bootstrap confidence intervals
```

Bootstrap genomes, not individual fragments, to avoid pseudo-replication.

---

## 9.2 Chimerism

Use:

```text
sensitivity
specificity
precision
F1
breakpoint localization error
```

Report stratified by:

```text
contig length
genus divergence
repeat density
N50
```

---

## 9.3 Gene-cluster search

Use:

```text
precision@1
precision@5
precision@10
recall@10
mean reciprocal rank
family-level accuracy
```

Use held-out MIBiG families and independent AMR datasets.

---

## 9.4 Protein structure correlation

Primary analysis:

```text
Spearman correlation between
Syn2b cluster similarity
and Foldseek structural similarity
```

Critical control:

```text
partial Spearman correlation controlling protein sequence identity
```

Also report:

```text
binomial sign test
stratified ANI analysis
TM-align validation subset
```

---

# 10. Engineering requirements

## Rust implementation

Use existing style:

```text
src/search/
src/cluster/
src/io/
```

Requirements:

- deterministic output;
- no panics in library code;
- explicit errors with `anyhow`;
- stable TSV headers;
- optional binary index;
- streaming query support;
- parallel search with Rayon;
- memory ceilings for large databases;
- no implicit annotation prediction inside core Syn2b.

---

## Testing

Add unit tests for:

```text
tag canonicalization
adjacency canonicalization
orientation handling
contig-boundary handling
gene-token parsing
GFF ingestion
KO/Pfam token extraction
cluster signature construction
chaining correctness
chimera segmentation
```

Add integration tests for:

```text
digest → contig-index → search-contig
features → cluster-build → cluster-index → search-cluster
fragmented genome placement
simulated chimera detection
simulated gene rearrangement
```

---

# 11. Benchmark datasets

## GTDB

Use GTDB-R207 multi-genome clusters for:

```text
contig placement simulation
genome-index scaling
quality-stratified benchmarking
```

## HROM

Use HROM species clusters for:

```text
real MAG fragmentation
strain-level placement
repair-gene candidates
```

## OAPGC

Use OAPGC oral/airway representatives for:

```text
site-specific placement
fragmented MAG simulation
functional-cluster search validation
```

## MIBiG / antiSMASH

Use MIBiG and antiSMASH clusters for:

```text
BGC retrieval
functional cluster benchmark
synthetic rearrangement
```

## AMRFinder / CARD / ResFinder

Use AMR loci for:

```text
AMR cluster retrieval
structural variant classification
independent validation
```

---

# 12. Priority sequence

## Immediate next step

Start with **Syn2b-Contig**.

Reasons:

1. closest to existing scaffold implementation;
2. clear benchmark truth;
3. immediate biological use in MAG workflows;
4. can use GTDB/HROM/OAPGC datasets already available.

---

## Second implementation

Implement **Syn2b-Cluster**.

Reasons:

1. expands Syn2b into functional search;
2. creates a natural companion to eggNOG/antiSMASH/AMRFinder;
3. supports BGC and AMR applications;
4. provides the data needed for protein-structure correlation later.

---

## Third implementation

Protein 3D correlation should be done after Syn2b-Cluster exists.

Reasons:

1. requires cluster definitions;
2. requires ortholog mapping;
3. requires large-scale protein structures/search;
4. should not be central to the first Syn2b extension paper.

---

# 13. Suggested paper framing

## Syn2b main paper

Keep primary claim:

```text
alignment-free ordered structural comparison of microbial genomes
```

Add three application demonstrations:

```text
1. GTDB/HROM/OAPGC within-species structural census
2. contig placement and chimerism detection
3. functional gene-cluster search
```

Do not place protein-structure correlation as a main claim unless the analysis is strong.

---

## Follow-up paper 1: Syn2b-Contig

Suggested title direction:

```text
Rapid contig placement and chimerism detection using ordered restriction-landmark synteny
```

Core claim:

```text
Syn2b-Contig places fragmented MAG contigs and detects chimeric contigs faster than alignment-based search while preserving structural context.
```

---

## Follow-up paper 2: Syn2b-Cluster

Suggested title direction:

```text
Syn2b-Cluster: alignment-free search for conserved and rearranged functional gene clusters
```

Core claim:

```text
Syn2b-Cluster detects conserved, rearranged, fragmented, and divergent functional gene clusters using ordered gene adjacency rather than sequence similarity alone.
```

---

# 14. Immediate handoff actions

For the separate Syn2b development conversation, start with these tasks:

## Task 1

Create module skeleton:

```text
src/search/mod.rs
src/search/index.rs
src/search/query.rs
src/search/chaining.rs
src/search/scoring.rs
src/search/chimera.rs
```

## Task 2

Move/refactor existing `run_scaffold` logic into `src/search/`.

## Task 3

Preserve old CLI:

```bash
syn2b scaffold
```

but add:

```bash
syn2b contig-index
syn2b search-contig
```

## Task 4

Implement `.s2cidx` version 1.

Minimal version can store:

```text
genome metadata
tag dictionary
tag postings
adjacency postings
```

## Task 5

Implement first benchmark:

```text
fragment 100 complete genomes at 10/25/50/100/250/500 kb
search fragments against full database
report top-1/top-5 accuracy and runtime
```

## Task 6

Compare against:

```bash
minimap2 -x asm20
skani search
mash dist
```

## Task 7

Only after Syn2b-Contig MVP passes benchmark, begin Syn2b-Cluster.

---

# 15. Current census context for the Syn2b paper

This handoff focuses on extensions, but the current census status should be preserved.

GTDB effective structural census:

```text
effective genomes:        261,737
effective clusters:        21,973
effective pairs:      653,584,719
reported pairs:       646,012,820
pairs with ≥2 breakpoints: 502,819,691
fraction ≥2:                   77.83%
```

QC exception:

```text
Escherichia coli
expected: 338,793,465
observed: 331,221,566
shortfall:  7,571,899
```

OAPGC:

```text
ANI clusters complete: 1,488 / 1,488
structural pairs: 30,560,247
oral–airway joint analysis complete
```

Main OAPGC result:

```text
200 eligible species
oral–airway pairs show SNP-distance enrichment
no consistent global breakpoint excess
```

HROM:

```text
62,872,901 pairs complete
27.89% with ≥2 breakpoints
ANI matched for all pairs
```

---

# 16. Final recommendation

Implement in this order:

```text
1. Syn2b-Contig MVP
2. Syn2b-Contig benchmark
3. Syn2b-Contig chimera detection
4. Syn2b-Cluster MVP
5. Syn2b-Cluster benchmark
6. Optional protein 3D correlation
```

The main Syn2b paper can present the first three as application extensions. The protein-structure analysis should be exploratory until strong sequence-controlled evidence is available.
