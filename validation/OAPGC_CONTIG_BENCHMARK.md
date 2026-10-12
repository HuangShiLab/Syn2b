# Syn2b-Contig real-data validation on OAPGC

Date: 2026-10-07  
Dataset: OAPGC TGT census, 98,229 genomes  
Metadata: `/Users/macstudio/Downloads/OAPGC/oapgc_genome_metadata.tsv`  
TGT directory: `/Users/macstudio/Downloads/OAPGC/census/tgt`  
Index format: `.s2cidx` v1  
Landmarks: OAPGC multi-enzyme Type IIB TGTs  
Selection seed: `20261006`

## Method

For each run, deterministic sampling selected high-quality species clusters
(completeness >= 95%, contamination <= 2%, N50 >= 50 kb). One genome per
selected cluster was chosen as the fragment source. Landmark records were cut
inside each source TGT contig at 500/250/100/50/25/10 kb windows and rewritten
as query TGTs.

The reference database consisted of:

- selected target genomes for the source-present fragmentation experiment;
- cluster mates from those selected clusters;
- a deterministic random OAPGC background.

The source-present experiment corresponds to the handoff's Benchmark A:
the source genome is present in the reference database and fragmented query
placement is evaluated against known truth. It is not a leave-one-genome-out
test.

## Results

### 1k-scale run

- Reference genomes: 2,861
- Source clusters: 30
- Queries: 3,958
- Index size: 1.232 GB
- Index build: 28.39 s
- Query time: 24.92 s
- Throughput: 158.8 queries/s

| Level (bp) | Queries | Top-1 species cluster | Top-5 species cluster | Top-1 genome | Forward orientation | Mean query coverage |
|---:|---:|---:|---:|---:|---:|---:|
| ALL | 3,958 | 99.39% | 99.65% | 93.81% | 97.32% | 0.9611 |
| 500,000 | 506 | 99.80% | 100.00% | 97.43% | 99.01% | 0.9993 |
| 250,000 | 523 | 99.81% | 100.00% | 97.51% | 99.04% | 0.9988 |
| 100,000 | 578 | 99.83% | 100.00% | 97.23% | 98.62% | 0.9950 |
| 50,000 | 659 | 99.85% | 100.00% | 96.97% | 98.79% | 0.9855 |
| 25,000 | 778 | 99.87% | 100.00% | 94.09% | 97.81% | 0.9619 |
| 10,000 | 914 | 97.92% | 98.47% | 85.01% | 93.11% | 0.8765 |

### 5k-scale run

- Reference genomes: 7,298
- Source clusters: 50
- Queries: 6,093
- Index size: 2.158 GB
- Index build: 51.06 s
- Query time: 34.04 s
- Throughput: 178.9 queries/s

| Level (bp) | Queries | Top-1 species cluster | Top-5 species cluster | Top-1 genome | Forward orientation | Mean query coverage |
|---:|---:|---:|---:|---:|---:|---:|
| ALL | 6,093 | 99.10% | 99.70% | 92.52% | 96.68% | 0.9616 |
| 500,000 | 774 | 99.35% | 100.00% | 96.64% | 98.45% | 0.9995 |
| 250,000 | 799 | 99.37% | 100.00% | 96.75% | 98.50% | 0.9990 |
| 100,000 | 884 | 99.43% | 100.00% | 96.49% | 98.30% | 0.9950 |
| 50,000 | 1,026 | 99.51% | 100.00% | 95.03% | 98.05% | 0.9847 |
| 25,000 | 1,206 | 99.59% | 100.00% | 91.63% | 96.68% | 0.9618 |
| 10,000 | 1,404 | 97.86% | 98.65% | 84.12% | 92.59% | 0.8803 |

## Acceptance check

At the 5k scale:

- species-cluster top-1 accuracy for >=100 kb fragments is
  **99.43%**, exceeding the 95% minimum and 98% strong target;
- genome top-1 accuracy for >=250 kb fragments is **96.62%**,
  exceeding the 90% minimum but not the 98% strong target;
- placement output is deterministic and the index is reusable.

The <10 kb tier remains exploratory, as specified in the handoff.

## Preliminary baseline timing

On the 1k-scale query set, `minimap2 -x asm20` with 8 threads took 125.8 s
for 3,958 coordinate-window FASTA queries; Syn2b took 24.92 s with one thread.
The local OAPGC FASTAs contain many `-` alignment-gap characters, so this is
reported as timing context only and not as a valid accuracy comparison. A
clean FASTA-derived baseline needs regenerated raw sequences or another cohort.

## Artifacts

- 1k summary: `validation/oapgc_fragment_1k/summary.tsv`
- 1k details: `validation/oapgc_fragment_1k/detail.tsv`
- 5k summary: `validation/oapgc_fragment_5k/summary.tsv`
- 5k details: `validation/oapgc_fragment_5k/detail.tsv`
- Run statistics: `run_stats.json` in each run directory
- Leave-one-genome-out exploratory run: `validation/oapgc_holdout_1k/summary.tsv`

## Real-data issue found and fixed

The first real OAPGC index build exposed a bug in the MVP: a multi-enzyme
Type IIB TGT was incorrectly rejected as a mixed landmark source. This was
fixed. Enzyme-panel TGTs are now indexed; only enzyme-vs-FracMinHash mixing is
rejected.

The initial global-coordinate fragmentation test also showed why fragment
construction must respect TGT contig boundaries. After switching to
within-contig windows, source-present placement reached the values above.


## Update after genome-level chaining

The original chaining restricted every block to one reference contig. It now
chains by genome-level landmark coordinate while reporting the dominant
reference contig. The source-present tables above were rerun after this change.
Small genome-top gains were observed and the acceptance conclusions are
unchanged.

The formal phase status, 10k runtime, chimera pilot, and source-excluded
failure mode are consolidated in `validation/IMPLEMENTATION_STATUS.md`.
