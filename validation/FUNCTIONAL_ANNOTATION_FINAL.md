# OAPGC functional annotation and cluster search — final run

Date: 2026-10-10  
Scope: 2,474 representative SGB annotation sets and 1,488 representative genomes.

## Inputs

- Proteins: 6,466,011
- Amino acids: 10,585,678,783
- Prodigal GFF features: 6,466,011
- Complete OAPGC clusters represented: 1,488

## eggNOG-mapper

- emapper: 2.1.9
- eggNOG DB: 5.0.2
- DIAMOND: 2.1.10
- CPU: 14
- Wall time: ~10.1h for DIAMOND plus annotation
- Peak RSS observed: ~13.25 GB
- Database storage: 58 GB
- Result: 6,120,468 emapper queries

Annotation coverage:

| Field | Genes |
|---|---:|
| Description | 5,676,354 |
| KEGG KO | 3,795,669 |
| Pfam | 2,280,291 |
| eggNOG OGs | 6,120,468 |

Converted Syn2b feature rows: 5,866,571.  
Gene-neighborhood clusters indexed: 202,470.

## antiSMASH

- Version: 8.0.4 Docker image
- Genomes: 1,488
- Complete: 1,488
- Regions: 6,276
- CDS features inside regions: 136,196

Top region classes included terpene-precursor (1,403), RiPP-like (783), NRPS (778), arylpolyene (455), terpene (443), T1PKS (368), T3PKS (345), NRPS-like (285), ranthipeptide (198), and NRP-metallophore (197).

## MIBiG annotation

MIBiG 4.0 proteins were annotated with eggNOG:

```text
queries             45,185
with description    42,958
with KEGG KO        21,468
with Pfam           11,436
with eggNOG OGs     45,185
```

## Searchable indexes

| Index | Items | Size |
|---|---:|---:|
| OAPGC antiSMASH regions | 6,276 | 51 MB |
| OAPGC eggNOG gene neighborhoods | 202,470 | 5.2 GB |
| MIBiG eggNOG-enriched clusters | 2,636 | available |
| CARD cluster families | 517 | available |
| AMRFinderPlus + ResFinder targets | 125 | available |

## Benchmarks

### MIBiG self-rearrangement

Source MIBiG clusters were perturbed by loss, duplication, inversion, or truncation.

| Event | Queries | Hit fraction | Top-1 | Top-5 |
|---|---:|---:|---:|---:|
| ALL | 800 | 97.25% | 93.00% | 97.00% |
| duplication | 200 | 98.00% | 95.00% | 97.00% |
| inversion | 200 | 97.00% | 95.50% | 97.00% |
| loss | 200 | 97.00% | 93.50% | 97.00% |
| truncation | 200 | 97.00% | 88.00% | 97.00% |

This tests robustness to known-source perturbation, not held-out family novelty.

### Strict OAPGC-to-MIBiG retrieval

Using antiSMASH-enriched OAPGC clusters and strict thresholds:

```text
query clusters            5,851
clusters with MIBiG hits     20
```

This is not a sensitivity benchmark; it measures only high-confidence cross-references under strict thresholds.

### MIBiG held-out family benchmark

A true held-out family test was run for 145 families where the held-out cluster
had another cluster with the same broad MIBiG biosynthetic class/subclass.

```text
queries                    145
hit fraction             68.97%
top-1 family accuracy      0.00%
top-5 family accuracy      0.00%
```

This benchmark FAILED. The current token/adjacency score can recover perturbed
source clusters, but cannot yet assign held-out MIBiG families reliably when the
exact source cluster is absent. A homology/LSH family model or protein-level
similarity stage is required before claiming MIBiG family retrieval.

### AMR cross-database pilot

CARD, AMRFinderPlus 4.2, and ResFinder 2.6.0 were imported and indexed. Direct
schema crosswalk produced only two class-level hits because the databases use
different gene-family nomenclature and lack shared KO identifiers. A robust AMR
validation requires a curated cross-database family map or protein-level
alignment stage.

## Final artifacts

```text
external_data/oapgc_functional/oapgc_representatives.emapper.annotations
external_data/oapgc_functional/oapgc_eggnog_features.tsv
external_data/oapgc_functional/oapgc_eggnog_genes.s2cldb
external_data/oapgc_antismash/*/<genome>.json/.gbk/.zip
external_data/oapgc_functional/oapgc_antismash_features.tsv
external_data/oapgc_functional/oapgc_antismash.s2cldb
external_data/oapgc_functional/functional_annotation_summary.json
validation/mibig_rearrangement/summary.tsv
validation/mibig_family_holdout_versioned/summary.tsv
```

## Storage

Final additional functional-annotation storage:

```text
eggNOG database                  58 GB
OAPGC functional outputs         22 GB
antiSMASH outputs                17 GB
MIBiG functional/database data    ~5 GB
```

Internal disk remains with ~370 GB free.

# Final completion update — 2026-10-11

## 100k cluster runtime: COMPLETE

The unique-seed 100k stress database completed:

```text
clusters                         100,000
query clusters                    10,000
feature rows                     500,000
index size                       345.5 MB
cluster build                    3.07s
index build                      3.56s
search                          1,295.51s
query latency                   0.1296s/query
```

The original handoff wording “index query <1s” is met on a per-query basis
(0.13s/query at 100k). It is not a 1s total search for all 10,000 queries.

## MIBiG held-out family rescue: PARTIAL

A MMseqs2 protein-homology rescue stage was added after the token search.

At broad biosynthetic family level (all semicolon-delimited class terms):

```text
queries       145
protein hits  100%
top-1 family  18.62%
top-5 family  30.35%
```

At primary biosynthetic class level:

```text
protein hits  100%
top-1 family  62.76%
top-5 family  85.52%
```

The exact held-out family benchmark remains unsolved. The protein rescue stage
recovers the primary BGC class for most entries but not the full MIBiG family
label.

## AMR cross-database validation: COMPLETE

MMseqs2 crosswalks were run between CARD, AMRFinderPlus, and ResFinder.

| Search | Queries with hits | >=80% identity | >=90% identity |
|---|---:|---:|---:|
| CARD proteins vs AMRFinderPlus proteins | 6,024 | 97.19% | 97.05% |
| CARD CDS vs ResFinder CDS | 5,585 | 95.42% | 92.77% |

Best-hit outputs:

```text
external_data/oapgc_functional/amr_crosswalk_best_hits.tsv
external_data/oapgc_functional/amr_crosswalk_summary.json
external_data/amr_validation/card_vs_amrfinder.m8
external_data/amr_validation/card_vs_resfinder.m8
```

## Final phase status

| Task | Status |
|---|---|
| OAPGC eggNOG annotation | COMPLETE |
| OAPGC antiSMASH | COMPLETE |
| MIBiG eggNOG annotation | COMPLETE |
| Functional feature/index construction | COMPLETE |
| MIBiG self-rearrangement benchmark | COMPLETE |
| 100k cluster runtime | COMPLETE |
| AMR cross-database validation | COMPLETE |
| MIBiG held-out exact-family assignment | NOT SOLVED |
| GTDB / HROM full fragmentation matrix | NOT RUN |
| Phase 5 protein structure correlation | NOT RUN |

# 2026-10-12 orthologue and HPC update

MIBiG orthologue mapping is now COMPLETE. The 18.33-million-row OAPGC-vs-MIBiG
MMseqs2 search was reduced to 582,042 protein best hits and integrated into
`oapgc_antismash_mibig_features.tsv`. It covers 6,271 antiSMASH regions, 2,287
MIBiG accessions, and 140 MIBiG families. Strict orthologues with identity >=70%
and reciprocal query/target coverage >=50% number 22,708; 3,376 are >=90%
identity and >=80% both-coverages. See
`validation/MIBIG_ORTHOLOGUE_MAPPING.md`.

A gene-normalized protein/KO/Pfam ensemble improved the held-out MIBiG
primary-class top-1 from 60.00% to 95.17% (top-5 99.31%). Exact MIBiG family
assignment remains unsolved at 19.31% top-1 / 34.48% top-5 for the best
ensemble, so it must not be reported as exact-family retrieval.

The GTDB/HROM HPC benchmark pipeline has been installed and launched. GTDB has
143,586 references/1,000 selected queries; HROM uses all 142,277 TGT references
and 500 selected queries. GTDB digest array `4257501` and HROM fragmentation
job `4257718` were submitted on HPC; GTDB search `4257792` is chained with
`afterok:4257501`. Final full-matrix results are pending. Foldseek is installed
on HPC, but FASTA-only structural search is invalid; structure generation or
PDB/AlphaFold mapping is still required before Phase 5.
