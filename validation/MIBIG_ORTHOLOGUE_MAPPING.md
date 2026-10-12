# MIBiG orthologue mapping and family-model update — 2026-10-12

## Orthologue mapping: COMPLETE

The existing 18,334,092-row OAPGC antiSMASH-protein vs MIBiG 4.0 MMseqs2
search was reduced deterministically to one best hit per query protein. Hits
were tied by descending bitscore, identity, and alignment length. Query
coverage and target coverage were computed from the FASTA records. Results were
then joined back to the 136,196 antiSMASH CDS rows by genome/record/locus.

```text
MMseqs2 raw hits                         18,334,092
query proteins available                  923,113
query proteins with a best hit            582,042
OAPGC regions with any hit                  6,271
unique MIBiG accessions hit                 2,287
unique MIBiG families hit                     140
best-hit identity >=50%                    134,488  (23.11%)
best-hit identity >=70%                     23,434   (4.03%)
best-hit identity >=90%                      3,572   (0.61%)
strict: ID>=70%, qcov>=50%, tcov>=50%       22,708   (3.90%)
strict: ID>=90%, qcov>=80%, tcov>=80%        3,376   (0.58%)
```

Best-hit distributions:

```text
identity:  Q25=32.1%, median=39.1%, Q75=48.8%, Q95=67.5%, Q99=84.8%
query cov: median=93.25%, Q95=99.83%
target cov: median=92.98%, Q95=99.77%
bitscore: median=166, Q95=610, Q99=982
```

Artifacts:

```text
external_data/oapgc_functional/oapgc_vs_mibig_best_hits.tsv
external_data/oapgc_functional/oapgc_vs_mibig_summary.json
external_data/oapgc_functional/oapgc_antismash_mibig_features.tsv
scripts/summarize_mibig_orthologues.py
```

The integrated feature table is 38 MB and has the same rows/columns as the
antiSMASH feature table plus `mibig_orthologue`, accession, family, identity,
query/target coverage, and bitscore. The raw 3.1-GB `.m8` remains outside git.

## Protein-family ensemble: exact-family remains difficult

A corrected weighted protein model and an IDF cosine model over KO/Pfam tokens
were evaluated on the existing 145-cluster MIBiG held-out benchmark. The
normalized protein score sums coverage-weighted bitscores after one best HSP
per query/target protein pair and divides by the geometric mean of matched
query and target gene counts.

| Model | Exact top-1 | Exact top-5 | Primary top-1 | Primary top-5 | Broad top-1 | Broad top-5 |
|---|---:|---:|---:|---:|---:|---:|
| Weighted protein best | 17.93% | 31.72% | 84.14% | 97.93% | 71.72% | 96.55% |
| Gene-normalized protein | 19.31% | 33.79% | **95.17%** | **99.31%** | 83.45% | 97.93% |
| KO/Pfam IDF cosine | 17.24% | 33.79% | 93.79% | 97.93% | 79.31% | 95.86% |
| Protein 0.8 + profile 0.2 | 19.31% | 34.48% | 95.17% | 99.31% | **84.14%** | 97.93% |
| Protein 0.6 + profile 0.4 | 17.24% | **35.86%** | 93.79% | 99.31% | 82.07% | 97.93% |

This substantially improves primary-class assignment over the prior protein
sum-bits model (60.00% -> 95.17% primary top-1). Exact MIBiG family assignment
still fails the original >=80% target (19.31% top-1), as expected for highly
granular MIBiG version/family labels and distant homologues.

Artifacts:

```text
validation/mibig_family_holdout_versioned/family_ensemble_summary.tsv
validation/mibig_family_holdout_versioned/family_ensemble_detail.*.tsv
scripts/evaluate_mibig_family_ensemble.py
```

## HPC GTDB/HROM benchmark

The benchmark code was synced to
`/lustre1/g/aos_shihuang/syn2b_benchmarks/Syn2b`, built there, and its full
152-test suite passed. A `flate2` FASTA reader was added so GTDB `.fna.gz`
inputs can be used directly. Deterministic metadata and query lists were
prepared:

```text
GTDB references      143,586; query genomes 1,000
HROM references      142,277; query genomes   500
random seed          20261012
```

HROM uses the existing 142,277 TGT records directly. GTDB digesting is running
as SLURM array `4257501` (29 x 5,000 genomes). The HROM source-present
fragmentation job was resubmitted as `4257552` after fixing multi-contig genome
IDs. The GTDB fragmentation job will be launched after its digest array and
index build finish.

Planned HPC resource envelope:

```text
GTDB TGT text                 40-70 GB
GTDB index                     25-40 GB
GTDB 1,000-genome fragments     5-9 GB FASTA + 2-4 GB TGT
GTDB peak search RAM           35-70 GB
HROM index                     25-40 GB
HROM 500-genome fragments      3-6 GB FASTA + 1-3 GB TGT
HROM peak search RAM           30-60 GB
requested jobs                 digest 29x4 cores/16 GB; search 16 cores/140-180 GB
```

These totals fit comfortably within the group's remaining ~5.35-TB quota and
the 166-TB Lustre filesystem free space.
