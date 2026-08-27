# Syn2b → Syn2bANI integration roadmap

**Goal.** Use Syn2b as the research sandbox for TGT-core synteny metrics, and feed only the metrics that are demonstrably useful into Syn2bANI. Syn2bANI must stay focused on fast, stable, production-grade genome search; Syn2b can explore richer (and possibly slower) models without destabilising the paper pipeline.

**Current state (2026-08-27).** Syn2b already implements the key insight:
- `structural_synteny()` separates **substitution load** from **structural change** by using canonical tag identity, shared-tag restriction, ordered adjacency, and repeat dropping.
- On simulated *E. coli* K-12 genomes it scores exactly 1.0 under substitutions-only and ~0.82 under substitutions + 400 kb inversion, with breakpoint_density separated ~50×.
- This is a stronger theoretical foundation for a synteny score than the current Syn2bANI `anchor_adjacency`, which is confounded by coverage.

## 1. Why a separate Syn2b project is useful

| Concern | Syn2bANI (production) | Syn2b (research) |
|---|---|---|
| Output stability | Columns must not change after publication | Metrics can be added/renamed freely |
| Speed | <100 ms per pair, <1 h per GTDB-scale batch | Can afford seconds per pair for validation |
| Scope | ANI + a few validated structural metrics | Explore all TGT-derived quantities |
| Calibration | Ridge/GBRT trained on GTDB truth | Theory + simulation + alignment truth |
| Code maturity | Optimised, well-tested, CLI-frozen | Prototype-friendly |

## 2. Metrics to explore in Syn2b

Each metric should be evaluated against **dnadiff/ANIm alignment-based synteny** on GTDB-R207 held-out pairs and against **exact-truth simulations** (inversion/indel ladders).

### 2.1 Coverage-aware synteny scores

The main weakness of Syn2bANI `anchor_adjacency` is that it ignores coverage. Syn2b can test:

1. **`structural_synteny().score`** (already implemented). The coverage is implicit in `shared_tags`, but the score itself is conditioned on shared tags.
2. **`structural_synteny().score × coverage_proxy`** where `coverage_proxy` is the fraction of query tags that are shared and unique.
3. **Expected-observed adjacency ratio.** Compute the expected number of conserved adjacencies under a random permutation of shared tags, then report observed / expected.
4. **Chain N50 / L50 on shared-tag order.** Longer collinear blocks → higher synteny.

### 2.2 Breakpoint-style metrics

These already correlate best with dnadiff breakpoints in Syn2bANI:

- `breakpoint_density` from `structural_synteny()` (normalised by shared tags).
- `breakpoint_count` from ordered adjacency difference.
- Weighted breakpoint count using inter-tag distance (large gaps count more).

### 2.3 Rank-correlation metrics

- **Kendall's τ on shared-tag order.** Already implemented as `kendall_tag_order()`. Should be compared with `structural_synteny().score` on exact-truth simulations.

### 2.4 Alignment-based synteny estimators

The long-term goal: estimate the **dnadiff synteny score** (base-pair fraction in collinear alignable blocks) from sparse TGT data.

Candidate approaches:

- **Calibrated linear model.** Train ridge regression on TGT features (`shared_tag_fraction`, `structural_synteny.score`, `breakpoint_density`, `kendall_tau`, `chain_N50`, `mean_gap`, `repeat_fraction`) to predict dnadiff synteny, analogous to Syn2bANI's ANI calibration.
- **Generative model.** Model the genome as a Poisson process of breakpoints; likelihood of observed tag order given breakpoint rate; infer expected syntenic coverage.
- **Hybrid coverage + order.** `synteny_est = f(coverage, order_conservation)` where coverage comes from shared-tag fraction and order conservation from `structural_synteny.score`.

## 3. Validation experiments

### 3.1 Exact-truth simulation (required for any metric)

Use the existing `scripts/simulate_rearrangement.py` framework:
- Fixed substitution rate (e.g. 0%, 0.1%, 1%, 5%) × 0–32 inversions (100 kb–1 Mb) × 0–10 large indels.
- For each condition, compute all candidate metrics and report:
  - Invariance to substitutions (structural metrics should be flat).
  - Monotonic response to inversion/indel count.
  - Dynamic range (can it separate 0 from 1 inversion?).

### 3.2 GTDB-R207 held-out benchmark

- Inputs: the same 43,334 pairs used in Syn2bANI (`results/gtdb50k/`).
- Truth: dnadiff/ANIm aligned fraction and dnadiff structural breakpoints.
- For each candidate metric report Pearson/Spearman correlation with dnadiff synteny and MAE if calibrated.

### 3.3 High-ANI discordant pairs

- Focus on pairs with ANI >99% but low alignment-based synteny.
- Question: does the new metric flag these as low-synteny more reliably than `anchor_adjacency`?

## 4. Criteria for promoting a metric to Syn2bANI

A metric graduates from Syn2b to Syn2bANI only if:

1. **It is computable in the same single pass as ANI** (no second alignment step).
2. **It is stable across the GTDB 50k benchmark** (no pathological outliers).
3. **It improves on `anchor_adjacency` for at least one concrete use case** (e.g. better separation of high-ANI rearranged pairs, or better correlation with dnadiff synteny).
4. **Its interpretation is simple enough for the manuscript** (one sentence, one figure).
5. **It has a clear failure-mode description** (when should users not trust it?).

## 5. Immediate next steps

1. **Export `structural_synteny()` to a pairwise CLI command** so it can be run on the GTDB 50k pairs and compared with Syn2bANI output. Currently `synteny` operates on a directory of TGTs; a `pairwise` subcommand or script wrapper would make benchmarking easier.
2. **Write a Python bridge** that digests a FASTA pair with Syn2b and outputs the full `StructuralSynteny` struct plus `kendall_tag_order` and `adjacency_jaccard`.
3. **Run the GTDB 50k benchmark** with this bridge and produce a report comparing every Syn2b metric against dnadiff truth.
4. **Test the calibrated synteny estimator** (ridge on TGT features → dnadiff synteny) and report MAE.
5. **If a metric passes the promotion criteria**, open an issue in Syn2bANI proposing a new output column (e.g. `synteny_index`) with a precise definition and validation numbers.

## 6. Relationship to the Syn2bANI paper

The paper should **not** wait for Syn2b results. The current Syn2bANI story is:
- Fast ANI estimation.
- One-pass structural signals: `af_query` (coverage), `anchor_adjacency` (anchor-order), `breakpoint_count` (rearrangement count).
- Honest caveat: `anchor_adjacency` is not alignment-based synteny; the diagnostic plot shows why.

Future Syn2b results can appear in:
- A revised manuscript if they significantly improve the story.
- A dedicated Syn2b paper if they become extensive.
- Syn2bANI v2 if a metric graduates.

## 7. Files and scripts

- Syn2b core: `src/synteny/scoring.rs`
- Existing validation: `validation/structural_invariance.sh`
- Simulation: `scripts/simulate_rearrangement.py`
- Proposed new benchmark script: `scripts/benchmark_vs_dnadiff.py` (to be written)
