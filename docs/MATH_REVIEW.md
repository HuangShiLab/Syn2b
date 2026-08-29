# Review: "Can TGT predict alignment-based synteny by analogy to Jaccard → ANI?"

Review of `TGT-synteny可预测性数学分析.docx` (2026-08-29), together with the
measurements that were run against it. Section numbers refer to the source
document.

---

## 0. Verdict

**The framework is sound, and the analogy is half-right in a way that matters.**

Jaccard → ANI works because k-mer survival is a *pointwise, independent* process:
each substitution independently destroys the k-mers covering it, so expected
containment is an analytic, invertible function of a single scalar `d`, and
`(1−d)^k` inverts cleanly.

Synteny is not pointwise. **A single inversion destroys exactly 2 adjacencies
regardless of its length.** The observable is therefore a function of the *number
of events*, not of their *extent*. This is not a nuisance — it is the structural
difference between the two problems:

| | Jaccard → ANI | TGT → synteny |
|---|---|---|
| Underlying process | pointwise, i.i.d. along the genome | sparse, extended events |
| Observable | fraction of surviving k-mers | count of broken adjacencies |
| Natural estimand | a rate (substitutions/site) | a **count** (events/genome) |
| What alignment methods report | ANI (a fraction) | aligned fraction / APSS (a **fraction**) |

So TGT gives a clean unbiased estimator of **event count**, and *not* of the
**conserved fraction**, unless the event-length distribution is separately
modelled. Alignment-based synteny reports a fraction. That gap — count → fraction
— is the whole problem, and §7's `Ĉ_bp` is exactly the right place to attack it.
Sections 1–4 and 5.1–5.3 are correct. Four specific defects follow.

---

## 1. What the measurements confirm

An inversion ladder was built on *E. coli* K-12 MG1655 (4,641,652 bp, verified —
**not** `data/e_coli_k12.fasta`, which is a truncated 4,500 KiB copy), R = 1…20
non-overlapping inversions, digested with the standard panel and scored with
`structural_synteny()` after the Phase-0 rewrite:

| R (true inversions) | junctions `b` | `scj_distance` | `breakpoint_density` | shared tags |
|---|---|---|---|---|
| 1  | 2  | 4  | 0.00071 | 2807 |
| 2  | 4  | 8  | 0.00143 | 2805 |
| 3  | 6  | 12 | 0.00214 | 2803 |
| 5  | 10 | 20 | 0.00357 | 2802 |
| 8  | 16 | 32 | 0.00572 | 2798 |
| 12 | 24 | 48 | 0.00861 | 2788 |
| 20 | 40 | 80 | 0.01444 | 2771 |

`b = 2R` and `scj_distance = 4R = 2b` hold **exactly**, with zero error, at every
point on the ladder. This is the doc's §5.2 claim (direction-free adjacency-set
symmetric difference *is* SCJ in the sense of Feijão & Meidanis, and SCJ = 2×
junction count in the balanced case) — confirmed empirically, not just cited
correctly. Junction coordinates come out in pairs, one pair per inversion
(e.g. R=3: 1084044/1154775, 2234213/2370384, 3374452/3551439).

The estimator is exact for the count. Everything downstream of it inherits that
exactness — which is precisely why the four defects below are worth fixing rather
than tolerating.

---

## 2. Defect 1 — §8's identifiability claim contradicts the document's own §5.3

§8 argues that the divergence of the *unchained* (unmatched) fraction can be
recovered from tag data. It cannot, and §5.3's own negative result for the spatial
model already says so.

A tag exists only where a restriction site survives. Once a site is destroyed the
locus becomes **invisible**, not "observed as diverged". The unmatched fraction
therefore contributes *zero* observations; it enters the likelihood only through
the count of missing tags. That count is confounded — a deleted region and a
hyper-diverged region both yield "no tag", and so does accessory content that was
never shared. Three distinct processes, one observable.

This is a genuine identifiability failure, not a small-sample problem: no amount
of data separates them without an external constraint (e.g. an assumed accessory
fraction, or a second enzyme panel with different site density used as a probe).
§8 should be rewritten as a **bound**, not an estimator: tag data constrains
`(deletion + divergence-beyond-site-loss)` jointly, and the split requires an
assumption that must be stated.

## 3. Defect 2 — §6.1's `1 − e^{−μs}` needs strand, which the pipeline destroys

§6.1's probability that a segment of length `s` contains ≥1 event is derived
assuming orientation is observable. It is not, for two independent reasons:

1. The TGT **text format does not serialise strand** — the record is
   (sequence, position, contig).
2. `structural_synteny()` canonicalises every tag to `min(seq, revcomp)`
   (`canonical_sequence()`), which erases orientation *by construction*. This is
   deliberate and correct for the adjacency metric — it is what makes the score
   invariant to reverse-complementing a whole genome — but it means the §6.1
   quantity is unobservable in the current pipeline.

Consequence: §6.1 is not wrong, it is **not yet measurable**. Phase 0 item 0.2
(add a `strand` field to the text format and derive an orientation-mismatch tag
count as a *separate* signal, leaving the canonical adjacency metric untouched)
is the prerequisite. Until then §6.1 should be marked as a proposed signal, not a
derived result.

## 4. Defect 3 — the chance correction `p₀ = 2/k` would erase the signal it corrects

The document imports a Mash-style chance-correction term, `p₀ = 2/k`. With
k ≈ 2807 landmarks this evaluates to **7.13 × 10⁻⁴**.

The measured `breakpoint_density` for a *single true inversion* is
**7.1 × 10⁻⁴**. Subtracting `p₀` would zero out one inversion exactly, and would
subtract ~5% of the signal at R = 20.

The correction is also unmotivated. In Mash it accounts for k-mers that collide by
chance in a small alphabet space. Here an adjacency's identity is an ordered pair
of canonical 32-mers; the chance-coincidence probability is on the order of 4⁻⁶⁴,
not 2/k. **Delete this term.**

The companion Jukes–Cantor-style multiple-hit correction is harmless but pointless:
across the entire ladder it inflates the estimate by 0.05% (R=1) to 0.98% (R=20).
At the densities where a tag-based method is usable, it is noise. Keep it only if
the method is ever pushed to densities where `b/k` approaches 0.1, and say so.

## 5. Defect 4 — §5.4's confidence interval models the wrong randomness

§5.4 places a binomial/Poisson interval around `b`, i.e. it models "which
adjacencies did we happen to sample". Under that model, at R = 20 with
`b ~ Bin(k, 2R/k)`, the standard deviation is **6.3 junctions**.

The measured standard deviation is **0** — `b = 2R` was exact at all seven ladder
points. A ±6-junction interval around a quantity measured without error is not
conservative, it is uninformative, and it hides the failure mode that actually
bites.

The real uncertainty is **structured, not stochastic**, and has two named sources:

- **Detection floor.** An event needs ≥2 landmarks inside it to be seen. With one
  landmark, canonicalisation maps the reverse-complemented tag to the same
  identity and the adjacency set is unchanged — the event is invisible. This sets
  a hard resolution limit of ≈ 8 kb (BcgI alone) / ≈ 4 kb (four-enzyme panel).
- **Repeat-mediated events.** Rearrangements whose boundaries fall inside repeated
  tags are dropped by the repeat filter, and are missed systematically rather than
  at random.

Both are *deterministic given the genome*, so the honest uncertainty statement is
a **detection-power curve as a function of event length**, plus a reported count
of dropped repeat tags — not a binomial CI. `repeats_dropped` and
`landmarks_collapsed` are already emitted for this purpose.

---

## 6. What is most valuable in the document: §7

§7's breakpoint-derived conserved-fraction estimator `Ĉ_bp` is the one section
that attacks the actual gap identified in §0 — converting a **count** into the
**fraction** that dnadiff/ANIm and SynTracker-APSS report. It is the only route by
which TGT numbers become directly comparable to alignment-based synteny, so it
deserves the effort.

Its validity, however, rests on an assumed event-length distribution: mapping
`b` events onto "fraction of genome in conserved blocks" requires knowing how long
those events typically are. The document assumes a form. That assumption is
testable and must be tested, not asserted — it is the single highest-value
experiment in the plan below.

---

## 7. Research plan

### Phase 0 — make the estimator well-defined (this commit)

- [x] **0.1** Direction-free adjacencies; respect contig boundaries; normalise
  circular origin. *(Ordered pairs gave 560 breakpoints for a 4-breakpoint event
  and scored 0.0000 on a whole-genome reverse-complement. Direction-free gives
  exactly 0 for all same-structure pairs and exactly 2 junctions per inversion at
  every substitution load.)*
- [x] **0.3** Collapse landmarks closer than 40 bp (overlapping cut sites from
  different enzymes at one locus; tag lengths differ by enzyme — BcgI/AlfI 32 bp,
  AloI/FalI 27 bp — so an inversion can swap them and fabricate a junction).
- [x] **0.4** Unify "breakpoint" semantics as **junction count** = |A \ B|, and
  report `scj_distance` = |A| + |B| − 2|A ∩ B| alongside it.
- [x] **0.5** Emit junction **coordinates**, not just counts
  (`<out>.junctions.tsv`).
- [ ] **0.2** Add `strand` to the TGT text format; derive an orientation-mismatch
  count as a **separate** signal. Prerequisite for §6.1.

### Phase 1 — invariance controls (must pass before any biology claim)

Each must return **exactly 0 junctions**, as repo tests:
self-vs-self; genome vs its reverse-complement; genome vs a fragmented assembly of
itself; circular genome vs a rotated origin. Both this repo and Syn2bANI have
failed at least one of these before; they are cheap and they are load-bearing.

### Phase 2 — detection power, measured not assumed

Simulate an event-length ladder (500 bp → 500 kb) × event-type (inversion,
translocation, deletion, insertion) × substitution load (0 → 5%) × enzyme panel.
Output the **power curve**: P(detected | length, type, panel). This replaces
§5.4's CI with the quantity that actually governs uncertainty, and it directly
tests whether the ≈4–8 kb floor derived above is right.

### Phase 3 — test §7's length-distribution assumption

On the same simulations, where the true conserved fraction is known by
construction, fit `Ĉ_bp` and measure its bias as the event-length distribution is
varied (fixed / exponential / heavy-tailed). If `Ĉ_bp` is robust only under one
family, that is a publishable limitation and must be stated, not buried.

### Phase 4 — external comparators on real closed genomes

No ground truth exists for real genomes; there are only comparators that measure
different things, and they must be treated as such:

| Comparator | What it actually measures | Use |
|---|---|---|
| `dnadiff` / ANIm aligned fraction | alignment-based conserved fraction | Phase-3 target for `Ĉ_bp` |
| SynTracker APSS | **divergence**, not order (measured Δ: −0.00016 structure / −0.08751 divergence) | *Not* a synteny target — a discriminant-validity control |
| Simulation | exact, by construction | primary truth |

The settled result stands: APSS responds to divergence, the corrected Syn2b
metric responds to structure (Δ −0.18053 structure / −0.00347 divergence) — they
are mirror images. **Agreement with APSS is not the goal, and the legacy metric's
r = +0.982 with APSS was evidence it was measuring the wrong thing.** The
defensible published claim is that two independent methods place breakpoints at
the same coordinates (≈1,544,640 / ≈1,945,050) at a ≈2790× cost difference.

### Known open issue

The `sub_2` residual: 6 false junctions on the four-enzyme panel, non-monotonic in
substitution load. Not explained by the 40 bp collapse. Must be resolved before
Phase 2 power curves are trusted.
