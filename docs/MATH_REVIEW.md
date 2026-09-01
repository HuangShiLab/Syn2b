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

So the adjacency metric gives a clean unbiased estimator of **event count**, and
*not* of the **conserved fraction**. Alignment-based synteny reports a fraction.
That gap — count → fraction — is the whole problem, and the document is right to
put its best section (§7) there.

The document's route to closing it is to infer extent from count under an assumed
event-length distribution. **There is a shorter one**: the orientation of each
tag relative to its canonical form is a second, independent channel that measures
inverted extent directly, with no length assumption. It is now implemented and
validated against exact truth at slope 1.0072, R² 0.9988 (§6). Junctions count
how *often* the genome moved; orientation measures how *much* of it moved.

Sections 1–4 and 5.1–5.3 are correct. Four specific defects follow, one of which
(§3) is now fixed — and §7 shows that three of the four were instances of one
general fact about statistics built on counts of transitions.

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

## 3. Defect 2 — §6.1 needed strand, which the pipeline destroyed *(now fixed)*

§6.1's probability that a segment of length `s` contains ≥1 event was derived
assuming orientation is observable. It was not, for two independent reasons:

1. The TGT **text format did not serialise strand** — the record was
   (sequence, position, contig), and the reader hard-coded every tag to
   `Strand::Forward`, so a text round-trip silently discarded the field.
2. `structural_synteny()` canonicalises every tag to `min(seq, revcomp)`, which
   erases orientation *by construction*. This is deliberate and correct for the
   adjacency metric — it is what makes the score invariant to
   reverse-complementing a whole genome.

Both are now addressed, and the fix turned out to be more informative than the
document anticipated:

- The digester previously hard-coded `Strand::Forward` for every tag, so the
  field was dead metadata occupying a byte of every binary record. It now
  records which recognition-site orientation actually matched (E. coli K-12,
  BcgI: 1461 forward sites, 1474 reverse). The text format serialises it as a
  `/+` or `/-` suffix, which older files simply lack and which defaults to
  forward, so every existing `.tgt` file still parses.
- **The inversion signal does not come from that field.** The digester always
  stores a tag as read off the forward strand of the assembly, so a locus inside
  an inverted segment is stored reverse-complemented. It still *matches*, because
  canonicalisation maps both forms to one identity, but the bit saying which of
  the two forms was stored has flipped. That bit is a function of the sequence
  alone, so it is recoverable from files written before the format change, and it
  is defined even for palindromic enzymes like AlfI, where site orientation is
  not.

The measured consequence is in §6 below.

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
of dropped repeat tags — not a binomial CI. That curve now exists and needs no
free parameters: see §8 Phase 2, and §7 for the fragmentation half of the same
statement. `repeats_dropped` and
`landmarks_collapsed` are already emitted for this purpose.

---

## 6. §7 asks the right question, but there is now a shorter answer

§7's breakpoint-derived conserved-fraction estimator `Ĉ_bp` is the section that
attacks the actual gap identified in §0 — turning a **count** into the
**fraction** that dnadiff/ANIm and SynTracker-APSS report. That is the right
target. Its route, however, is indirect: mapping `b` events onto a fraction of
genome requires assuming how long those events typically are, and the document
assumes a form rather than estimating it.

The orientation channel measures the fraction **directly**, with no
event-length assumption at all. Every landmark inside an inversion flips and
every landmark outside it does not, so the flipped share of shared landmarks
*is* the inverted share of the genome.

Measured against a ladder of R exactly-100 kb inversions on E. coli K-12
MG1655 (R = 1, 2, 3, 5, 8, 12, 20; BcgI; every interval recorded at construction
time), regressing the reported fraction on the true inverted base-pair fraction:

| R | true bp fraction | reported | ratio | junctions |
|---|---|---|---|---|
| 1  | 0.02154 | 0.02423 | 1.125 | 2  |
| 2  | 0.04309 | 0.04277 | 0.993 | 4  |
| 3  | 0.06463 | 0.06130 | 0.948 | 6  |
| 5  | 0.10772 | 0.11083 | 1.029 | 10 |
| 8  | 0.17235 | 0.17463 | 1.013 | 16 |
| 12 | 0.25853 | 0.24982 | 0.966 | 24 |
| 20 | 0.43088 | 0.43799 | 1.016 | 40 |

**slope 1.0072, intercept −0.00073, R² 0.9988.** The residual is landmark-
sampling noise and shrinks as 1/√(landmarks inside the event): 12.5% off at
R = 1 (68 landmarks), 1.6% at R = 20 (1229).

So the two channels are complementary and both exact in their own currency:
**junctions count how often the genome moved; the orientation fraction measures
how much of it moved.** Together they answer what neither answers alone, and
they do it without §7's length-distribution assumption.

Two limits, both stated rather than hidden:

- Flips are counted against the **majority** orientation, so reverse-
  complementing a whole assembly reads as 0.0 rather than 1.0. The price is a
  genuine identifiability limit: past 50% inversion the minority frame becomes
  the majority one and the fraction saturates. The junction count does not
  saturate, so the pair still separates those cases.
- A tag that is its own reverse complement reads the same in both orientations.
  These are counted and reported (`orientation_uninformative`) rather than
  silently dropped. On E. coli K-12 with BcgI there are none.

This does not make §7 worthless — `Ĉ_bp` remains the only route available to a
method that has counts but no orientation, and it is worth keeping as a baseline
to measure the orientation channel against. But it is no longer the primary
route, and the effort budgeted for it in Phase 3 should shrink accordingly.

## 7. A general result the specific critiques were instances of

Added 2026-08-30, after the Phase 2 grid and the GTDB50k external comparison.
Three of the four defects above, and two independent findings since, are the same
fact in different clothes. It is worth stating once, in general form, because it
decides which quantities the method should report.

### Setup

Every observation process here partitions a genome into contiguous observed
segments: contigs (assembly), 1-to-1 blocks (nucmer), chains (Syn2bANI anchors),
landmark runs (Syn2b). Let a genome carry `S` landmarks and let the process yield
`K` segments, hence `K − 1` internal boundaries at which adjacency information is
simply absent.

### Transition counts acquire a term linear in K

The observable adjacencies drop from `S` (circular) to `S − K`. So any statistic
that counts transitions — junctions, breakpoints, blocks, chain ends — has

    E[T] = T_true + c·(K − 1) + …

and `c` is decided by one design choice: **is an absent adjacency counted as a
contradicted one?**

| implementation | c | measured |
|---|---:|---|
| absence counted as a junction (Syn2b before the fix) | 1 | 119 false junctions on a 120-contig assembly of the genome itself |
| dnadiff `Breakpoints` (every 1-to-1 block boundary) | 1 | intercept **290** in `dnadiff_breakpoints = 5.35·b + 290`, and a median of 92 where Syn2bANI reports zero |
| segment count subtracted on one side only (Syn2bANI: query yes, reference no) | 1 for the other side | `+ (n_ref − 1)` exactly: 10 → 29 → 207 at n_ref = 1, 20, 200 |
| positive contradiction required (Syn2b now) | **0** | 0 junctions at every K up to 1000 |

Note the third row is not a bug so much as an asymmetry, and the fourth shows the
term is removable in principle, not merely reducible.

### Length-weighted ratios are invariant to K

Define `F = Σ_{i∈P} ℓ_i / Σ_i ℓ_i`, where `P` is a property that splitting
preserves — "this segment is inverted relative to the other genome" is one.
Splitting segment `i` into `i₁, i₂` with `ℓ_i = ℓ_{i₁} + ℓ_{i₂}` leaves both sums
exactly unchanged, so **`F` does not depend on `K` at all**. No correction term
exists to get wrong.

Measured: `inverted_fraction` reads **0.000** on a 120-contig assembly of the
genome itself, and regresses on true inverted base-pair fraction at slope 0.968,
R² 0.9993 across a 512× range of event lengths and 0–5% divergence.

The price is real and should be stated with it: a ratio carries no event count,
and it saturates — past 50% the minority frame becomes the majority one. So the
two are **complementary with disjoint failure modes**: counts are exact in event
number and fragile to fragmentation; ratios are robust to fragmentation and blind
to event number. Report both; do not derive one from the other.

When the goal is comparison with a fixed-reference alignment method (e.g.
dnadiff), the saturation can be avoided by reporting the orientation mismatch
fraction **relative to the chosen reference genome** instead of relative to the
majority frame. This `raw_inverted_fraction` ranges in [0, 1], loses whole-
genome reverse-complement invariance (a genome and its complement read as 1.0,
exactly as dnadiff reports), and becomes the direct analog of the alignment-
based inverted aligned fraction.

### Corollary: the power discount, in closed form

If the count is corrected rather than contaminated (`c = 0`), what remains is a
loss of *sensitivity*: junctions falling at segment boundaries are invisible.
That share is exactly the observable adjacency fraction,

    observable_fraction ≈ 1 − (K − 1) / S

exact to four decimals up to K = 300 on E. coli/BcgI, drifting below the formula
once segments hold fewer than about two landmarks. It is an unbiased predictor of
how many true junctions survive: `10 × observable_fraction` against observed, over
16 fragmentation levels, gives mean error −0.31, sd 1.26 — against a binomial
sampling sd of ≈1.5 at a truth of 10.

Practical form, scale-free in the panel: recovery is essentially complete while
segments hold **≳10 landmarks**, i.e. contig N50 ≳ 10× the landmark spacing.

### What this changes about §7 of the source document

The earlier argument against `Ĉ_bp` was that it needs an event-length prior it
cannot estimate (measured cost: 42× worse than the orientation channel across a
512× length range). The result above is a **second, independent** argument, and a
stronger one, because it bites even when the length prior happens to be right:
`Ĉ_bp` is built on a count, so it inherits the `c·(K − 1)` term from whatever
fragmented the input. The orientation channel is a ratio and does not.

### A prediction, stated before the data

The dnadiff `.1coords` files needed for this are already on the HPC. Computing the
inverted aligned fraction — `Σ|E2 − S2|` over reversed-query blocks over `Σ|E2 −
S2|` over all blocks — should, if the above is right, show:

1. **no material intercept** against Syn2b's `inverted_fraction`, where the count
   comparison showed 290; and
2. **correlation that does not improve when contig count is controlled**, since
   neither side carries a `K` term — unlike `breakpoint_count ~
   dnadiff_breakpoints`, which rose 0.465 → 0.534 under exactly that control.

In addition, the fixed-reference `raw_inverted_fraction` should correlate with
dnadiff across the full [0, 1] range, because it no longer flips the reference
frame at 50%.

A materially non-zero intercept would mean the invariance argument is wrong
somewhere, and that is worth finding out.

## 8. Research plan

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
- [x] **0.2** Serialise recognition-site orientation in the text format
  (`/+`, `/-`; absent in older files and defaulting to forward), make the
  digester record it instead of hard-coding forward, and add an
  orientation-mismatch signal derived from the sequence itself. See §6.

### Phase 1 — invariance controls *(done; one real bug found)*

Each transformation changes nothing biological, so each must return exactly 0
junctions and 0 inverted fraction. All four are now asserted as repo tests and
verified at genome scale on E. coli K-12 MG1655 (BcgI):

| Control | junctions | inverted_fraction | shared landmarks |
|---|---|---|---|
| self vs self | 0 | 0.00000 | 2806 |
| vs whole-genome reverse complement | 0 | 0.00000 | 2806 |
| vs 120-contig fragmented assembly | 0 | 0.00000 | 2804 |
| vs rotated circular origin | 0 | — | (unit test) |

Two defects were found by writing these rather than assuming them:

1. **Fragmentation inflated the junction count by one per contig break** — 119
   false junctions on a 120-contig draft. The cause is that an absent adjacency
   was being read as a broken one. A contig boundary *hides* an adjacency; it
   does not contradict it. A junction now requires positive contradiction: B must
   place other landmarks on both sides of one of the two partners. `scj_distance`
   is deliberately left uncorrected, since it is a published distance on
   adjacency sets, so it still reads 119 there and the two columns agree only for
   two closed genomes.
2. **The 40 bp overlap collapse was not reverse-complement symmetric** — it kept
   the first landmark of each run by position, and reversing the genome reverses
   which end a run starts from, so a genome and its own complement kept different
   survivors and stopped matching. That cost 64 of 2807 shared landmarks (2.3%).
   Runs now chain against the previous landmark rather than the run's first, and
   the representative is chosen by smallest canonical sequence, which is
   reversal-invariant. Shared landmarks are now identical (2806) across all three
   genome-scale controls.

### Phase 2 — detection power, measured not assumed *(done)*

Full results in [PHASE2_DETECTION_POWER.md](PHASE2_DETECTION_POWER.md); 72 cells
over event type × length (500 bp–256 kb) × substitution load (0–5%) × enzyme
panel, 40 events per cell, truth recorded at construction.

The result that replaces §5.4's confidence interval: **detection power is
predicted by landmark spacing alone**, with no free parameters — an inversion
needs ≥2 shared landmarks inside it, a translocation ≥1, and
`P(≥ m landmarks in a window)` fits the whole grid at **mean error +0.0051,
sd 0.0594**. So the resolution limit of any panel on any genome is *computable*:
digest once, count landmarks, read off the curve. L50 is predicted well; L95 is
under-predicted by ~1.5× because restriction sites clump rather than being
Poisson, so quote the measured L95.

Measured floors: inversion 2.6 kb (L50) / 12.7 kb (L95) with BcgI, 1.2 kb / 5.0 kb
with the four-enzyme panel; translocations at roughly half those sizes.
Specificity across the entire grid: **31 false-positive junctions out of 5,310
(0.58%)**. Events ≥32 kb are detected with power 1.000 at every substitution load
on both panels, except BcgI at 5% where tag survival falls to 20%.

### Phase 3 — `Ĉ_bp` as a baseline, not the primary route *(first result in)*

Quantified on the Phase 2 grid (see PHASE2_DETECTION_POWER.md §6). Any map from a
junction count to a genome fraction needs a mean event length λ, which counts
cannot supply. Against a fixed λ = 50 kb prior, the count-only estimate runs from
**12× over** (500 bp events) to **5× under** (256 kb events), for a mean absolute
error **42× worse** than the orientation channel — which needs no prior at all
and holds slope 0.968, R² 0.9993 across the same 512× range of event lengths and
0–5% divergence.

`Ĉ_bp` stays useful as the baseline a count-only method is stuck with, but it
should not be the primary estimator, and calibrating its length distribution is
not where the effort belongs.

**Still open:** the orientation channel measures *inverted* extent only. A
segment that moves without flipping flips no landmark, so translocation extent
remains count-only. Whether a comparable direct measure exists is the next
question worth attacking.

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
