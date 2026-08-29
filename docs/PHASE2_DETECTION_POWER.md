# Phase 2 — detection power, measured

Answers the question §5.4 of the math analysis got wrong: what is the actual
uncertainty in a Syn2b call? Not a binomial interval around the junction count —
that quantity is measured with zero error (see [MATH_REVIEW.md](MATH_REVIEW.md)
§5). The uncertainty is **structured**: an event is either large enough to
contain landmarks or it is not.

Harness: [`scripts/detection_power.py`](../scripts/detection_power.py).
Raw grid: [`detection_power.tsv`](detection_power.tsv).
Reference: *E. coli* K-12 MG1655 (4,641,652 bp). 72 cells, 40 events per cell,
every event interval recorded at construction time.

---

## 1. The headline result: the detection floor is landmark spacing, and nothing else

An inversion needs **≥2** shared landmarks inside it to be seen — with one
landmark, canonicalisation maps the flipped tag back to the same identity and
both flanking adjacencies survive. A translocation needs only **≥1**, because the
segment moves away and breaks both its flanks regardless.

If that is the whole story, then power is just
`P(≥ m landmarks in a window of length x)` for Poisson landmarks at the observed
spacing — a formula with **no free parameters**, since spacing is
`genome length / shared landmarks`, both of which are known before any
simulation.

Tested against all 72 cells (including every substitution load, where spacing
rises because tags are lost):

> **mean error +0.0051, sd 0.0594** over the 35 cells where either the prediction
> or the observation is informative (not pinned at 0 or 1). All 37 saturated
> cells are correctly predicted saturated.

The practical consequence: **the resolution limit of any panel on any genome can
be computed, not simulated.** Digest once, count landmarks, read off the curve.

## 2. Power vs event length (no substitutions)

| event length | inversion BcgI | inversion panel4 | translocation BcgI | translocation panel4 |
|---:|---:|---:|---:|---:|
| 500 bp | 0.075 | 0.125 | 0.350 | 0.475 |
| 1 kb | 0.075 | 0.375 | 0.375 | 0.600 |
| 2 kb | 0.375 | 0.775 | 0.600 | 0.975 |
| 4 kb | 0.700 | 0.925 | 0.825 | 0.975 |
| 8 kb | 0.900 | **1.000** | 0.975 | **1.000** |
| 16 kb | 0.975 | 1.000 | **1.000** | 1.000 |
| 32 kb – 256 kb | **1.000** | 1.000 | 1.000 | 1.000 |

`panel4` = BcgI,AlfI,AloI,FalI (5724 shared landmarks, 811 bp spacing);
`BcgI` alone = 2805 landmarks, 1655 bp spacing.

| channel | spacing | L50 obs | L50 predicted | L95 obs | L95 predicted |
|---|---:|---:|---:|---:|---:|
| inversion, BcgI | 1655 | 2,611 | 2,777 | 12,699 | 7,850 |
| inversion, panel4 | 811 | 1,242 | 1,361 | 5,040 | 3,847 |
| translocation, BcgI | 1655 | 1,470 | 1,147 | 7,127 | 4,957 |
| translocation, panel4 | 811 | 574 | 562 | 1,910 | 2,429 |

L50 is predicted well. **L95 is under-predicted by roughly 1.5×**, and that gap is
real rather than noise: restriction sites are not Poisson-distributed along a
genome — GC content and repeats clump them — so the spacing distribution has a
heavier tail than the model, and the last 5% of events fall in the sparse
stretches. Quote the *measured* L95, and use the Poisson form only for L50 and
for ranking panels.

Translocations are detected at roughly half the size of inversions
(2,611/1,470 = 1.78; 1,242/574 = 2.16), which is the ≥2-vs-≥1 landmark
requirement showing up exactly where the mechanism predicts it.

## 3. Specificity

**31 false-positive junctions out of 5,310 across the entire grid (0.58%)** —
counting every cell, including 5% substitution load. A junction that matches no
true event boundary is rare enough that a single call is meaningful.

## 4. Substitution load

Divergence does not attack the metric directly; it removes tags, which widens
the spacing, which raises the floor. That is the only route, and it is visible:

**BcgI**

| substitutions | 2 kb | 8 kb | 32 kb | 128 kb | shared | spacing |
|---:|---:|---:|---:|---:|---:|---:|
| 0% | 0.375 | 0.900 | 1.000 | 1.000 | 2805 | 1655 |
| 0.5% | 0.400 | 0.875 | 1.000 | 1.000 | 2414 | 1923 |
| 1% | 0.325 | 0.925 | 1.000 | 1.000 | 2044 | 2271 |
| 2% | 0.150 | 0.775 | 1.000 | 1.000 | 1458 | 3184 |
| 5% | 0.025 | 0.200 | 0.875 | 1.000 | 572 | 8115 |

**panel4**

| substitutions | 2 kb | 8 kb | 32 kb | 128 kb | shared | spacing |
|---:|---:|---:|---:|---:|---:|---:|
| 0% | 0.775 | 1.000 | 1.000 | 1.000 | 5724 | 811 |
| 0.5% | 0.700 | 1.000 | 1.000 | 1.000 | 4893 | 949 |
| 1% | 0.650 | 0.975 | 1.000 | 1.000 | 4218 | 1100 |
| 2% | 0.400 | 0.950 | 1.000 | 1.000 | 3085 | 1505 |
| 5% | 0.100 | 0.600 | 1.000 | 1.000 | 1222 | 3798 |

**Events ≥32 kb are detected with power 1.000 at every substitution load tested,
on both panels** — except BcgI at 5%, where 20% tag survival is not enough. The
four-enzyme panel is worth its cost precisely here: at 5% divergence it still
holds 1.000 at 32 kb where BcgI has fallen to 0.875, and it recovers the 8 kb
band (0.600 vs 0.200).

## 5. The orientation channel across the whole grid

Reported inverted fraction regressed on true inverted base-pair fraction, over
event lengths spanning 500 bp to 256 kb and substitution loads 0–5%:

| cells | slope | intercept | R² |
|---|---:|---:|---:|
| all 52 inversion cells, to 5% divergence | 0.9679 | +0.00255 | 0.9993 |
| the 20 at 0% divergence | 0.9724 | +0.00094 | 0.9998 |

The ~3% shortfall is the small-event end: an event holding a handful of landmarks
quantises. It is systematic and correctable if it ever matters.

---

## 6. Phase 3 — what a count-only estimator pays for its length prior

Any map from a junction count to a genome fraction must supply a mean event
length: `fraction ≈ (b/2)·λ/L`. λ is not observable from counts, so it is
assumed. Using a fixed prior of λ = 50 kb against the same panel4 cells:

| true event length | orientation, obs/true | count + 50 kb prior, obs/true |
|---:|---:|---:|
| 500 bp | 0.93 | 12.5× |
| 2 kb | 1.06 | 19.4× |
| 8 kb | 1.01 | 6.3× |
| 32 kb | 1.02 | 1.6× |
| 128 kb | 0.97 | 0.39× |
| 256 kb | 0.98 | 0.20× |

The count-only route is exact only where the prior happens to match, and the grid
spans 512× in event length. Mean absolute error over the panel4 cells is **42×
worse** than the orientation channel, which has no such parameter.

This does not retire `Ĉ_bp` — it is the only route open to a method that has
counts and no orientation, and it stays useful as a baseline. But it should not
be the primary estimator, and the effort budgeted for calibrating its length
distribution is better spent elsewhere.

**Still open:** the orientation channel measures *inverted* extent. A segment
that moves without flipping does not flip any landmark, so translocation extent
is still count-only. Whether a comparable direct measure exists for it is the
next question.

---

## 7. How this was measured — four traps

All four were mistakes in the harness, not in the method, and each one produced a
plausible-looking wrong answer first. Recorded because anyone reproducing this
will hit them.

1. **Applying translocations one at a time corrupts the truth.** Each insertion
   shifts every coordinate to its right, so later excisions cut the wrong bases.
   The symptom was detection power *falling* with event size — 1.000 at 16 kb
   down to 0.175 at 256 kb — which is the opposite of what a resolution limit
   looks like. Build the variant as a single block permutation.
2. **A fixed detection window discards a predictable share of true calls.** A
   junction is reported at the left landmark of the broken adjacency, so it lands
   in `[boundary − gap, boundary]` with the gap exponentially distributed. A
   window of k × spacing therefore misses `e^−k` of real detections: at k = 3 the
   plateau reads 0.95 instead of 1.000, and looks like a method ceiling. Match
   each junction to its nearest unclaimed boundary instead.
3. **Count events placed, not events requested.** Packing constraints silently
   reduce how many events fit; scoring detections against the requested number
   understates power for exactly the middle band where packing is tightest.
4. **Packing events past 50% inverted hits the documented saturation.** The
   fraction is scored against the majority orientation, so a genome more than
   half inverted reports the minority. Cap the per-genome inverted fraction at
   ~0.30 when measuring the fraction estimator.
