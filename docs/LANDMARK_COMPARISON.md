# 2bRAD digestion vs FracMinHash: what changes, what does not

Everything downstream of landmark extraction in `synteny::scoring` consumes only
`(canonical identity, position, contig, orientation)`. Syn2b now supplies that from
either a Type IIB restriction digest or a FracMinHash sketch, so the two can be
measured against each other with the whole rest of the pipeline held fixed. This
records what that measurement found.

All numbers are from *E. coli* K-12 (4,543,028 bp) or from controlled synthetic
sequence where a variable had to be isolated. Every one is reproducible from the
commands in this document.

**Scope note.** The structural (SV) comparison below is measured. **The ANI
comparison is not yet run** — it is queued as task 8 in the paper repo's
`results/gtdb50k/HPC_TODO.md`, over the same 43k GTDB pairs the enzyme panel was
validated on. Nothing here should be read as an ANI result.

---

## 1. The two selection rules

| | 2bRAD (`--mode 2brad`) | FracMinHash (`--mode fracminhash`) |
|---|---|---|
| rule | sequence matches a Type IIB recognition site | `h(canonical(kmer)) < u64::MAX / scale` |
| context | fixed anchors in the sequence | the k-mer alone |
| density control | discrete: 1, 2, 4 or 16 enzymes | continuous, via `--scale` |
| wet-lab realisable | **yes** — this is a sequencing protocol | no, in silico only |
| needs | reads, from a mixed sample | an assembled genome |

The last two rows are the ones that decide which tool is for which job, and no
measurement below changes them.

Both are **context-free**: whether a locus is selected depends on that locus alone,
never on its neighbours. This is the property minimizers lack — minimizer selection
is window-relative, so one substitution re-selects a whole neighbourhood and landmark
identity stops being stable across genomes, which an adjacency-based metric cannot
tolerate. FracMinHash is additionally *genome-independent*, unlike bottom-*s* MinHash
whose cutoff is the *s*-th smallest hash of that particular genome.

---

## 2. Where they are identical

At matched density (FracMinHash `--scale 750` gives 6,034 landmarks against the
four-enzyme panel's 6,079), every structural control returns the same answer:

| control | 4-enzyme | FracMinHash | truth |
|---|---|---|---|
| self-comparison | 0 bp, 0 SCJ, obs 1.0000 | 0 bp, 0 SCJ, obs 1.0000 | collinear |
| origin rotation, 1.2 Mb | 0 bp, 0 SCJ, f 0.0000 | 0 bp, 0 SCJ, f 0.0000 | no change |
| one 500 kb inversion | 2 bp, 4 SCJ, f 0.10919 | 2 bp, 4 SCJ, f 0.11224 | 2, 4, f 0.1101 |
| 1 / 2 / 3 / 5 inversions | 2/4, 4/8, 6/12, 10/20 | identical | 2R, 4R |
| 40-contig shatter | 0 bp, obs 0.9930 | 0 bp, obs 0.9934 | 0 bp, K = 40 |

`breakpoints = 2R` and `scj_distance = 4R` hold exactly under both. The unit test
`fracminhash_landmarks_drive_the_same_structural_metric` asserts equality of every
output field on the same landmark set.

**This is the load-bearing result of the whole exercise.** The structural
mathematics — the fragmentation theorem, the direction-free adjacency sets, SCJ, the
orientation channel, the `>= 2`-landmark relocation rule, the error model's form —
is a property of the metric, not of restriction enzymes.

---

## 3. Where they differ

### 3.1 GC dependence — the largest difference, ~10x versus ~1x

Landmark count on 2 Mb of synthetic sequence at controlled GC:

| GC | BcgI | 4-enzyme | FMH s=750 | FMH s=1582 |
|---|---|---|---|---|
| 0.25 | 141 | 870 | 2,617 | 1,224 |
| 0.35 | 413 | 1,446 | 2,649 | 1,271 |
| 0.45 | 746 | 1,955 | 2,667 | 1,245 |
| 0.50 | 965 | 2,168 | 2,717 | 1,290 |
| 0.55 | 1,163 | 2,421 | 2,680 | 1,258 |
| 0.65 | 1,373 | 2,531 | 2,695 | 1,298 |
| 0.75 | 1,240 | 2,063 | 2,635 | 1,242 |
| **max/min** | **9.7x** | **2.9x** | **1.04x** | **1.06x** |

At GC 0.25 BcgI yields one landmark per 14.2 kb, which is below the useful floor for
this metric. The four-enzyme panel reduces the dependence but does not remove it, and
its response is non-monotonic, peaking at GC 0.65. FracMinHash varies by 4% across
the whole range, because its selection is a hash threshold and carries no sequence
preference.

Real bacteria span this range — *Buchnera* near 0.25, *Streptomyces* near 0.72 — so
this is not a synthetic edge case. It is the reason a fixed enzyme panel cannot be a
uniform instrument across a phylogenetically broad set.

```bash
syn2b digest -i gc25.fna -o /dev/null -e BcgI -f text
syn2b digest -i gc25.fna -o /dev/null --mode fracminhash --kmer 31 --scale 750 -f text
```

### 3.2 Density is a continuous knob on one side and four steps on the other

Observed against the prediction `4.54 Mb / scale`:

| scale | landmarks | expected | ratio |
|---|---|---|---|
| 250 | 18,196 | 18,172 | 1.001 |
| 500 | 9,120 | 9,086 | 1.004 |
| 750 | 6,034 | 6,057 | 0.996 |
| 1,000 | 4,539 | 4,543 | 0.999 |
| 2,000 | 2,234 | 2,272 | 0.983 |

The enzyme path offers 1, 2, 4 or 16 enzymes and nothing between. This matters
directly for the error model: `Var(err) = 1.504 p(1-p)/m + 0.0205^2` is a function of
landmark count `m`, and FracMinHash is the only source that can sweep `m` across a
decade with everything else held fixed. See task 8 in the paper repo — the
`--scale 1582` run, matched to BcgI's density, separates "fewer landmarks" from
"differently distributed landmarks" as explanations for that model's failure to
transfer to BcgI (z SD 1.08 vs 2.88).

### 3.3 Near-duplicate landmarks — the `sub_2` mechanism

The `>= 2`-landmark relocation rule exists because a landmark can appear to teleport
across megabases without anything being rearranged. The mechanism needs a unique
landmark sitting one substitution from a **multi-copy family**: the family is dropped
by the per-genome uniqueness filter, but in a diverged genome, once enough of its
copies are destroyed the survivor becomes unique and collides with the other locus.

| source | unique landmarks | multi-copy families | at risk | share |
|---|---|---|---|---|
| BcgI | 2,809 | 13 | 7 | 0.249% |
| four-enzyme panel | 5,889 | 38 | 20 | 0.340% |
| FracMinHash s=1582 | 2,776 | 17 | **0** | **0.000%** |
| FracMinHash s=750 | 5,880 | 39 | **0** | **0.000%** |

Read the third column before the fourth: **FracMinHash carries just as many genuine
multi-copy families.** Repeats are a property of the genome, not of the selection
rule. What differs is that none of its unique landmarks sits one substitution from
one. Enzyme landmarks must contain a recognition motif, so several thousand of them
are crammed into a small region of sequence space and near-collisions are
correspondingly likelier; FracMinHash k-mers are drawn from the whole 4^31 space with
no shared constraint.

The consequence is visible in the raw metric. On a substitution ladder with **no
rearrangement at all**, so every junction is false:

| substitution | 4-enzyme bp | 4-enzyme SCJ | FMH bp | FMH SCJ |
|---|---|---|---|---|
| 0.1% – 2% | 0 | 0 | 0 | 0 |
| 3% | 0 | **6** | 0 | **0** |
| 5% | 0 | **18** | 0 | **0** |

`breakpoints` is 0 for both because the relocation rule filters it. `scj_distance` is
the *unfiltered* symmetric difference, so it still shows the damage the rule is
hiding — on the enzyme path only. **The relocation rule protects the enzyme path from
a failure mode the sketch path does not have.**

### 3.4 Breakpoint localisation — spacing, and the gaps spacing hides

Three known 200 kb inversions on a closed genome, six true breakpoints:

| source | landmarks | spacing | median error | max error |
|---|---|---|---|---|
| BcgI | 2,872 | 1,582 bp | 614 bp | 3,671 bp |
| 4-enzyme | 6,079 | 747 bp | 273 bp | **3,671 bp** |
| FMH s=750 | 6,034 | 753 bp | 456 bp | **1,031 bp** |
| FMH s=200 | 22,708 | 200 bp | 248 bp | 783 bp |
| FMH s=50 | 90,394 | 50 bp | **44 bp** | **71 bp** |

The median tracks spacing on both, as it must — Syn2b reports the left landmark of
the broken adjacency, so its error is the gap to the next landmark.

The **maximum** is the interesting column. It does not improve from BcgI to the
four-enzyme panel: 3,671 bp both times. Adding enzymes cannot fill a gap that has no
sites, because they all cluster on their recognition motifs. FracMinHash at matched
density already reaches 1,031 bp, and at `--scale 50` localises to 71 bp — gene-level
resolution, which the enzyme path cannot reach at all, since 16 enzymes is the
ceiling.

The practical pattern is two-stage: screen at `--scale 750` or the four-enzyme panel,
then re-run only the interesting pairs at `--scale 50`.

### 3.5 Landmark retention under divergence

Share of landmarks surviving as shared, on the same ladder:

| substitution | 4-enzyme | FracMinHash s=750 |
|---|---|---|
| 0.1% | 89.5% | 94.6% |
| 1% | 67.7% | 71.1% |
| 3% | 36.4% | 37.6% |
| 5% | 19.1% | 19.0% |

FracMinHash retains more at low divergence and the two converge by 5%. A 31-mer is
destroyed only by a substitution inside it; an enzyme tag is destroyed by a
substitution in its recognition site as well, which is a larger effective target.

---

## 4. What has not been measured

- **ANI.** Task 8 in `Syn2bANI-paper/results/gtdb50k/HPC_TODO.md`, over the same 43k
  GTDB pairs. Until it runs there is no FracMinHash ANI result, and the enzyme-side
  ANI numbers must not be quoted as if they applied to both.
- **Multi-genome near-duplicate rates.** Section 3.3 is one genome. The at-risk
  fraction should be measured across taxa before it is stated as a general property.
- **Behaviour on reads rather than assemblies.** FracMinHash needs an assembly.
  Whether a read-level sketch reproduces these results is untested and is the
  question that decides whether it can ever substitute for 2bRAD in the lab.

---

## 5. What follows

**They are not competitors.** 2bRAD is a wet-lab protocol that delivers landmarks
directly from a mixed sample without assembly; FracMinHash is a computational
selection rule that needs an assembled genome. Nothing measured here lets one replace
the other.

What the comparison establishes is narrower and more useful: **the structural
mathematics is landmark-agnostic, and demonstrably so.** Every metric behaves
identically from two selection rules that share no mechanism. Where they differ, the
differences are properties of the *selection*, not of the metric, and each one is
explained by a stated mechanism rather than observed and left standing:

| difference | mechanism |
|---|---|
| GC dependence 9.7x vs 1.04x | motif frequency depends on base composition; a hash threshold does not |
| discrete vs continuous density | enzyme count vs a real-valued scale |
| 0.34% vs 0.00% at risk | motif constraint compresses the sequence space landmarks occupy |
| max localisation error 3,671 vs 1,031 bp | motif clustering leaves gaps that adding enzymes cannot fill |
| retention 89.5% vs 94.6% at 0.1% | the recognition site is part of the target |
