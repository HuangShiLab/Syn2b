//! Synteny scoring functions
//!
//! Provides metrics for evaluating synteny conservation between genomes,
//! including Jaccard similarity on tag adjacencies, Kendall's tau rank
//! correlation on tag order, breakpoint counting, and path-level scoring.

use crate::synteny::graph::TagAdjacencyGraph;
use crate::tgt::record::TgtRecord;
use crate::tgt::tag::Tag;
use std::collections::{HashMap, HashSet};

// ─────────────────────────────────────────────────────────────────────────────
// Synteny score for a linear path
// ─────────────────────────────────────────────────────────────────────────────

/// Score synteny conservation for a linear path.
///
/// The score is the average edge weight along the path divided by the
/// maximum possible weight (the number of genomes in the graph), adjusted
/// by a length bonus that rewards longer conserved segments:
///
///   score = (mean_weight / num_genomes) * min(1.0, sqrt(path_len / 10.0))
///
/// Higher score = more conserved synteny.
/// - A score of 1.0 means all genomes agree on every adjacency in the path
///   and the path is sufficiently long.
/// - A score of 0.0 means no synteny conservation.
pub fn synteny_score(path: &[u64], graph: &TagAdjacencyGraph) -> f64 {
    if path.len() < 2 {
        return 0.0;
    }

    let mut total_weight = 0u32;
    let mut edge_count = 0usize;

    // Sum edge weights along the path
    for window in path.windows(2) {
        let src = window[0];
        let tgt = window[1];
        // Try both edge directions (undirected path)
        let weight = graph
            .edges
            .get(&(src, tgt))
            .or_else(|| graph.edges.get(&(tgt, src)))
            .map_or(0, |e| e.weight);
        total_weight += weight;
        edge_count += 1;
    }

    if edge_count == 0 || graph.num_genomes == 0 {
        return 0.0;
    }

    let mean_weight = total_weight as f64 / edge_count as f64;
    let max_weight = graph.num_genomes as f64;
    let weight_ratio = mean_weight / max_weight;

    // Length bonus: longer paths score slightly higher, saturating at len=10
    let path_len = path.len() as f64;
    let length_factor = (path_len / 10.0).sqrt().min(1.0);

    weight_ratio * length_factor
}

// ─────────────────────────────────────────────────────────────────────────────
// Pairwise synteny matrix
// ─────────────────────────────────────────────────────────────────────────────

/// Compute pairwise synteny matrix for all genomes in the graph.
///
/// For each pair of genomes (g_i, g_j), computes the Jaccard similarity of
/// their adjacency sets: the set of (tag_a, tag_b) pairs that are adjacent
/// in each genome.
///
/// Returns a map from (genome_i, genome_j) to synteny score in [0, 1].
pub fn pairwise_synteny_matrix(
    graph: &TagAdjacencyGraph,
) -> HashMap<(String, String), f64> {
    let mut matrix = HashMap::new();

    if graph.num_genomes < 2 {
        return matrix;
    }

    // Build adjacency sets for each genome
    let genome_ids: Vec<String> = graph
        .nodes
        .values()
        .flat_map(|n| n.positions.keys().cloned())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    for i in 0..genome_ids.len() {
        for j in (i + 1)..genome_ids.len() {
            let gi = &genome_ids[i];
            let gj = &genome_ids[j];

            let score = pairwise_score(graph, gi, gj);
            matrix.insert((gi.clone(), gj.clone()), score);
            matrix.insert((gj.clone(), gi.clone()), score);
        }
    }

    matrix
}

/// Compute the pairwise synteny score between two genomes using their
/// shared adjacencies in the graph.
fn pairwise_score(graph: &TagAdjacencyGraph, genome_a: &str, genome_b: &str) -> f64 {
    // Collect tag adjacencies for genome A
    let mut adj_a: HashSet<(u64, u64)> = HashSet::new();
    let mut adj_b: HashSet<(u64, u64)> = HashSet::new();

    // Build adjacency sets by traversing the genome order
    if let Some(order_a) = genome_order_in_graph(graph, genome_a) {
        for w in order_a.windows(2) {
            let (u, v) = (w[0], w[1]);
            // Normalize: smaller id first
            adj_a.insert(if u < v { (u, v) } else { (v, u) });
        }
    }

    if let Some(order_b) = genome_order_in_graph(graph, genome_b) {
        for w in order_b.windows(2) {
            let (u, v) = (w[0], w[1]);
            adj_b.insert(if u < v { (u, v) } else { (v, u) });
        }
    }

    if adj_a.is_empty() && adj_b.is_empty() {
        return 0.0;
    }

    let intersection: HashSet<_> = adj_a.intersection(&adj_b).collect();
    let union: HashSet<_> = adj_a.union(&adj_b).collect();

    intersection.len() as f64 / union.len() as f64
}

/// Reconstruct the tag ID ordering for a genome from the graph's node data.
/// Structural synteny between two genomes, isolated from substitution load.
///
/// # The problem this replaces
///
/// [`pairwise_score`] builds adjacency sets over *all* consecutive tags of each
/// genome and compares them as **unordered** pairs. Both choices are wrong for a
/// structural metric, and measured on E. coli K-12 (BcgI, 2935 tags,
/// substitutions only, no structural variation at all) the consequence is total:
///
/// ```text
///   popANI    pairwise_score    predicted by tag loss alone
///  100.00%          1.0000                    1.0000
///   99.90%          0.8678                    0.8832
///   99.00%          0.3438                    0.3565
///   95.00%          0.0110                    0.0191
/// ```
///
/// The score *is* the tag-survival curve: a 32 bp tag survives 1% divergence
/// with probability 0.99^32 = 0.725, an adjacency needs both flanking tags, so
/// 0.725^2 = 0.526 survive and the Jaccard of the two sets is
/// 0.526/(2 - 0.526) = 0.357. It measures substitution load, not structure.
///
/// # What this function does differently
///
/// Three changes, each necessary and none sufficient alone — verified by
/// measuring all four combinations:
///
/// 1. **Canonical tag identity.** Without it, tags inside an inverted segment
///    are read from the other strand and drop out of the shared set entirely, so
///    the inversion signal travels the same path as substitution loss and the two
///    cannot be separated.
/// 2. **Restriction to shared tags.** Tag *presence* is a function of sequence
///    divergence; tag *order* among survivors is a function of structure. Only
///    the second belongs in a structural metric. ([`breakpoint_count`] already
///    does this; `pairwise_score` does not, which is why the two disagree on the
///    same input.)
/// 3. **Direction-free adjacency, one contig at a time.** An adjacency is the
///    unordered pair {a, b}. A chromosome read backwards is the same chromosome,
///    and in a draft assembly every contig's orientation is arbitrary. The
///    ordered variant this function used to build scored a genome against its own
///    reverse complement as 0.0000, and reported 560 junctions for a 400 kb
///    inversion whose true junction count is 2 — it was counting every adjacency
///    *inside* the inverted segment. Adjacencies are never formed across a contig
///    boundary, where no adjacency exists.
/// 4. **Overlapping cut sites collapsed.** Different enzymes cut the same locus:
///    with the four-enzyme panel, 305 adjacent landmark pairs on E. coli K-12 sit
///    under 40 bp apart while the tags themselves are 27–32 bp long. Their
///    relative order is an artifact, and since tag length varies by enzyme an
///    inversion shifts them by different amounts and can swap them. Collapsing
///    them took a 400 kb inversion from 8 reported junctions to the correct 4
///    (total absolute error over 12 test genomes: 14 → 6).
/// 5. **Circular origin normalised.** Bacterial chromosomes are circular and
///    assemblies begin at an arbitrary base, so a single-contig genome is closed
///    into a cycle. Without this, two identical genomes differing only in where
///    the assembly starts report two junctions.
///
/// Measured on genome suites built from E. coli K-12 MG1655 (ENA U00096.3) with
/// every event recorded at construction time:
///
/// ```text
///   construction                       junctions   scj_distance
///   substitutions only, 0.5% to 5%             0              0
///   one 400 kb inversion                       2              4
///   five inversions                           10             20
///   twenty inversions                         40             80
///   one 100 kb translocation                   3              6
/// ```
///
/// The junction count is exact and follows the rearrangement-theory convention
/// that one inversion cuts the chromosome in two places, so it agrees with
/// Syn2bANI's `breakpoint_count` and with nucmer/SyRI. It stays at exactly zero
/// under substitution loads up to 5%, and stays exact up to at least twenty
/// events. `scj_distance` is the symmetric difference of the two adjacency sets —
/// the published single-cut-or-join distance — which counts each junction twice
/// because every cut destroys one adjacency and creates another.
///
/// # Extent, as opposed to count
///
/// The junction count deliberately cannot see how *large* an event is: one
/// inversion breaks exactly two adjacencies whether it spans 5 kb or 500 kb.
/// That is a problem when comparing against alignment-based synteny, which
/// reports a **fraction** of the genome rather than a count of events.
///
/// [`StructuralSynteny::inverted_fraction`] supplies the missing axis. The
/// digester always stores a tag as read off the forward strand, so a locus
/// inside an inverted segment is stored reverse-complemented; it still matches,
/// because canonicalisation maps both to one identity, but the bit saying which
/// of the two forms was stored has flipped. Every landmark inside the inversion
/// flips and every landmark outside it does not, so counting flips measures how
/// much of the genome moved while the junction count measures how often.
///
/// Measured against a ladder of R exactly-100 kb inversions on E. coli K-12
/// (R = 1, 2, 3, 5, 8, 12, 20), regressing the reported fraction on the true
/// inverted base-pair fraction: **slope 1.0072, intercept −0.00073, R² 0.9988**.
/// The residual is landmark-sampling noise and shrinks as 1/sqrt(landmarks
/// inside): 12.5% off at R = 1 (68 landmarks), 1.6% at R = 20 (1229).
///
/// Flips are counted against the *majority* orientation, so reverse-
/// complementing a whole assembly — a strand convention, not biology — reads as
/// 0.0 rather than 1.0. The price is a real identifiability limit: past 50%
/// inversion the minority frame becomes the majority one and the fraction
/// saturates. The junction count does not saturate, so the pair still separates
/// those cases.
///
/// # Resolution limit
///
/// An event is invisible to the *junction* count unless at least **two**
/// landmarks fall inside it: with one landmark, the inversion
/// reverse-complements that tag, canonicalisation maps it back to the same
/// identity, and both flanking adjacencies are unchanged. Measured
/// 95%-detection sizes: ~8 kb for BcgI alone, ~4 kb for the four-enzyme panel.
/// This is a sampling limit, not an implementation defect. Note that the
/// orientation signal has a *lower* floor — a single landmark inside an
/// inversion still flips — at the cost of not localising the event.
///
/// A tag that is its own reverse complement reads the same in both orientations
/// and is counted in [`StructuralSynteny::orientation_uninformative`] rather
/// than silently dropped. On E. coli K-12 with BcgI there are none.
///
/// Returns `None` when fewer than two tags are shared, since no adjacency exists
/// to compare.
/// Minimum separation between consecutive landmarks, in base pairs.
///
/// Type IIB enzymes cut at overlapping loci. On E. coli K-12 with the
/// BcgI,AlfI,AloI,FalI panel, 305 adjacent landmark pairs sit under 40 bp apart
/// while the tags are 27–32 bp long, so they physically overlap and describe one
/// locus. Single-enzyme panels have uniform tag length and are unaffected.
const MIN_TAG_SEPARATION: u64 = 40;

/// One landmark: a tag reduced to what the order metric needs, plus the
/// orientation bit that the canonical form would otherwise hide.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Landmark {
    seq: [u8; 32],
    pos: u64,
    contig: u16,
    /// The stored forward-strand window is the reverse complement of `seq`.
    rc: bool,
    /// The tag is its own reverse complement, so `rc` is uninformative here.
    palindromic: bool,
}

/// An adjacency, keyed so that {a, b} and {b, a} are the same entry.
type Adjacency = ([u8; 32], [u8; 32]);

fn undirected(a: [u8; 32], b: [u8; 32]) -> Adjacency {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Landmarks in genome order, one per tag. Sorting is defensive: callers should
/// not have to guarantee it.
fn raw_landmarks(record: &TgtRecord) -> Vec<Landmark> {
    let mut all: Vec<Landmark> = record
        .tags
        .iter()
        .map(|t| Landmark {
            seq: t.canonical_sequence(),
            pos: t.position,
            contig: t.contig_id,
            rc: t.is_revcomp_of_canonical(),
            palindromic: t.is_palindromic(),
        })
        .collect();
    all.sort_by_key(|lm| (lm.contig, lm.pos));
    all
}

/// Collapse overlapping cut sites to one representative per run.
///
/// **Called after the shared-tag restriction, not before.** The run boundaries
/// are found by chaining against the previous landmark, and the representative is
/// the smallest canonical sequence in the run, so both depend on which members
/// are present. Collapsing first meant each genome collapsed its *own* tag set:
/// a tag the other genome had lost could still decide the representative, and the
/// two genomes then disagreed about a locus they both carried. Restricting first
/// makes the run structure agree by construction.
///
/// Measured on E. coli K-12 against a 5%-substituted copy, four-enzyme panel:
/// shared landmarks 1,230 -> 1,276 (+3.7%), and +1% to +2% across the rest of the
/// ladder. `landmarks_collapsed` now counts collapses among shared landmarks only,
/// so it reads much lower than before (610 -> 50 at 5%); it is a diagnostic, not
/// a rate.
///
/// This is **not** the `sub_2` residual. That hypothesis was tested by A/B against
/// a46badb over a six-point substitution ladder and rejected: the residual is 6
/// false junctions before and after, unchanged. See MATH_REVIEW.md for what it
/// actually is.
fn collapse_runs(all: Vec<Landmark>) -> Vec<Landmark> {
    // Group into runs of landmarks that chain together within MIN_TAG_SEPARATION,
    // then keep one representative per run. Two choices here have to be made
    // reverse-complement symmetric, or the same locus yields different survivors
    // in a genome and its complement and the two stop matching:
    //   - chain against the *previous* landmark, not the run's first, because
    //     reversing the coordinate order reverses which end the run starts from;
    //   - pick the representative by smallest canonical sequence, not by
    //     position, because the run's canonical set is reversal-invariant while
    //     its order is not.
    // Keeping the first by position instead cost 64 of 2807 shared landmarks
    // (2.3%) on an E. coli K-12 genome-versus-its-own-reverse-complement control.
    let mut kept: Vec<Landmark> = Vec::with_capacity(all.len());
    let mut run_start = 0usize;
    for i in 0..=all.len() {
        let breaks = i == all.len()
            || i == run_start
            || all[i].contig != all[i - 1].contig
            || all[i].pos.saturating_sub(all[i - 1].pos) >= MIN_TAG_SEPARATION;
        if breaks && i > run_start {
            let rep = all[run_start..i]
                .iter()
                .min_by(|x, y| x.seq.cmp(&y.seq))
                .copied()
                .expect("non-empty run");
            kept.push(rep);
            run_start = i;
        }
    }
    kept
}

/// Adjacency set of a landmark series, each mapped to the position of its left
/// landmark so a broken adjacency can be reported as a coordinate.
///
/// Closure is decided **per contig**, not per genome. A long-read assembly is
/// typically a closed chromosome plus one or two closed plasmids, and a
/// genome-wide rule leaves all of them open, so the origin of each replicon stops
/// being normalised.
///
/// Measured on E. coli K-12 split into a 4.0 Mb "chromosome" and a 0.64 Mb
/// "plasmid", both closed, compared against itself with each replicon's origin
/// rotated: with the topology declared, SCJ 0 and `observable_fraction` 1.0000;
/// without it, SCJ 4 (two moved seams, counted in both genomes) and 0.9993 over
/// 2,802 of 2,806 adjacencies. The junction count survives either way — a moved
/// seam lands on a contig end in the other genome, where the
/// positive-contradiction rule declines to call it broken — so the cost falls on
/// `scj_distance`, which is deliberately uncorrected, and on the power discount.
fn adjacency_set(series: &[Landmark], circular: &HashSet<u16>) -> HashMap<Adjacency, u64> {
    let mut adj = HashMap::with_capacity(series.len());
    for w in series.windows(2) {
        if w[0].contig != w[1].contig {
            continue; // consecutive in the file, but on different contigs
        }
        adj.insert(undirected(w[0].seq, w[1].seq), w[0].pos);
    }
    // `series` is sorted by (contig, pos), so each contig is one run.
    let mut i = 0usize;
    while i < series.len() {
        let contig = series[i].contig;
        let mut j = i;
        while j + 1 < series.len() && series[j + 1].contig == contig {
            j += 1;
        }
        if circular.contains(&contig) && j - i + 1 >= 3 {
            adj.insert(undirected(series[j].seq, series[i].seq), series[j].pos);
        }
        i = j + 1;
    }
    adj
}

fn contigs_of(series: &[Landmark]) -> HashSet<u16> {
    series.iter().map(|lm| lm.contig).collect()
}

/// The contigs each genome should close, for one pairwise comparison.
///
/// Two regimes, and the split is deliberate.
///
/// When **both** records carry declared topology, each genome closes exactly the
/// contigs its assembler called circular. This is the long-read case and the one
/// the per-contig closure exists for.
///
/// When **either** record is silent — every TGT written before the field existed
/// — both fall back to the original rule: close only when each genome is a single
/// contig. Deciding per genome instead would change legacy results, because a
/// closed genome compared against a draft would gain its seam adjacency and shift
/// `observable_fraction` from `1 − (K−1)/S` to `1 − K/S`. That is arguably the
/// truer number, but it is not a change to make silently on the strength of an
/// absent field.
fn circular_contigs(
    record_a: &TgtRecord,
    kept_a: &[Landmark],
    record_b: &TgtRecord,
    kept_b: &[Landmark],
) -> (HashSet<u16>, HashSet<u16>) {
    let declared = |record: &TgtRecord, series: &[Landmark]| -> HashSet<u16> {
        contigs_of(series)
            .into_iter()
            .filter(|&c| {
                // contig_id 0 means "unspecified / single contig"; 1+ indexes the table
                let idx = if c == 0 { 0 } else { c as usize - 1 };
                record.contig_circular.get(idx).copied().unwrap_or(false)
            })
            .collect()
    };

    if !record_a.contig_circular.is_empty() && !record_b.contig_circular.is_empty() {
        return (declared(record_a, kept_a), declared(record_b, kept_b));
    }

    let a = contigs_of(kept_a);
    let b = contigs_of(kept_b);
    if a.len() == 1 && b.len() == 1 {
        (a, b)
    } else {
        (HashSet::new(), HashSet::new())
    }
}

pub fn structural_synteny(record_a: &TgtRecord, record_b: &TgtRecord) -> Option<StructuralSynteny> {
    let series_a = raw_landmarks(record_a);
    let series_b = raw_landmarks(record_b);

    // A canonical sequence occurring at several loci cannot be assigned to one
    // of them, so it carries no usable order information and every copy
    // contributes a spurious adjacency. Repeats are dropped from BOTH genomes.
    //
    // Measured on E. coli K-12 (BcgI): 13 canonical sequences are multi-copy,
    // 63 tag instances or 2.1% of the total, one of them at 11 loci.
    let count = |v: &[Landmark]| -> HashMap<[u8; 32], usize> {
        let mut m = HashMap::new();
        for lm in v {
            *m.entry(lm.seq).or_insert(0) += 1;
        }
        m
    };
    let ca = count(&series_a);
    let cb = count(&series_b);

    let unique_a: HashSet<[u8; 32]> =
        ca.iter().filter(|(_, &n)| n == 1).map(|(s, _)| *s).collect();
    let unique_b: HashSet<[u8; 32]> =
        cb.iter().filter(|(_, &n)| n == 1).map(|(s, _)| *s).collect();
    let shared: HashSet<[u8; 32]> = unique_a.intersection(&unique_b).copied().collect();
    if shared.len() < 2 {
        return None;
    }

    // Keep only shared landmarks, preserving each genome's own order — and do it
    // *before* collapsing overlapping cut sites, so both genomes collapse the
    // same tag set and cannot disagree about run structure. See `collapse_runs`.
    let filtered_a: Vec<Landmark> =
        series_a.into_iter().filter(|l| shared.contains(&l.seq)).collect();
    let filtered_b: Vec<Landmark> =
        series_b.into_iter().filter(|l| shared.contains(&l.seq)).collect();
    let before_collapse = filtered_a.len() + filtered_b.len();

    let collapsed_a = collapse_runs(filtered_a);
    let collapsed_b = collapse_runs(filtered_b);
    let collapsed = before_collapse - (collapsed_a.len() + collapsed_b.len());

    // Indels shift positions, so a run that chains in one genome can fall apart
    // in the other and pick a different representative. Re-intersect to restore
    // the invariant every downstream step assumes: the two series carry exactly
    // the same set of sequences, and any difference between them is order.
    let surviving: HashSet<[u8; 32]> = collapsed_a
        .iter()
        .map(|l| l.seq)
        .collect::<HashSet<_>>()
        .intersection(&collapsed_b.iter().map(|l| l.seq).collect::<HashSet<_>>())
        .copied()
        .collect();
    if surviving.len() < 2 {
        return None;
    }
    let kept_a: Vec<Landmark> =
        collapsed_a.into_iter().filter(|l| surviving.contains(&l.seq)).collect();
    let kept_b: Vec<Landmark> =
        collapsed_b.into_iter().filter(|l| surviving.contains(&l.seq)).collect();

    // Closure is per contig, from the topology the assembler recorded.
    let (circ_a, circ_b) = circular_contigs(record_a, &kept_a, record_b, &kept_b);
    // Reported as a single flag for compatibility: true when every contig that
    // carries landmarks, in both genomes, is closed.
    let circular = !circ_a.is_empty()
        && circ_a.len() == contigs_of(&kept_a).len()
        && !circ_b.is_empty()
        && circ_b.len() == contigs_of(&kept_b).len();

    let adj_a = adjacency_set(&kept_a, &circ_a);
    let adj_b = adjacency_set(&kept_b, &circ_b);
    if adj_a.is_empty() && adj_b.is_empty() {
        return None;
    }

    let conserved = adj_a.keys().filter(|k| adj_b.contains_key(*k)).count();
    let union = adj_a.len() + adj_b.len() - conserved;

    // A junction is an adjacency of A that B *contradicts* — not merely one it
    // fails to show. The difference is the whole behaviour on draft assemblies:
    // each contig boundary in B hides one adjacency, so calling every absence a
    // junction would report one spurious junction per contig break, 99 of them
    // on a 100-contig draft.
    //
    // B contradicts {a, b} when it puts something else on both sides of a (or of
    // b): a landmark with two neighbours in B has no room left for the partner
    // it has in A. A landmark at a contig end has one neighbour, so the missing
    // adjacency is unobserved rather than broken, and is not counted.
    let mut degree_b: HashMap<[u8; 32], u8> = HashMap::with_capacity(adj_b.len());
    for (x, y) in adj_b.keys() {
        *degree_b.entry(*x).or_insert(0) += 1;
        *degree_b.entry(*y).or_insert(0) += 1;
    }
    let saturated = |seq: &[u8; 32]| degree_b.get(seq).copied().unwrap_or(0) >= 2;

    let mut junctions: Vec<u64> = adj_a
        .iter()
        .filter(|(k, _)| !adj_b.contains_key(*k))
        .filter(|((x, y), _)| saturated(x) || saturated(y))
        .map(|(_, &pos)| pos)
        .collect();
    junctions.sort_unstable();

    // How much of A's order B is in a position to judge at all. An adjacency
    // whose partners both sit at contig ends in B can never be contradicted, so
    // a junction there is invisible rather than absent. This is where contig
    // count belongs: as the denominator of a power statement, never as a divisor
    // of the junction count. The contamination from fragmentation is additive
    // (roughly one hidden adjacency per contig break), so dividing a count by
    // contig number shrinks the real signal as 1/n while the artifact term tends
    // to 1 -- the normalised statistic converges to the same value whether or not
    // the genome is rearranged.
    let observable = adj_a
        .keys()
        .filter(|(x, y)| saturated(x) || saturated(y))
        .count();

    // Orientation. The adjacency metric deliberately cannot see the *extent* of
    // an event: one inversion breaks exactly two adjacencies whether it spans
    // 5 kb or 500 kb. The orientation bit supplies what is missing, because
    // every landmark inside an inversion flips while every landmark outside it
    // does not. Counting flips therefore measures how much of the genome moved,
    // and pairs with the junction count, which measures how many times it moved.
    //
    // Palindromic tags are its blind spot: they equal their own reverse
    // complement, so they read the same in both orientations. They are reported
    // rather than silently dropped, since they bound the signal's resolution.
    let orient_a: HashMap<[u8; 32], Landmark> =
        kept_a.iter().map(|lm| (lm.seq, *lm)).collect();
    let mut orientation_mismatches = 0usize;
    let mut orientation_uninformative = 0usize;
    for lm_b in &kept_b {
        let Some(lm_a) = orient_a.get(&lm_b.seq) else { continue };
        if lm_a.palindromic || lm_b.palindromic {
            orientation_uninformative += 1;
        } else if lm_a.rc != lm_b.rc {
            orientation_mismatches += 1;
        }
    }
    let informative = kept_b.len().saturating_sub(orientation_uninformative);

    // Reverse-complementing a whole assembly flips every landmark, which is a
    // strand convention rather than biology. Score against the majority frame,
    // so "all flipped" reads as zero, exactly as "none flipped" does. The price
    // is a genuine identifiability limit: a genome inverted over more than half
    // its length is indistinguishable from its complement, and the fraction
    // saturates at 0.5. The junction count does not saturate, so the two
    // together still separate the cases.
    //
    // For comparison with alignment-based methods that use a fixed reference
    // (e.g. dnadiff), also report the raw mismatch fraction relative to the
    // first genome. This is orientation_mismatches / informative and ranges in
    // [0, 1] with no saturation.
    let minority = orientation_mismatches.min(informative - orientation_mismatches);

    Some(StructuralSynteny {
        score: conserved as f64 / union as f64,
        shared_tags: surviving.len(),
        repeats_dropped: (ca.len() - unique_a.len()) + (cb.len() - unique_b.len()),
        landmarks_collapsed: collapsed,
        circular,
        conserved_adjacencies: conserved,
        breakpoints: junctions.len(),
        scj_distance: adj_a.len() + adj_b.len() - 2 * conserved,
        // Normalised so genomes of different landmark density stay comparable.
        breakpoint_density: junctions.len() as f64 / surviving.len() as f64,
        junctions,
        orientation_uninformative,
        observable_adjacencies: observable,
        observable_fraction: if adj_a.is_empty() {
            0.0
        } else {
            observable as f64 / adj_a.len() as f64
        },
        inverted_fraction: if informative == 0 {
            0.0
        } else {
            minority as f64 / informative as f64
        },
        // Raw orientation mismatch fraction relative to genome A (the first
        // argument). This matches fixed-reference alignment methods like
        // dnadiff and does not saturate at 0.5.
        raw_inverted_fraction: if informative == 0 {
            0.0
        } else {
            orientation_mismatches as f64 / informative as f64
        },
        orientation_mismatches: minority,
        orientation_mismatches_raw: orientation_mismatches,
    })
}

/// Output of [`structural_synteny`].
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralSynteny {
    /// Conserved adjacencies over their union. 1.0 = collinear. Reported for
    /// continuity, but it divides a discrete event by the landmark count: one
    /// inversion among ~2,900 landmarks moves it by 0.0014. Prefer
    /// [`Self::breakpoints`].
    pub score: f64,
    /// Landmarks present and unique in both genomes — the basis of the score.
    pub shared_tags: usize,
    pub conserved_adjacencies: usize,
    /// Junctions: adjacencies of A that B contradicts by placing other
    /// landmarks on both sides. Exactly 2 per inversion and 3 per translocation,
    /// matching nucmer/SyRI and Syn2bANI's `breakpoint_count`. Zero under
    /// substitution loads up to at least 5%, and zero when B is merely a
    /// fragmented assembly of A.
    pub breakpoints: usize,
    /// Symmetric difference of the two adjacency sets: the single-cut-or-join
    /// distance as published, so it is left uncorrected and counts an adjacency
    /// hidden by a contig break the same as one genuinely broken. It is twice
    /// [`Self::breakpoints`] for two closed genomes, and larger than that for a
    /// draft assembly; prefer [`Self::breakpoints`] there.
    pub scj_distance: usize,
    /// Reference positions of the broken adjacencies. This is the output worth
    /// consuming: on E. coli K-12 a 400 kb inversion is localised to within one
    /// landmark spacing of the true junctions.
    pub junctions: Vec<u64>,
    /// Junctions per shared landmark, for comparing across landmark densities.
    pub breakpoint_density: f64,
    /// Multi-copy canonical sequences excluded, summed over both genomes. They
    /// are ambiguous for order and would each add a spurious adjacency.
    pub repeats_dropped: usize,
    /// Landmarks dropped as overlapping cut sites, summed over both genomes.
    pub landmarks_collapsed: usize,
    /// Whether both genomes were single-contig and the series were closed into
    /// cycles. When false the comparison is origin-sensitive at the ends.
    pub circular: bool,
    /// Adjacencies of A that B is in a position to contradict — both partners
    /// are not stranded at contig ends. Junctions can only ever be found here.
    pub observable_adjacencies: usize,
    /// [`Self::observable_adjacencies`] over all of A's adjacencies: the share of
    /// A's order that this comparison can judge. **This is where contig count
    /// belongs** — as a discount on detection power, never as a divisor of
    /// [`Self::breakpoints`]. Expect to recover roughly this fraction of the true
    /// junctions; 1.0 for two closed genomes.
    pub observable_fraction: f64,
    /// Shared landmarks whose stored window is reverse-complemented in one
    /// genome relative to the other: the landmarks that lie inside an inverted
    /// segment. Counted against the majority orientation, so a whole-genome
    /// reverse complement scores 0. Independent of [`Self::breakpoints`], which
    /// counts events rather than their extent.
    pub orientation_mismatches: usize,
    /// Raw count of orientation mismatches relative to genome A (the first
    /// argument), before scoring against the majority frame. This is the
    /// numerator of [`Self::raw_inverted_fraction`].
    pub orientation_mismatches_raw: usize,
    /// Shared landmarks that are their own reverse complement and so cannot
    /// report orientation. Reported because they bound the resolution of
    /// [`Self::orientation_mismatches`].
    pub orientation_uninformative: usize,
    /// [`Self::orientation_mismatches`] over the orientation-informative shared
    /// landmarks: the fraction of the shared genome that sits in inverted
    /// orientation. This is a fraction, so unlike the junction count it is
    /// directly comparable with alignment-based synteny measures. Saturates at
    /// 0.5, since past that point the minority frame becomes the majority one.
    pub inverted_fraction: f64,
    /// Orientation mismatch fraction relative to genome A (the first argument),
    /// i.e. [`Self::orientation_mismatches_raw`] / informative. This matches
    /// fixed-reference alignment methods such as dnadiff and ranges in [0, 1]
    /// without saturation. It is not invariant to whole-genome reverse
    /// complement: a genome and its reverse complement read as 1.0, just as
    /// dnadiff reports.
    pub raw_inverted_fraction: f64,
}

fn genome_order_in_graph(graph: &TagAdjacencyGraph, genome_id: &str) -> Option<Vec<u64>> {
    // Collect all (tag_id, position) pairs for this genome
    let mut tagged: Vec<(u64, u64)> = graph
        .nodes
        .values()
        .filter_map(|node| {
            node.positions
                .get(genome_id)
                .map(|&(pos, _)| (node.tag_id, pos))
        })
        .collect();

    if tagged.is_empty() {
        return None;
    }

    // Sort by genomic position
    tagged.sort_by_key(|&(_, pos)| pos);
    Some(tagged.into_iter().map(|(id, _)| id).collect())
}

// ─────────────────────────────────────────────────────────────────────────────
// Adjacency Jaccard
// ─────────────────────────────────────────────────────────────────────────────

/// Jaccard similarity of adjacent tag pairs between two TGT records.
///
/// For each record, builds a set of (tag_sequence_i, tag_sequence_{i+1}) pairs.
/// The Jaccard similarity is |A ∩ B| / |A ∪ B|.
///
/// This operates on tag *sequences* (not IDs), so it can compare records
/// that have not been inserted into a graph.
///
/// Returns a value in [0, 1] where 1.0 means identical adjacencies.
pub fn adjacency_jaccard(record_a: &TgtRecord, record_b: &TgtRecord) -> f64 {
    if record_a.adjacency_set.is_empty() && record_b.adjacency_set.is_empty() {
        return 0.0;
    }

    let intersection: HashSet<_> = record_a.adjacency_set.intersection(&record_b.adjacency_set).collect();
    let union: HashSet<_> = record_a.adjacency_set.union(&record_b.adjacency_set).collect();

    intersection.len() as f64 / union.len() as f64
}



/// Kendall's tau correlation on tag presence order between two records.
///
/// 1. Find common tags between the two records (by exact sequence match).
/// 2. For each common tag, record its rank (position index) in each genome.
/// 3. Compute Kendall's tau_b: (concordant_pairs - discordant_pairs) / total_pairs.
///
/// Returns a value in [-1, 1] where 1.0 means identical relative order,
/// -1.0 means completely reversed order, and 0.0 means uncorrelated.
pub fn kendall_tag_order(record_a: &TgtRecord, record_b: &TgtRecord) -> f64 {
    let rank_a = build_rank_map(&record_a.tags);

    // Collect common tags' ranks in A, ordered by B's tag order
    let mut common_ranks: Vec<usize> = Vec::new();
    for tag in &record_b.tags {
        if let Some(&rank_a) = rank_a.get(&tag.sequence) {
            common_ranks.push(rank_a);
        }
    }

    let n = common_ranks.len();
    if n < 2 {
        return 0.0;
    }

    let inversions = count_inversions(&mut common_ranks);
    let total_pairs = n * (n - 1) / 2;
    if total_pairs == 0 {
        return 0.0;
    }

    1.0 - 2.0 * inversions as f64 / total_pairs as f64
}

/// Count inversions in a slice using merge-sort (O(n log n)).
fn count_inversions(arr: &mut [usize]) -> usize {
    let n = arr.len();
    if n < 2 {
        return 0;
    }
    let mut temp = vec![0; n];
    merge_sort_count(arr, &mut temp, 0, n)
}

fn merge_sort_count(arr: &mut [usize], temp: &mut [usize], left: usize, right: usize) -> usize {
    if right - left <= 1 {
        return 0;
    }
    let mid = left + (right - left) / 2;
    let mut inv_count = merge_sort_count(arr, temp, left, mid);
    inv_count += merge_sort_count(arr, temp, mid, right);

    let mut i = left;
    let mut j = mid;
    let mut k = left;
    while i < mid && j < right {
        if arr[i] <= arr[j] {
            temp[k] = arr[i];
            i += 1;
        } else {
            temp[k] = arr[j];
            inv_count += mid - i;
            j += 1;
        }
        k += 1;
    }
    while i < mid {
        temp[k] = arr[i];
        i += 1;
        k += 1;
    }
    while j < right {
        temp[k] = arr[j];
        j += 1;
        k += 1;
    }
    arr[left..right].copy_from_slice(&temp[left..right]);
    inv_count
}

/// Build a map from tag sequence to its position index (rank) in the genome.
fn build_rank_map(tags: &[Tag]) -> HashMap<[u8; 32], usize> {
    tags.iter()
        .enumerate()
        .map(|(idx, tag)| (tag.sequence, idx))
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Breakpoint count
// ─────────────────────────────────────────────────────────────────────────────

/// Breakpoint count: number of positions where tag adjacency changes.
///
/// Two genomes share an adjacency if they have the same pair of consecutive
/// tags (by sequence). A breakpoint is a position in one genome where the
/// next tag differs from what the other genome has at the corresponding
/// location.
///
/// Algorithm:
/// 1. Find common tags between the two records.
/// 2. Extract the subsequence of common tags in each genome (preserving order).
/// 3. Count how many adjacent pairs in genome A are *not* adjacent in genome B.
///
/// This is a symmetric measure: breakpoint_count(A, B) == breakpoint_count(B, A).
pub fn breakpoint_count(record_a: &TgtRecord, record_b: &TgtRecord) -> usize {
    // Get ranks of common tags in each genome
    let rank_a = build_rank_map(&record_a.tags);
    let rank_b = build_rank_map(&record_b.tags);

    let common_seqs: Vec<[u8; 32]> = rank_a
        .keys()
        .filter(|seq| rank_b.contains_key(*seq))
        .copied()
        .collect();

    if common_seqs.len() < 2 {
        return 0;
    }

    // Build adjacency sets of common tags for each genome
    let mut adj_a: HashSet<([u8; 32], [u8; 32])> = HashSet::new();
    let mut adj_b: HashSet<([u8; 32], [u8; 32])> = HashSet::new();

    // For genome A: scan through tags and record adjacencies between common tags
    for w in record_a.tags.windows(2) {
        let (s1, s2) = (w[0].sequence, w[1].sequence);
        if rank_a.contains_key(&s1)
            && rank_a.contains_key(&s2)
            && rank_b.contains_key(&s1)
            && rank_b.contains_key(&s2)
        {
            let key = normalize_pair(s1, s2);
            adj_a.insert(key);
        }
    }

    for w in record_b.tags.windows(2) {
        let (s1, s2) = (w[0].sequence, w[1].sequence);
        if rank_a.contains_key(&s1)
            && rank_a.contains_key(&s2)
            && rank_b.contains_key(&s1)
            && rank_b.contains_key(&s2)
        {
            let key = normalize_pair(s1, s2);
            adj_b.insert(key);
        }
    }

    // Breakpoints = adjacencies present in one but not the other
    let symmetric_diff: HashSet<_> = adj_a.symmetric_difference(&adj_b).collect();
    symmetric_diff.len()
}

/// Normalize a tag pair for undirected comparison.
fn normalize_pair(a: [u8; 32], b: [u8; 32]) -> ([u8; 32], [u8; 32]) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzyme::enzyme::EnzymeType;
    use crate::tgt::record::TgtRecord;
    use crate::tgt::tag::{Strand, Tag};
    use crate::synteny::graph::TagAdjacencyGraph;

    /// Build a record from explicit 32 bp sequences at 1 kb spacing.
    fn record_from(id: &str, seqs: &[String]) -> TgtRecord {
        let mut r = TgtRecord::new(id, (seqs.len() as u64 + 1) * 1000);
        for (i, s) in seqs.iter().enumerate() {
            let mut buf = [0u8; 32];
            let b = s.as_bytes();
            buf[..b.len()].copy_from_slice(b);
            r.add_tag(Tag::new(buf, (i as u64 + 1) * 1000, EnzymeType::BcgI, Strand::Forward, 0));
        }
        r
    }

    /// Deterministic distinct 32-mers, so tests do not depend on an RNG.
    fn seqs(n: usize) -> Vec<String> {
        (0..n)
            .map(|i| {
                let mut v = i as u64 + 1;
                (0..32)
                    .map(|_| {
                        v = v.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                        match (v >> 33) % 4 {
                            0 => 'A',
                            1 => 'C',
                            2 => 'G',
                            _ => 'T',
                        }
                    })
                    .collect()
            })
            .collect()
    }

    fn revcomp(s: &str) -> String {
        s.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                o => o,
            })
            .collect()
    }

    #[test]
    fn canonical_sequence_agrees_across_strands() {
        let s = seqs(1).remove(0);
        let a = record_from("a", &[s.clone()]);
        let b = record_from("b", &[revcomp(&s)]);
        assert_eq!(
            a.tags[0].canonical_sequence(),
            b.tags[0].canonical_sequence(),
            "a tag and its reverse complement must share one identity"
        );
    }

    #[test]
    fn structural_synteny_is_one_for_collinear_genomes() {
        let s = seqs(40);
        let a = record_from("a", &s);
        let b = record_from("b", &s);
        let r = structural_synteny(&a, &b).expect("40 shared tags");
        assert_eq!(r.score, 1.0);
        assert_eq!(r.breakpoints, 0);
        assert_eq!(r.shared_tags, 40);
    }

    /// The property the old metric fails. Losing tags to substitutions must not
    /// move a structural score, because presence is a divergence signal and only
    /// order is a structural one.
    #[test]
    fn structural_synteny_is_invariant_to_tag_loss() {
        let s = seqs(40);
        let a = record_from("a", &s);
        // Genome B has lost every third tag, as substitutions in the recognition
        // site or tag body would cause. Order is otherwise untouched.
        let kept: Vec<String> = s
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 3 != 0)
            .map(|(_, x)| x.clone())
            .collect();
        let b = record_from("b", &kept);

        let r = structural_synteny(&a, &b).expect("shared tags remain");
        assert_eq!(
            r.score, 1.0,
            "losing a third of the tags must leave a structural score untouched"
        );
        assert_eq!(r.breakpoints, 0);
        assert_eq!(r.shared_tags, kept.len());

        // Contrast: the legacy adjacency metric collapses on the same input.
        let legacy = adjacency_jaccard(&a, &b);
        assert!(
            legacy < 0.4,
            "legacy adjacency_jaccard should collapse here, got {legacy}"
        );
    }

    /// An inversion reverse-complements its segment and reverses tag order.
    /// Canonicalisation keeps the tags matchable; ordered adjacency sees the
    /// reversal.
    #[test]
    fn structural_synteny_detects_an_inversion() {
        let s = seqs(40);
        let a = record_from("a", &s);

        let mut inv: Vec<String> = s.clone();
        inv[10..30].reverse();
        for x in inv[10..30].iter_mut() {
            *x = revcomp(x);
        }
        let b = record_from("b", &inv);

        let r = structural_synteny(&a, &b).expect("all tags still shared");
        assert_eq!(
            r.shared_tags, 40,
            "canonicalisation must keep inverted tags in the shared set"
        );
        assert!(
            r.score < 0.95,
            "a 20-tag inversion must be visible, got score {}",
            r.score
        );
        assert!(r.breakpoints >= 2, "got {} breakpoints", r.breakpoints);
    }

    // ── Invariance controls ────────────────────────────────────────────────
    //
    // Each of these is a transformation that changes nothing biological, so
    // every one must report exactly zero junctions and zero inverted fraction.
    // Both this crate and Syn2bANI have shipped metrics that failed at least
    // one of them, so they are asserted rather than assumed.

    #[test]
    fn invariance_self_comparison() {
        let s = seqs(40);
        let r = structural_synteny(&record_from("a", &s), &record_from("b", &s))
            .expect("identical genomes share every tag");
        assert_eq!(r.breakpoints, 0);
        assert_eq!(r.scj_distance, 0);
        assert_eq!(r.inverted_fraction, 0.0);
        assert_eq!(r.score, 1.0);
    }

    #[test]
    fn invariance_whole_genome_reverse_complement() {
        // A chromosome read off the other strand is the same chromosome. Every
        // landmark flips, which is exactly why the fraction is scored against
        // the majority orientation.
        let s = seqs(40);
        let mut rc: Vec<String> = s.iter().map(|x| revcomp(x)).collect();
        rc.reverse();

        let r = structural_synteny(&record_from("a", &s), &record_from("b", &rc))
            .expect("canonicalisation must keep every tag shared");
        assert_eq!(r.shared_tags, 40, "reverse complement must not lose tags");
        assert_eq!(r.breakpoints, 0, "got {} junctions", r.breakpoints);
        assert_eq!(
            r.inverted_fraction, 0.0,
            "a strand convention is not a rearrangement"
        );
    }

    #[test]
    fn invariance_circular_origin_rotation() {
        // Assemblies of a circular chromosome start at an arbitrary base.
        let s = seqs(40);
        let mut rotated = s[13..].to_vec();
        rotated.extend_from_slice(&s[..13]);

        let r = structural_synteny(&record_from("a", &s), &record_from("b", &rotated))
            .expect("rotation keeps every tag");
        assert!(r.circular, "single-contig genomes must be closed into cycles");
        assert_eq!(r.breakpoints, 0, "got {} junctions", r.breakpoints);
    }

    #[test]
    fn invariance_closed_multi_replicon_origin_rotation() {
        // A long-read assembly is normally a closed chromosome plus closed
        // plasmids. Rotating each replicon's origin must change nothing. A
        // genome-wide rule sees two contigs and declines to close either, which
        // leaves the moved seams in the symmetric difference: on real E. coli
        // split into two closed replicons that costs SCJ 4 and drops
        // observable_fraction to 0.9993. Junctions are protected either way by
        // the positive-contradiction rule, so SCJ is what this asserts on.
        let s = seqs(40);
        let chrom: Vec<String> = s[..28].to_vec();
        let plasmid: Vec<String> = s[28..].to_vec();

        let build = |id: &str, chrom_off: usize, plasmid_off: usize| {
            let mut r = TgtRecord::new(id, 60_000);
            r.contig_names = vec!["chromosome".into(), "plasmid".into()];
            r.contig_circular = vec![true, true];
            let mut push = |r: &mut TgtRecord, seqs: &[String], off: usize, contig: u16| {
                for i in 0..seqs.len() {
                    let seq = &seqs[(i + off) % seqs.len()];
                    let mut buf = [0u8; 32];
                    buf[..seq.len()].copy_from_slice(seq.as_bytes());
                    r.add_tag(Tag::new(
                        buf,
                        (i as u64 + 1) * 1000,
                        EnzymeType::BcgI,
                        Strand::Forward,
                        contig,
                    ));
                }
            };
            push(&mut r, &chrom, chrom_off, 1);
            push(&mut r, &plasmid, plasmid_off, 2);
            r
        };

        let r = structural_synteny(&build("a", 0, 0), &build("b", 9, 4))
            .expect("rotation keeps every tag");
        assert!(r.circular, "both replicons are declared circular");
        assert_eq!(r.breakpoints, 0, "got {} junction(s)", r.breakpoints);
        assert_eq!(
            r.scj_distance, 0,
            "rotating a closed replicon must not move its seam, got {}",
            r.scj_distance
        );
        assert!(
            (r.observable_fraction - 1.0).abs() < 1e-9,
            "closed replicons hide no adjacency, got {}",
            r.observable_fraction
        );

        // Without the declaration the same pair falls back to the genome-wide
        // rule, and the seams reappear. This is the defect, held in place.
        let (mut a, mut b) = (build("a", 0, 0), build("b", 9, 4));
        a.contig_circular.clear();
        b.contig_circular.clear();
        let undeclared = structural_synteny(&a, &b).expect("same tags");
        assert!(!undeclared.circular);
        assert!(
            undeclared.scj_distance > 0,
            "undeclared topology should leave the seams in the symmetric difference"
        );
    }

    #[test]
    fn unshared_tag_does_not_split_a_collapse_run() {
        // Overlapping cut sites are collapsed to one representative per run, and
        // the run boundaries are found by chaining against the previous landmark.
        // A tag lost to a substitution or a sequencing error can therefore split
        // one run into two, yielding two representatives where the other genome
        // yields one — an extra landmark, an extra adjacency, and a junction that
        // no rearrangement produced. This is the `sub_2` residual.
        //
        // Collapsing after the shared-tag restriction removes the mechanism.
        // Checked for every choice of which run member is missing, since the
        // damaging case is the one where the lost tag was the representative.
        let s = seqs(12);
        let run_positions = [10_000u64, 10_020, 10_040]; // chained under 40 bp
        let solo_positions = [1_000u64, 2_000, 3_000, 20_000, 21_000];

        let build = |id: &str, skip: Option<usize>| {
            let mut r = TgtRecord::new(id, 40_000);
            let mut add = |r: &mut TgtRecord, seq: &str, pos: u64| {
                let mut buf = [0u8; 32];
                buf[..seq.len()].copy_from_slice(seq.as_bytes());
                r.add_tag(Tag::new(buf, pos, EnzymeType::BcgI, Strand::Forward, 1));
            };
            for (i, pos) in solo_positions.iter().take(3).enumerate() {
                add(&mut r, &s[i], *pos);
            }
            for (i, pos) in run_positions.iter().enumerate() {
                if Some(i) == skip {
                    continue; // the tag a substitution destroyed
                }
                add(&mut r, &s[3 + i], *pos);
            }
            for (i, pos) in solo_positions.iter().skip(3).enumerate() {
                add(&mut r, &s[6 + i], *pos);
            }
            r
        };

        let a = build("a", None);
        for missing in 0..run_positions.len() {
            let b = build("b", Some(missing));
            let r = structural_synteny(&a, &b)
                .unwrap_or_else(|| panic!("missing={missing} left too few shared tags"));
            assert_eq!(
                r.breakpoints, 0,
                "losing run member {missing} manufactured {} junction(s)",
                r.breakpoints
            );
        }
    }

    #[test]
    fn invariance_fragmented_assembly() {
        // Splitting one contig into three must not invent junctions across the
        // breaks, and must not close a cycle that no longer has an origin.
        let s = seqs(40);
        let a = record_from("a", &s);

        let mut b = TgtRecord::new("b", 41_000);
        for (i, seq) in s.iter().enumerate() {
            let mut buf = [0u8; 32];
            buf[..seq.len()].copy_from_slice(seq.as_bytes());
            let contig = (i / 14 + 1) as u16;
            b.add_tag(Tag::new(
                buf,
                (i as u64 + 1) * 1000,
                EnzymeType::BcgI,
                Strand::Forward,
                contig,
            ));
        }
        b.contig_names = vec!["c1".into(), "c2".into(), "c3".into()];

        let r = structural_synteny(&a, &b).expect("fragmentation keeps every tag");
        assert!(!r.circular, "a multi-contig genome has no origin to normalise");
        assert_eq!(
            r.breakpoints, 0,
            "contig boundaries are missing adjacencies, not broken ones; got {:?}",
            r.junctions
        );
        // Three contigs hide two adjacencies out of 39, and that is what
        // observable_fraction is for. Contig count is a discount on power, never
        // a divisor of the junction count.
        assert_eq!(r.observable_adjacencies, 37);
        assert!(
            (r.observable_fraction - 37.0 / 39.0).abs() < 1e-9,
            "got {}",
            r.observable_fraction
        );
    }

    #[test]
    fn observable_fraction_is_one_for_two_closed_genomes() {
        let s = seqs(40);
        let r = structural_synteny(&record_from("a", &s), &record_from("b", &s))
            .expect("shared");
        assert_eq!(r.observable_fraction, 1.0);
    }

    // ── Orientation signal ─────────────────────────────────────────────────

    #[test]
    fn orientation_measures_the_extent_of_an_inversion() {
        // The junction count is blind to how large an event is; the orientation
        // signal is exactly the missing axis. Two inversions of very different
        // size must give the same junction count and different fractions.
        let s = seqs(100);

        let invert = |from: usize, to: usize| -> Vec<String> {
            let mut v = s.clone();
            v[from..to].reverse();
            for x in v[from..to].iter_mut() {
                *x = revcomp(x);
            }
            v
        };

        let a = record_from("a", &s);
        let small = structural_synteny(&a, &record_from("small", &invert(20, 30)))
            .expect("shared");
        let large = structural_synteny(&a, &record_from("large", &invert(20, 70)))
            .expect("shared");

        assert_eq!(small.breakpoints, 2, "one inversion is two junctions");
        assert_eq!(large.breakpoints, 2, "size must not change the count");

        assert_eq!(small.orientation_mismatches, 10, "10 landmarks moved");
        assert_eq!(large.orientation_mismatches, 50, "50 landmarks moved");
        assert!((small.inverted_fraction - 0.10).abs() < 1e-9);
        assert!((large.inverted_fraction - 0.50).abs() < 1e-9);

        // raw_inverted_fraction is the same as inverted_fraction when the
        // minority frame is genome B's, but does not saturate past 50%.
        assert_eq!(small.orientation_mismatches_raw, 10);
        assert_eq!(large.orientation_mismatches_raw, 50);
        assert!((small.raw_inverted_fraction - 0.10).abs() < 1e-9);
        assert!((large.raw_inverted_fraction - 0.50).abs() < 1e-9);
    }

    #[test]
    fn raw_inverted_fraction_matches_fixed_reference() {
        // A fixed-reference metric should report 1.0 when genome B is the
        // reverse complement of genome A, exactly as dnadiff would.
        let s = seqs(40);
        let mut rc: Vec<String> = s.iter().map(|x| revcomp(x)).collect();
        rc.reverse();

        let r = structural_synteny(&record_from("a", &s), &record_from("b", &rc))
            .expect("canonicalisation must keep every tag shared");
        assert!((r.raw_inverted_fraction - 1.0).abs() < 1e-9,
            "whole-genome reverse complement should read as 1.0 for fixed-reference metric");
    }

    #[test]
    fn orientation_sees_an_inversion_too_small_for_a_junction() {
        // A single landmark inside an inversion leaves both flanking
        // adjacencies intact, so the junction count reads zero. The orientation
        // bit still flips, which is why the two signals have different floors.
        let s = seqs(40);
        let mut v = s.clone();
        v[17] = revcomp(&v[17]);

        let r = structural_synteny(&record_from("a", &s), &record_from("b", &v))
            .expect("shared");
        assert_eq!(r.breakpoints, 0, "below the junction resolution limit");
        assert_eq!(r.orientation_mismatches, 1, "but the orientation bit flips");
    }

    /// The two signals must be separable: the inversion's cost must not depend on
    /// how many tags substitutions removed.
    #[test]
    fn inversion_signal_survives_tag_loss() {
        let s = seqs(60);
        let a = record_from("a", &s);

        let mut inv: Vec<String> = s.clone();
        inv[20..40].reverse();
        for x in inv[20..40].iter_mut() {
            *x = revcomp(x);
        }

        let clean = structural_synteny(&a, &record_from("b", &inv)).unwrap();

        // Now drop a quarter of the tags from the inverted genome as well.
        let lossy: Vec<String> = inv
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 4 != 0)
            .map(|(_, x)| x.clone())
            .collect();
        let degraded = structural_synteny(&a, &record_from("b", &lossy)).unwrap();

        assert!(
            (clean.score - degraded.score).abs() < 0.25,
            "the inversion signal must not be washed out by tag loss: \
             clean {} vs degraded {}",
            clean.score,
            degraded.score
        );
    }

    #[test]
    fn structural_synteny_declines_when_nothing_is_shared() {
        let a = record_from("a", &seqs(10));
        let b = record_from("b", &seqs(10).iter().map(|s| revcomp(s)).collect::<Vec<_>>());
        // Reverse complements canonicalise to the same identity, so this pair IS
        // shared; use genuinely different sequences instead.
        assert!(structural_synteny(&a, &b).is_some());

        let c = record_from("c", &["A".repeat(32)]);
        let d = record_from("d", &["C".repeat(32)]);
        assert!(
            structural_synteny(&c, &d).is_none(),
            "fewer than two shared tags leaves no adjacency to compare"
        );
    }

    /// Helper: create a tag with index-encoded sequence
    fn make_tag(idx: u8, position: u64, enzyme: EnzymeType, strand: Strand) -> Tag {
        let mut seq = [b'A'; 32];
        seq[0] = idx;
        seq[1] = idx.wrapping_add(1);
        Tag::new(seq, position, enzyme, strand, 0)
    }

    /// Helper: create a record with evenly spaced tags
    fn make_record(genome_id: &str, tag_count: usize) -> TgtRecord {
        let mut record = TgtRecord::new(genome_id, 1_000_000);
        for i in 0..tag_count {
            let tag = make_tag(i as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward);
            record.add_tag(tag);
        }
        record
    }

    // ─────────────────────────────────────────────────────────────────────────
    // synteny_score
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_synteny_score_perfect() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 10);
        let record_b = make_record("genome_a", 10);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let paths = graph.linear_paths();
        assert!(!paths.is_empty(), "Should have at least one path");

        // With 2 identical genomes, edges have weight 2 (out of max 2)
        // Score should be close to 1.0
        let score = synteny_score(&paths[0], &graph);
        assert!(
            score >= 0.9 && score <= 1.0,
            "Two identical genomes should have near-perfect synteny score, got {}",
            score
        );
    }

    #[test]
    fn test_synteny_score_empty_path() {
        let graph = TagAdjacencyGraph::new();
        let score = synteny_score(&[], &graph);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_synteny_score_singleton_path() {
        let graph = TagAdjacencyGraph::new();
        let score = synteny_score(&[42], &graph);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_synteny_score_single_genome() {
        let mut graph = TagAdjacencyGraph::new();
        let record = make_record("genome_a", 10);
        graph.add_genome("genome_a", &record);
        graph.build_edges();
        graph.simplify(1);

        let paths = graph.linear_paths();
        // With only 1 genome, edge weight = 1, max = 1, ratio = 1.0
        let score = synteny_score(&paths[0], &graph);
        assert!(
            score >= 0.9 && score <= 1.0,
            "Single genome should have perfect score, got {}",
            score
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // adjacency_jaccard
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_adjacency_jaccard_identical() {
        let record_a = make_record("a", 5);
        let record_b = make_record("a", 5);
        let jaccard = adjacency_jaccard(&record_a, &record_b);
        assert!(
            (jaccard - 1.0).abs() < 1e-9,
            "Identical records should have Jaccard = 1.0, got {}",
            jaccard
        );
    }

    #[test]
    fn test_adjacency_jaccard_completely_different() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags 0, 1, 2
        for i in 0..3 {
            record_a.add_tag(make_tag(i, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: tags 10, 11, 12 (no overlap with A)
        for i in 10..13usize {
            record_b.add_tag(make_tag(i as u8, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        let jaccard = adjacency_jaccard(&record_a, &record_b);
        assert_eq!(jaccard, 0.0, "No common tags means Jaccard = 0");
    }

    #[test]
    fn test_adjacency_jaccard_partial() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags 0, 1, 2, 3 (adjacencies: (0,1), (1,2), (2,3))
        for i in 0..4usize {
            record_a.add_tag(make_tag(i as u8, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: tags 0, 1, 2, 4 (shared adjacencies: (0,1), (1,2); unique: (2,4))
        record_b.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(1, 200, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(2, 300, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(4, 400, EnzymeType::BcgI, Strand::Forward));

        // intersection = {(0,1), (1,2)} → 2
        // union = {(0,1), (1,2), (2,3), (2,4)} → 4
        // Jaccard = 2/4 = 0.5
        let jaccard = adjacency_jaccard(&record_a, &record_b);
        assert!(
            (jaccard - 0.5).abs() < 1e-9,
            "Partial overlap should give Jaccard = 0.5, got {}",
            jaccard
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // kendall_tag_order
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_kendall_identical_order() {
        let record_a = make_record("a", 5);
        let record_b = make_record("a", 5);
        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            (tau - 1.0).abs() < 1e-9,
            "Identical order should give tau = 1.0, got {}",
            tau
        );
    }

    #[test]
    fn test_kendall_reversed_order() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags in order 0, 1, 2, 3
        for i in 0..4u8 {
            record_a.add_tag(make_tag(i, (i as u64) * 100, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: same tags in reversed order 3, 2, 1, 0
        for i in (0..4u8).rev() {
            record_b.add_tag(make_tag(i, ((3 - i) as u64) * 100, EnzymeType::BcgI, Strand::Forward));
        }

        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            (tau - (-1.0)).abs() < 1e-9,
            "Reversed order should give tau = -1.0, got {}",
            tau
        );
    }

    #[test]
    fn test_kendall_no_common_tags() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        record_a.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(10, 100, EnzymeType::BcgI, Strand::Forward));

        let tau = kendall_tag_order(&record_a, &record_b);
        assert_eq!(tau, 0.0, "No common tags should give tau = 0");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // breakpoint_count
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_breakpoint_count_identical() {
        let record_a = make_record("a", 5);
        let record_b = make_record("a", 5);
        let bp = breakpoint_count(&record_a, &record_b);
        assert_eq!(bp, 0, "Identical genomes should have 0 breakpoints");
    }

    #[test]
    fn test_breakpoint_count_with_rearrangement() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);

        // Genome A: tags 0, 1, 2, 3, 4 (adjacencies: 0-1, 1-2, 2-3, 3-4)
        for i in 0..5usize {
            record_a.add_tag(make_tag(i as u8, (i * 100) as u64, EnzymeType::BcgI, Strand::Forward));
        }

        // Genome B: tags 0, 1, 3, 2, 4 (breakpoint at 1-3 and 3-2 vs 2-3)
        // Adjacencies in B: 0-1, 1-3, 3-2, 2-4
        // Common adjacencies: 0-1
        // Breakpoints: (1-2 in A but not B), (2-3 in A, 3-2 in B = same undirected),
        //              (3-4 in A but not B)
        record_b.add_tag(make_tag(0, 100, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(1, 200, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(3, 300, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(2, 400, EnzymeType::BcgI, Strand::Forward));
        record_b.add_tag(make_tag(4, 500, EnzymeType::BcgI, Strand::Forward));

        let bp = breakpoint_count(&record_a, &record_b);
        assert!(bp > 0, "Rearranged genome should have breakpoints, got {}", bp);
    }

    #[test]
    fn test_breakpoint_count_empty() {
        let record_a = TgtRecord::new("a", 1_000_000);
        let record_b = TgtRecord::new("b", 1_000_000);
        let bp = breakpoint_count(&record_a, &record_b);
        assert_eq!(bp, 0, "Empty records should have 0 breakpoints");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // pairwise_synteny_matrix
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_pairwise_matrix_identical_genomes() {
        let mut graph = TagAdjacencyGraph::new();
        let record_a = make_record("genome_a", 10);
        let record_b = make_record("genome_a", 10);

        graph.add_genome("genome_a", &record_a);
        graph.add_genome("genome_b", &record_b);
        graph.build_edges();
        graph.simplify(1);

        let matrix = pairwise_synteny_matrix(&graph);
        let score = matrix
            .get(&("genome_a".to_string(), "genome_b".to_string()))
            .copied()
            .unwrap_or(0.0);
        assert!(
            score >= 0.9,
            "Identical genomes should have high pairwise score, got {}",
            score
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // count_inversions / kendall_tag_order
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_count_inversions_sorted() {
        let mut arr = vec![1, 2, 3, 4, 5];
        assert_eq!(count_inversions(&mut arr), 0, "Sorted array should have 0 inversions");
    }

    #[test]
    fn test_count_inversions_reversed() {
        let mut arr = vec![5, 4, 3, 2, 1];
        assert_eq!(count_inversions(&mut arr), 10, "Reversed array of 5 should have 10 inversions");
    }

    #[test]
    fn test_kendall_tag_order_identical() {
        let record_a = make_record("a", 10);
        let record_b = make_record("a", 10);
        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            tau >= 0.99,
            "Identical genomes should have tau ≈ 1.0, got {}",
            tau
        );
    }

    #[test]
    fn test_kendall_tag_order_reversed() {
        let mut record_a = TgtRecord::new("a", 1_000_000);
        let mut record_b = TgtRecord::new("b", 1_000_000);
        for i in 0..10 {
            record_a.add_tag(make_tag(i as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward));
        }
        for i in 0..10 {
            record_b.add_tag(make_tag((9 - i) as u8, (i * 1000) as u64, EnzymeType::BcgI, Strand::Forward));
        }
        let tau = kendall_tag_order(&record_a, &record_b);
        assert!(
            tau <= -0.7,
            "Reversed genomes should have strongly negative tau, got {}",
            tau
        );
    }
}

/// Compute Pearson correlation of shared tag positions between two genomes.
/// Returns absolute value (inversions produce negative correlation but still syntenic).
pub fn position_correlation(record_a: &TgtRecord, record_b: &TgtRecord) -> f64 {
    let mut pos_a = Vec::new();
    let mut pos_b = Vec::new();

    // Build sequence->position map for record_b
    let pos_map_b: std::collections::HashMap<[u8; 32], u64> = record_b
        .tags
        .iter()
        .map(|t| (t.sequence, t.position))
        .collect();

    // Find shared tags and collect positions
    for tag_a in &record_a.tags {
        if let Some(&pos_b_val) = pos_map_b.get(&tag_a.sequence) {
            pos_a.push(tag_a.position as f64);
            pos_b.push(pos_b_val as f64);
        }
    }

    let n = pos_a.len();
    if n < 2 {
        return 0.0;
    }

    let mean_a: f64 = pos_a.iter().sum::<f64>() / n as f64;
    let mean_b: f64 = pos_b.iter().sum::<f64>() / n as f64;

    let cov: f64 = pos_a.iter().zip(pos_b.iter())
        .map(|(a, b)| (a - mean_a) * (b - mean_b))
        .sum();
    let var_a: f64 = pos_a.iter().map(|a| (a - mean_a).powi(2)).sum();
    let var_b: f64 = pos_b.iter().map(|b| (b - mean_b).powi(2)).sum();

    if var_a <= 0.0 || var_b <= 0.0 {
        return 0.0;
    }

    let corr = cov / (var_a.sqrt() * var_b.sqrt());
    corr.abs()
}

/// Compute windowed position correlation: split shared tags into windows
/// and compute the maximum Pearson correlation across windows.
/// Captures local synteny even when global shared ratio is low.
pub fn windowed_position_correlation(record_a: &TgtRecord, record_b: &TgtRecord, windows: usize) -> f64 {
    let mut pos_a = Vec::new();
    let mut pos_b = Vec::new();

    let pos_map_b: std::collections::HashMap<[u8; 32], u64> = record_b
        .tags
        .iter()
        .map(|t| (t.sequence, t.position))
        .collect();

    for tag_a in &record_a.tags {
        if let Some(&pos_b_val) = pos_map_b.get(&tag_a.sequence) {
            pos_a.push(tag_a.position);
            pos_b.push(pos_b_val);
        }
    }

    let n = pos_a.len();
    if n < windows * 2 {
        return position_correlation(record_a, record_b);
    }

    // Sort by position in genome A
    let mut pairs: Vec<(u64, u64)> = pos_a.into_iter().zip(pos_b.into_iter()).collect();
    pairs.sort_by_key(|(a, _)| *a);

    let window_size = n / windows;
    let mut max_corr = 0.0f64;

    for w in 0..windows {
        let start = w * window_size;
        let end = if w == windows - 1 { n } else { (w + 1) * window_size };
        
        let wa: Vec<f64> = pairs[start..end].iter().map(|(a, _)| *a as f64).collect();
        let wb: Vec<f64> = pairs[start..end].iter().map(|(_, b)| *b as f64).collect();

        let n_w = wa.len();
        if n_w < 2 {
            continue;
        }

        let mean_a: f64 = wa.iter().sum::<f64>() / n_w as f64;
        let mean_b: f64 = wb.iter().sum::<f64>() / n_w as f64;

        let cov: f64 = wa.iter().zip(wb.iter())
            .map(|(a, b)| (a - mean_a) * (b - mean_b))
            .sum();
        let var_a: f64 = wa.iter().map(|a| (a - mean_a).powi(2)).sum();
        let var_b: f64 = wb.iter().map(|b| (b - mean_b).powi(2)).sum();

        if var_a > 0.0 && var_b > 0.0 {
            let corr = (cov / (var_a.sqrt() * var_b.sqrt())).abs();
            if corr > max_corr {
                max_corr = corr;
            }
        }
    }

    max_corr
}

/// Compute tag density correlation across genomic windows.
/// Divides each genome into `windows` bins and computes the Pearson correlation
/// of tag counts per bin. Captures synteny even when exact tag order is disrupted
/// by rearrangements, because syntenic regions share similar tag density patterns.
pub fn windowed_density_correlation(record_a: &TgtRecord, record_b: &TgtRecord, windows: usize) -> f64 {
    if record_a.total_length == 0 || record_b.total_length == 0 || windows == 0 {
        return 0.0;
    }

    let bin_size_a = (record_a.total_length as f64 / windows as f64).ceil() as u64;
    let bin_size_b = (record_b.total_length as f64 / windows as f64).ceil() as u64;

    let mut bins_a = vec![0.0f64; windows];
    for tag in &record_a.tags {
        let bin = (tag.position / bin_size_a).min((windows - 1) as u64) as usize;
        bins_a[bin] += 1.0;
    }

    let mut bins_b = vec![0.0f64; windows];
    for tag in &record_b.tags {
        let bin = (tag.position / bin_size_b).min((windows - 1) as u64) as usize;
        bins_b[bin] += 1.0;
    }

    // Compute Pearson correlation
    let n = windows as f64;
    let mean_a: f64 = bins_a.iter().sum::<f64>() / n;
    let mean_b: f64 = bins_b.iter().sum::<f64>() / n;

    if mean_a == 0.0 || mean_b == 0.0 {
        return 0.0;
    }

    let cov: f64 = bins_a.iter().zip(bins_b.iter())
        .map(|(a, b)| (a - mean_a) * (b - mean_b))
        .sum();
    let var_a: f64 = bins_a.iter().map(|a| (a - mean_a).powi(2)).sum();
    let var_b: f64 = bins_b.iter().map(|b| (b - mean_b).powi(2)).sum();

    if var_a <= 0.0 || var_b <= 0.0 {
        return 0.0;
    }

    let corr = cov / (var_a.sqrt() * var_b.sqrt());
    corr.abs()
}

// ─────────────────────────────────────────────────────────────────────────────
// Approximate tag matching (Hamming distance ≤ 1)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute Hamming distance between two 32-byte tag sequences.
/// Only counts non-null bytes (stops at first null byte in either sequence).
pub fn hamming_distance(a: &[u8; 32], b: &[u8; 32]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Count approximate common tags between two records, allowing Hamming distance ≤ max_dist.
///
/// For max_dist=1, uses a fast variant-generation approach:
/// For each tag in r1, generate all 1bp variants and check against r2's exact set.
/// Complexity: O(n * 93) for 32bp tags, where n = r1.tag_count().
pub fn count_approximate_common_tags(r1: &TgtRecord, r2: &TgtRecord, max_dist: usize) -> usize {
    if max_dist == 0 {
        return r1.sequence_set.intersection(&r2.sequence_set).count();
    }

    if max_dist == 1 {
        return count_approximate_common_tags_h1(r1, r2);
    }

    // For max_dist > 1, fall back to brute-force (slow for large datasets)
    let mut count = 0;
    for seq1 in &r1.sequence_set {
        for seq2 in &r2.sequence_set {
            if hamming_distance(seq1, seq2) <= max_dist {
                count += 1;
                break;
            }
        }
    }
    count
}

/// Fast Hamming-1 approximate matching using variant generation.
fn count_approximate_common_tags_h1(r1: &TgtRecord, r2: &TgtRecord) -> usize {
    let bases = [b'A', b'T', b'C', b'G'];
    let mut count = 0;

    for seq1 in &r1.sequence_set {
        // Check exact match first
        if r2.sequence_set.contains(seq1) {
            count += 1;
            continue;
        }

        // Generate all 1bp variants
        let mut found = false;
        for pos in 0..32 {
            let original = seq1[pos];
            for &base in &bases {
                if base == original {
                    continue;
                }
                let mut variant = *seq1;
                variant[pos] = base;
                if r2.sequence_set.contains(&variant) {
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }

        if found {
            count += 1;
        }
    }

    count
}

// ─────────────────────────────────────────────────────────────────────────────
// 5kb local window comparison (SynTracker-inspired)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute synteny score by comparing 5kb non-overlapping windows.
///
/// For each window, compute the Jaccard similarity of tag sequences present
/// in both genomes. The final score is the mean Jaccard across all windows
/// that contain at least one tag in either genome. Windows with no tags in
/// both genomes are skipped (not counted as 0).
///
/// This approach is robust to rearrangements and indels because it only
/// compares local tag presence within each window, not global tag order.
pub fn windowed_synteny_score(
    record_a: &TgtRecord,
    record_b: &TgtRecord,
    window_size: u64,
) -> f64 {
    if record_a.total_length == 0 || record_b.total_length == 0 || window_size == 0 {
        return 0.0;
    }

    let windows_a = ((record_a.total_length + window_size - 1) / window_size) as usize;
    let windows_b = ((record_b.total_length + window_size - 1) / window_size) as usize;
    let n_windows = windows_a.max(windows_b);

    // Build window -> tag set for genome A
    let mut window_tags_a: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_a.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_a[w].insert(tag.sequence);
    }

    // Build window -> tag set for genome B
    let mut window_tags_b: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_b.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_b[w].insert(tag.sequence);
    }

    // Compute mean Jaccard across windows with at least one tag
    let mut total_jaccard = 0.0;
    let mut valid_windows = 0usize;

    for i in 0..n_windows {
        let set_a = &window_tags_a[i];
        let set_b = &window_tags_b[i];

        if set_a.is_empty() && set_b.is_empty() {
            continue;
        }

        let intersection: HashSet<_> = set_a.intersection(set_b).collect();
        let union: HashSet<_> = set_a.union(set_b).collect();

        if !union.is_empty() {
            total_jaccard += intersection.len() as f64 / union.len() as f64;
            valid_windows += 1;
        }
    }

    if valid_windows == 0 {
        return 0.0;
    }

    total_jaccard / valid_windows as f64
}

/// Compute approximate windowed synteny score with Hamming distance ≤ 1.
///
/// For each window, tags from genome A are matched against tags from genome B
/// using approximate matching (1bp tolerance). This increases sensitivity
/// for divergent strains where SNPs disrupt exact tag matches.
pub fn windowed_synteny_score_approx(
    record_a: &TgtRecord,
    record_b: &TgtRecord,
    window_size: u64,
) -> f64 {
    if record_a.total_length == 0 || record_b.total_length == 0 || window_size == 0 {
        return 0.0;
    }

    let windows_a = ((record_a.total_length + window_size - 1) / window_size) as usize;
    let windows_b = ((record_b.total_length + window_size - 1) / window_size) as usize;
    let n_windows = windows_a.max(windows_b);

    // Build window -> tag set for genome A
    let mut window_tags_a: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_a.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_a[w].insert(tag.sequence);
    }

    // Build window -> tag set for genome B
    let mut window_tags_b: Vec<HashSet<[u8; 32]>> = vec![HashSet::new(); n_windows];
    for tag in &record_b.tags {
        let w = (tag.position / window_size).min((n_windows - 1) as u64) as usize;
        window_tags_b[w].insert(tag.sequence);
    }

    let mut total_score = 0.0;
    let mut valid_windows = 0usize;

    for i in 0..n_windows {
        let set_a = &window_tags_a[i];
        let set_b = &window_tags_b[i];

        if set_a.is_empty() && set_b.is_empty() {
            continue;
        }

        // Approximate matching: count A tags that match B tags within Hamming-1
        let mut matched = 0usize;
        for seq_a in set_a {
            if set_b.contains(seq_a) {
                matched += 1;
                continue;
            }
            // Generate 1bp variants
            let bases = [b'A', b'T', b'C', b'G'];
            let mut found = false;
            for pos in 0..32 {
                let orig = seq_a[pos];
                for &base in &bases {
                    if base == orig {
                        continue;
                    }
                    let mut variant = *seq_a;
                    variant[pos] = base;
                    if set_b.contains(&variant) {
                        found = true;
                        break;
                    }
                }
                if found {
                    break;
                }
            }
            if found {
                matched += 1;
            }
        }

        let union_size = set_a.len() + set_b.len();
        if union_size > 0 {
            // Score = 2 * matched / union_size (similar to shared_ratio)
            total_score += 2.0 * matched as f64 / union_size as f64;
            valid_windows += 1;
        }
    }

    if valid_windows == 0 {
        return 0.0;
    }

    total_score / valid_windows as f64
}
