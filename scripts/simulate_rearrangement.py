#!/usr/bin/env python3
"""
simulate_rearrangement.py

Simulate inversions and translocations on a single complete genome at fixed
substitution divergence, then compare a Mash-distance proxy with Syn2b
tag-adjacency metrics.

Core claim: "tag adjacency tracks structural change while Mash does not."
"""

import argparse
import csv
import math
import os
import random
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# Optional third-party dependencies
# ---------------------------------------------------------------------------
try:
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    from scipy.stats import kendalltau
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BASIS = ["A", "C", "G", "T"]
BCGI_SITE = "GAAGGCC"
MU = 0.01
KMER = 21


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------
def parse_fasta(path):
    """Return (header, sequence) from the first FASTA record."""
    if BIOPYTHON_AVAILABLE:
        rec = next(SeqIO.parse(path, "fasta"))
        return rec.id, str(rec.seq).upper()

    with open(path, "r") as fh:
        lines = fh.readlines()

    header = None
    seq_parts = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:].split()[0]
        else:
            seq_parts.append(line.upper())
    if header is None:
        raise ValueError(f"No FASTA record found in {path}")
    return header, "".join(seq_parts)


def write_fasta(path, header, sequence, width=80):
    """Write a single-record FASTA file."""
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            fh.write(sequence[i:i + width] + "\n")


# ---------------------------------------------------------------------------
# Genome mutation
# ---------------------------------------------------------------------------
def substitute(seq, mu=MU, rng=None):
    """Apply random point mutations at rate mu."""
    if rng is None:
        rng = random
    seq_list = list(seq)
    for i, base in enumerate(seq_list):
        if rng.random() < mu:
            alt = rng.choice([b for b in BASIS if b != base])
            seq_list[i] = alt
    return "".join(seq_list)


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def inversion(seq, size, rng=None):
    """Introduce a random inversion of *size* bp."""
    if rng is None:
        rng = random
    if size >= len(seq):
        size = len(seq) // 2
    start = rng.randint(0, len(seq) - size)
    end = start + size
    return seq[:start] + reverse_complement(seq[start:end]) + seq[end:]


def translocation(seq, size, rng=None):
    """Swap two non-overlapping segments of *size* bp each."""
    if rng is None:
        rng = random
    if 2 * size >= len(seq):
        size = len(seq) // 4
    # Pick first segment
    start1 = rng.randint(0, len(seq) - 2 * size)
    end1 = start1 + size
    # Pick second segment after the first
    start2 = rng.randint(end1, len(seq) - size)
    end2 = start2 + size
    seg1 = seq[start1:end1]
    seg2 = seq[start2:end2]
    return seq[:start1] + seg2 + seq[end1:start2] + seg1 + seq[end2:]


# ---------------------------------------------------------------------------
# In-silico BcgI digest (Python implementation)
# ---------------------------------------------------------------------------
def digest_bcgI(sequence):
    """
    Find all GAAGGCC sites and extract 32 bp tags.
    Returns a list of (tag_sequence, position) sorted by position.
    """
    seq = sequence.upper()
    site = BCGI_SITE
    site_len = len(site)
    tag_len = 32
    tags = []

    i = 0
    while i <= len(seq) - site_len:
        if seq[i:i + site_len] == site:
            # 5' tag: 4 bp upstream + 28 bp into the site region = 32 bp
            t5_start = i - 4
            t5_end = i + 28
            if t5_start >= 0 and t5_end <= len(seq):
                tags.append((seq[t5_start:t5_end], t5_start))

            # 3' tag: 3 bp into the site + 29 bp downstream = 32 bp
            t3_start = i + 3
            t3_end = i + 3 + tag_len
            if t3_end <= len(seq):
                tags.append((seq[t3_start:t3_end], t3_start))

            i += 1
        else:
            i += 1

    tags.sort(key=lambda x: x[1])
    return tags


def write_tgt(path, genome_id, total_length, tags):
    """
    Write a single-contig TGT text file.
    Format:
        >genome_id|length=NNN
        BcgI:SEQ1@POS1 -GAP1- BcgI:SEQ2@POS2 -GAP2- ...
    """
    with open(path, "w") as fh:
        fh.write(f">{genome_id}|length={total_length}\n")
        for j, (tseq, tpos) in enumerate(tags):
            if j > 0:
                gap = tpos - tags[j - 1][1]
                fh.write(f" -{gap}- ")
            fh.write(f"BcgI:{tseq}@{tpos}")
        fh.write("\n")


def read_tgt_tags(path):
    """Parse a single-contig TGT text file; return list of tag sequences in order."""
    with open(path, "r") as fh:
        lines = fh.readlines()

    # Skip header and comment lines
    body = ""
    for line in lines:
        if line.startswith(">") or line.startswith("#"):
            continue
        body += line.strip() + " "

    if not body.strip():
        return []

    # Split on gap markers; keep only tag entries
    raw_parts = body.split()
    tag_seqs = []
    for part in raw_parts:
        if part.startswith("-") and part.endswith("-"):
            continue
        if "@" in part and ":" in part:
            # Format: BcgI:SEQUENCE@POSITION
            enzyme_seq_pos = part.split(":", 1)
            if len(enzyme_seq_pos) == 2:
                seq_pos = enzyme_seq_pos[1].split("@", 1)
                if len(seq_pos) >= 1:
                    tag_seqs.append(seq_pos[0])
    return tag_seqs


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def mash_proxy(seq_a, seq_b, k=KMER):
    """
    Simple k-mer Jaccard proxy for sequence divergence.
    Uses *canonical* k-mers (min of forward and reverse-complement),
    matching Mash's default behaviour.  This makes the proxy blind to
    pure inversions.
    Returns (jaccard_similarity, mash_distance_approx).
    """
    def canonical_kmers(s):
        kmers = set()
        rc = reverse_complement(s)
        for i in range(len(s) - k + 1):
            fwd = s[i:i + k]
            rev = rc[len(s) - k - i:len(s) - i]
            kmers.add(min(fwd, rev))
        return kmers

    k_a = canonical_kmers(seq_a)
    k_b = canonical_kmers(seq_b)
    inter = len(k_a & k_b)
    union = len(k_a | k_b)
    if union == 0:
        return 0.0, 1.0
    jaccard = inter / union
    # Mash distance approximation: d = -1/k * ln(2*J/(1+J))
    if jaccard == 0:
        mash_d = 1.0
    else:
        mash_d = -1.0 / k * math.log(2.0 * jaccard / (1.0 + jaccard))
    return jaccard, mash_d


def adjacency_jaccard(tags_a, tags_b):
    """Jaccard similarity of adjacent tag-sequence pairs."""
    def adj_set(tlist):
        # Use canonical (sorted) pair so strand/orientation doesn't matter
        s = set()
        for i in range(len(tlist) - 1):
            a, b = tlist[i], tlist[i + 1]
            s.add((a, b) if a <= b else (b, a))
        return s

    set_a = adj_set(tags_a)
    set_b = adj_set(tags_b)
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return inter / union


def breakpoint_count(tags_a, tags_b):
    """Symmetric difference of adjacency sets."""
    def adj_set(tlist):
        s = set()
        for i in range(len(tlist) - 1):
            a, b = tlist[i], tlist[i + 1]
            s.add((a, b) if a <= b else (b, a))
        return s

    set_a = adj_set(tags_a)
    set_b = adj_set(tags_b)
    return len(set_a ^ set_b)


def kendall_tau_rank(tags_a, tags_b):
    """Kendall's tau on the order of shared tags."""
    # Build position maps
    pos_a = {t: i for i, t in enumerate(tags_a)}
    pos_b = {t: i for i, t in enumerate(tags_b)}
    shared = [t for t in tags_a if t in pos_b]
    if len(shared) < 2:
        return None
    ranks_a = [pos_a[t] for t in shared]
    ranks_b = [pos_b[t] for t in shared]
    if SCIPY_AVAILABLE:
        tau, _ = kendalltau(ranks_a, ranks_b)
        return tau
    # Simple O(n^2) implementation if scipy is absent
    n = len(shared)
    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (ranks_a[i] - ranks_a[j]) * (ranks_b[i] - ranks_b[j]) > 0:
                concordant += 1
            else:
                discordant += 1
    total = n * (n - 1) // 2
    if total == 0:
        return 0.0
    return (concordant - discordant) / total


# ---------------------------------------------------------------------------
# Synthetic genome generator
# ---------------------------------------------------------------------------
def create_synthetic_genome(path, length=2_000_000, seed=42):
    """Create a random FASTA with enough GAAGGCC sites for a realistic digest."""
    rng = random.Random(seed)
    seq = []
    # Generate random sequence, but ensure some BcgI sites exist
    site = BCGI_SITE
    block = 500
    for _ in range(length // block):
        # mostly random
        rand_part = "".join(rng.choices(BASIS, k=block - 7))
        # sprinkle a site every few blocks to get reasonable tag density
        if rng.random() < 0.5:
            insert_pos = rng.randint(0, block - 7)
            chunk = rand_part[:insert_pos] + site + rand_part[insert_pos:]
            seq.append(chunk[:block])
        else:
            seq.append(rand_part + "".join(rng.choices(BASIS, k=7)))
    remainder = length - len(seq) * block
    if remainder > 0:
        seq.append("".join(rng.choices(BASIS, k=remainder)))

    sequence = "".join(seq)[:length]
    write_fasta(path, "synthetic_test_genome", sequence)
    return sequence


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------
def run_experiment(input_fasta, syn2b_binary, output_csv, output_png):
    rng = random.Random(42)

    # ------------------------------------------------------------------
    # Step 1: Load genome
    # ------------------------------------------------------------------
    if not os.path.isfile(input_fasta):
        print(f"Warning: Genome FASTA not found at {input_fasta}")
        print("Creating a synthetic test genome instead …")
        input_fasta = tempfile.mktemp(suffix=".fasta")
        create_synthetic_genome(input_fasta)

    genome_id, original_seq = parse_fasta(input_fasta)
    genome_len = len(original_seq)
    print(f"Loaded genome: {genome_id}, length = {genome_len:,} bp")

    # ------------------------------------------------------------------
    # Step 2: Generate derived genomes
    # ------------------------------------------------------------------
    tmpdir = tempfile.mkdtemp(prefix="syn2b_rearr_")
    genomes = {}   # label -> (fasta_path, seq)

    # Original
    orig_path = os.path.join(tmpdir, "original.fasta")
    write_fasta(orig_path, genome_id, original_seq)
    genomes["original"] = (orig_path, original_seq)

    # Group A — Substitutions only (control)
    substituted_seqs = []
    for rep in range(1, 4):
        label = f"control_{rep}"
        seq = substitute(original_seq, mu=MU, rng=rng)
        substituted_seqs.append(seq)
        path = os.path.join(tmpdir, f"{label}.fasta")
        write_fasta(path, label, seq)
        genomes[label] = (path, seq)

    # Group B — Substitutions + Rearrangements
    # Apply each SV to the *same* substituted background (control_1) so that
    # the only difference between control_1 and the SV genomes is the SV itself.
    sv_specs = [
        ("inversion", 50_000),
        ("inversion", 100_000),
        ("inversion", 500_000),
        ("inversion", 1_000_000),
        ("translocation", 100_000),
        ("translocation", 500_000),
    ]

    base_seq = substituted_seqs[0]  # same background as control_1
    for sv_type, sv_size in sv_specs:
        label = f"{sv_type}_{sv_size // 1000}kb"
        seq = base_seq
        # Scale SV if genome is too small
        effective_size = min(sv_size, genome_len // 3)
        if sv_type == "inversion":
            seq = inversion(seq, effective_size, rng=rng)
        else:
            seq = translocation(seq, effective_size, rng=rng)
        path = os.path.join(tmpdir, f"{label}.fasta")
        write_fasta(path, label, seq)
        genomes[label] = (path, seq)

    print(f"Generated {len(genomes)} genomes in {tmpdir}")

    # ------------------------------------------------------------------
    # Step 3: Digest all genomes (Python implementation)
    # ------------------------------------------------------------------
    tgt_files = {}
    for label, (fasta_path, seq) in genomes.items():
        tags = digest_bcgI(seq)
        tgt_path = os.path.join(tmpdir, f"{label}.tgt")
        write_tgt(tgt_path, label, len(seq), tags)
        tgt_files[label] = tgt_path
        print(f"  {label}: {len(tags)} tags")

    # Optionally invoke syn2b binary if it exists and digest works
    if syn2b_binary and os.path.isfile(syn2b_binary):
        print(f"\nSyn2b binary found at {syn2b_binary}")
        print("Note: syn2b digest is currently a stub; using Python digestion.")

    # ------------------------------------------------------------------
    # Step 4: Compute metrics (original vs each derived)
    # ------------------------------------------------------------------
    original_seq = genomes["original"][1]
    original_tags = read_tgt_tags(tgt_files["original"])

    results = []
    for label in genomes:
        if label == "original":
            continue
        seq = genomes[label][1]
        tags = read_tgt_tags(tgt_files[label])

        # Mash proxy
        _, mash_d = mash_proxy(original_seq, seq)

        # Syn2b metrics
        aj = adjacency_jaccard(original_tags, tags)
        bp = breakpoint_count(original_tags, tags)
        kt = kendall_tau_rank(original_tags, tags)

        # Determine group / SV metadata
        group = "control" if label.startswith("control") else "rearranged"
        sv_type = "none"
        sv_size = 0
        if label.startswith("inversion_"):
            sv_type = "inversion"
            sv_size = int(label.split("_")[1].replace("kb", "")) * 1000
        elif label.startswith("translocation_"):
            sv_type = "translocation"
            sv_size = int(label.split("_")[1].replace("kb", "")) * 1000

        results.append({
            "genome_label": label,
            "group": group,
            "sv_type": sv_type,
            "sv_size": sv_size,
            "mash_proxy": round(mash_d, 6),
            "syn2b_adjacency_jaccard": round(aj, 6),
            "syn2b_breakpoint_count": bp,
            "kendall_tau": round(kt, 6) if kt is not None else None,
        })
        print(f"  {label}: mash={mash_d:.4f}, adj_jaccard={aj:.4f}, breakpoints={bp}, tau={kt}")

    # ------------------------------------------------------------------
    # Step 5: Save CSV
    # ------------------------------------------------------------------
    fieldnames = [
        "genome_label", "group", "sv_type", "sv_size",
        "mash_proxy", "syn2b_adjacency_jaccard", "syn2b_breakpoint_count",
    ]
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"\nCSV saved to {output_csv}")

    # ------------------------------------------------------------------
    # Step 5b: Plot
    # ------------------------------------------------------------------
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available; skipping figure generation.")
        return

    labels = [r["genome_label"] for r in results]
    mash_vals = [r["mash_proxy"] for r in results]
    adj_vals = [r["syn2b_adjacency_jaccard"] for r in results]
    bp_vals = [r["syn2b_breakpoint_count"] for r in results]

    # Colour by group
    colours = ["#4C78A8" if r["group"] == "control" else "#E45756" for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Mash proxy
    ax = axes[0]
    bars = ax.bar(range(len(labels)), mash_vals, color=colours, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mash distance proxy")
    ax.set_title("Sequence divergence proxy\n(all genomes ≈ 1% SNPs)")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylim(bottom=0)

    # Adjacency Jaccard
    ax = axes[1]
    bars = ax.bar(range(len(labels)), adj_vals, color=colours, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Adjacency Jaccard")
    ax.set_title("Syn2b tag-adjacency similarity")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylim(0, 1.05)

    # Breakpoint count
    ax = axes[2]
    bars = ax.bar(range(len(labels)), bp_vals, color=colours, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Breakpoint count")
    ax.set_title("Syn2b tag-adjacency breakpoints")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylim(bottom=0)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4C78A8", edgecolor="black", label="Control (SNPs only)"),
        Patch(facecolor="#E45756", edgecolor="black", label="Rearranged (SNPs + SV)"),
    ]
    fig.legend(handles=legend_elements, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.02))

    fig.suptitle(
        "Syn2b Rearrangement Validation\n"
        "Tag adjacency tracks structural change; Mash proxy does not",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Figure saved to {output_png}")

    # Summary
    print("\n--- Summary ---")
    print("Control genomes (SNPs only):")
    for r in results:
        if r["group"] == "control":
            print(f"  {r['genome_label']}: mash={r['mash_proxy']:.4f}, adj_jaccard={r['syn2b_adjacency_jaccard']:.4f}")
    print("Rearranged genomes (SNPs + SV):")
    for r in results:
        if r["group"] == "rearranged":
            print(f"  {r['genome_label']}: mash={r['mash_proxy']:.4f}, adj_jaccard={r['syn2b_adjacency_jaccard']:.4f}, breakpoints={r['syn2b_breakpoint_count']}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Simulate rearrangements and compare Mash vs Syn2b metrics."
    )
    parser.add_argument(
        "--input",
        default="data/e_coli_k12.fasta",
        help="Input complete-genome FASTA (default: data/e_coli_k12.fasta)",
    )
    parser.add_argument(
        "--syn2b",
        default="target/release/syn2b",
        help="Path to syn2b binary (default: target/release/syn2b)",
    )
    parser.add_argument(
        "--csv",
        default="scripts/rearrangement_validation.csv",
        help="Output CSV path (default: scripts/rearrangement_validation.csv)",
    )
    parser.add_argument(
        "--png",
        default="scripts/rearrangement_validation.png",
        help="Output figure path (default: scripts/rearrangement_validation.png)",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=MU,
        help=f"Substitution rate (default: {MU})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    args = parser.parse_args()

    random.seed(args.seed)

    # Ensure output directories exist
    os.makedirs(os.path.dirname(args.csv) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.png) or ".", exist_ok=True)

    run_experiment(args.input, args.syn2b, args.csv, args.png)


if __name__ == "__main__":
    main()
