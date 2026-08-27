#!/usr/bin/env python3
"""Pairwise structural-synteny wrapper for Syn2b.

Digests two FASTA genomes with a Type IIB enzyme and returns the Syn2b
structural-synteny metrics (canonical-tag, ordered-adjacency, substitution-
invariant) plus legacy adjacency Jaccard.

Usage:
    python3 pairwise_structural_synteny.py genome_A.fasta genome_B.fasta \
        [--enzyme BcgI] [--tmpdir /tmp/syn2b_pair]

Output (TSV to stdout):
    query   reference   enzyme   structural   breakpoints   breakpoint_density
    shared_tags   repeats_dropped   legacy_adjacency
"""

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{r.stderr}")
    return r


def main():
    p = argparse.ArgumentParser(description="Pairwise Syn2b structural synteny")
    p.add_argument("query", help="query FASTA")
    p.add_argument("reference", help="reference FASTA")
    p.add_argument("--enzyme", default="BcgI")
    p.add_argument("--syn2b", default="syn2b",
                   help="path to syn2b binary (default: $PATH)")
    p.add_argument("--tmpdir", default=None,
                   help="working directory (default: tempfile)")
    args = p.parse_args()

    tmp = Path(args.tmpdir) if args.tmpdir else Path(tempfile.mkdtemp(prefix="syn2b_pair_"))
    tmp.mkdir(parents=True, exist_ok=True)

    q_tgt = tmp / "query.tgt"
    r_tgt = tmp / "reference.tgt"

    run([args.syn2b, "digest", "-i", args.query, "-o", str(q_tgt),
         "-e", args.enzyme, "-f", "text"])
    run([args.syn2b, "digest", "-i", args.reference, "-o", str(r_tgt),
         "-e", args.enzyme, "-f", "text"])

    # Synteny expects a directory of .tgt files
    csv_out = tmp / "synteny.csv"
    run([args.syn2b, "synteny", "-i", str(tmp), "-o", str(csv_out)])

    with open(csv_out) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["genome_A"] in (Path(args.query).stem, Path(args.reference).stem):
                print("query\treference\tenzyme\tstructural\tbreakpoints\t"
                      "breakpoint_density\tshared_tags\trepeats_dropped\tlegacy_adjacency")
                print("\t".join([
                    Path(args.query).stem,
                    Path(args.reference).stem,
                    args.enzyme,
                    row["structural"],
                    row["breakpoints"],
                    row["breakpoint_density"],
                    row["shared_tags"],
                    row["repeats_dropped"],
                    row["legacy_adjacency"],
                ]))
                return

    raise RuntimeError("no matching row in synteny output")


if __name__ == "__main__":
    main()
