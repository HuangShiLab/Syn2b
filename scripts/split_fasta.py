#!/usr/bin/env python3
"""Create a synthetic draft genome by splitting a FASTA into contigs and optionally reversing some."""

import sys
import random


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def split_fasta_into_contigs(input_path, output_path, n_contigs=10, reverse_indices=None, seed=42):
    random.seed(seed)
    
    with open(input_path, 'r') as f:
        lines = f.readlines()
    
    # Read sequence
    header = lines[0].strip()
    seq = ''.join(line.strip() for line in lines[1:])
    
    # Split into contigs
    contig_size = len(seq) // n_contigs
    contigs = []
    for i in range(n_contigs):
        start = i * contig_size
        end = start + contig_size if i < n_contigs - 1 else len(seq)
        contig_seq = seq[start:end]
        contigs.append((f"contig_{i+1}", contig_seq))
    
    # Shuffle order
    original_order = list(range(n_contigs))
    random.shuffle(contigs)
    shuffled_order = [int(c[0].split('_')[1]) for c in contigs]
    
    # Reverse specified contigs
    if reverse_indices is None:
        reverse_indices = [1, 3, 5, 7]
    
    print(f"Original order: {original_order}")
    print(f"Shuffled order: {shuffled_order}")
    print(f"Reversing contigs at positions: {reverse_indices}")
    
    for idx in reverse_indices:
        if idx < len(contigs):
            name, seq = contigs[idx]
            contigs[idx] = (name + "_rev", reverse_complement(seq))
            print(f"  Reversed {name} -> {name}_rev")
    
    # Write output FASTA
    with open(output_path, 'w') as f:
        for name, seq in contigs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    
    total_len = sum(len(s) for _, s in contigs)
    print(f"Wrote {len(contigs)} contigs, total {total_len} bp to {output_path}")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: split_fasta.py <input.fasta> <output.fasta> <n_contigs> [seed]")
        sys.exit(1)
    
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42
    split_fasta_into_contigs(sys.argv[1], sys.argv[2], int(sys.argv[3]), seed=seed)
