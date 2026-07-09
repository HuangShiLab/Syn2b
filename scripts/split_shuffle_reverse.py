#!/usr/bin/env python3
"""Split a single-contig TGT file into multiple contigs, shuffle, and reverse some."""

import sys
import re
import random


def split_shuffle_reverse_tgt(input_path, output_path, n_contigs=10, seed=42, reverse_indices=None):
    """Split TGT into n_contigs, shuffle, and reverse specified contigs."""
    random.seed(seed)
    
    with open(input_path, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip()
    match = re.match(r'>([^|]+)\|length=(\d+)', header)
    if not match:
        print(f"Cannot parse header: {header}")
        sys.exit(1)

    genome_id = match.group(1)
    total_length = int(match.group(2))

    all_tags = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        for i in range(0, len(parts), 2):
            tag_part = parts[i]
            gap = int(parts[i+1].strip('-')) if i+1 < len(parts) else 0
            enzyme_seq_pos = tag_part.split(':')
            enzyme = enzyme_seq_pos[0]
            seq_pos = enzyme_seq_pos[1].split('@')
            seq = seq_pos[0]
            pos_str = seq_pos[1]
            if ':' in pos_str:
                pos_str = pos_str.split(':')[0]
            pos = int(pos_str)
            all_tags.append((enzyme, seq, pos, gap))

    all_tags.sort(key=lambda x: x[2])

    contig_size = total_length // n_contigs
    contigs = [[] for _ in range(n_contigs)]

    for enzyme, seq, pos, gap in all_tags:
        contig_idx = min(pos // contig_size, n_contigs - 1)
        contigs[contig_idx].append((enzyme, seq, pos, gap))

    contigs = [(i, c) for i, c in enumerate(contigs) if c]
    original_order = [i+1 for i, _ in contigs]
    
    random.shuffle(contigs)
    shuffled_order = [i+1 for i, _ in contigs]
    
    # Reverse specified contigs (reverse tag order within contig)
    if reverse_indices is None:
        reverse_indices = [1, 3, 5, 7]  # reverse half of contigs
    
    reversed_contigs = []
    for idx, (orig_idx, tags) in enumerate(contigs):
        if idx in reverse_indices:
            # Reverse tag order (simulates reverse complement)
            tags = list(reversed(tags))
        reversed_contigs.append((orig_idx, tags))
    
    print(f"Original order: {original_order}")
    print(f"Shuffled order: {shuffled_order}")
    print(f"Reversed contig positions: {reverse_indices}")

    with open(output_path, 'w') as f:
        f.write(f">{genome_id}|length={total_length}\n")

        for idx, (orig_idx, tags) in enumerate(reversed_contigs):
            contig_name = f"contig_{orig_idx+1}"
            
            for j, (enzyme, seq, pos, gap) in enumerate(tags):
                f.write(f"{enzyme}:{seq}@{pos}:{contig_name}")
                if j < len(tags) - 1:
                    f.write(f" -{gap}- ")
            f.write("\n")

    total_tags = sum(len(c) for _, c in reversed_contigs)
    print(f"Wrote {total_tags} tags into {len(reversed_contigs)} contigs")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: split_shuffle_reverse.py <input.tgt> <output.tgt> <n_contigs> [seed]")
        sys.exit(1)

    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42
    split_shuffle_reverse_tgt(sys.argv[1], sys.argv[2], int(sys.argv[3]), seed)
