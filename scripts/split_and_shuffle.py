#!/usr/bin/env python3
"""Split a single-contig TGT file into multiple contigs and shuffle their order."""

import sys
import re
import math
import random


def split_and_shuffle_tgt(input_path, output_path, n_contigs=10, seed=42):
    """Split a TGT file into n_contigs and shuffle their order."""
    random.seed(seed)
    
    with open(input_path, 'r') as f:
        lines = f.readlines()

    # Parse header
    header = lines[0].strip()
    match = re.match(r'>([^|]+)\|length=(\d+)', header)
    if not match:
        print(f"Cannot parse header: {header}")
        sys.exit(1)

    genome_id = match.group(1)
    total_length = int(match.group(2))

    # Parse all tags
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

    # Sort by position
    all_tags.sort(key=lambda x: x[2])

    # Divide into contigs
    contig_size = total_length // n_contigs
    contigs = [[] for _ in range(n_contigs)]

    for enzyme, seq, pos, gap in all_tags:
        contig_idx = min(pos // contig_size, n_contigs - 1)
        contigs[contig_idx].append((enzyme, seq, pos, gap))

    # Filter out empty contigs
    contigs = [(i, c) for i, c in enumerate(contigs) if c]
    
    # Remember original order
    original_order = [i+1 for i, _ in contigs]
    
    # Shuffle contigs (but keep internal tag order)
    random.shuffle(contigs)
    shuffled_order = [i+1 for i, _ in contigs]
    
    print(f"Original contig order: {original_order}")
    print(f"Shuffled contig order: {shuffled_order}")

    # Write shuffled multi-contig TGT
    with open(output_path, 'w') as f:
        f.write(f">{genome_id}|length={total_length}\n")

        for idx, (orig_idx, contig_tags) in enumerate(contigs):
            contig_name = f"contig_{orig_idx+1}"
            
            for j, (enzyme, seq, pos, gap) in enumerate(contig_tags):
                f.write(f"{enzyme}:{seq}@{pos}:{contig_name}")
                if j < len(contig_tags) - 1:
                    f.write(f" -{gap}- ")
            f.write("\n")

    # Report
    total_tags = sum(len(c) for _, c in contigs)
    print(f"Split {total_tags} tags into {len(contigs)} contigs (shuffled)")
    for idx, (orig_idx, c) in enumerate(contigs):
        print(f"  position {idx}: contig_{orig_idx+1} ({len(c)} tags)")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: split_and_shuffle.py <input.tgt> <output.tgt> <n_contigs> [seed]")
        sys.exit(1)

    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42
    split_and_shuffle_tgt(sys.argv[1], sys.argv[2], int(sys.argv[3]), seed)
