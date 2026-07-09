#!/usr/bin/env python3
"""Split a single-contig TGT file into multiple contigs to simulate a draft genome."""

import sys
import re
import math


def split_tgt_into_contigs(input_path, output_path, n_contigs=10):
    """Split a TGT file into n_contigs roughly equal-sized contigs."""
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
        # Parse format: Enzyme:SEQ@POS:contig_name or Enzyme:SEQ@POS
        parts = line.split()
        for i in range(0, len(parts), 2):
            tag_part = parts[i]
            gap = int(parts[i+1].strip('-')) if i+1 < len(parts) else 0

            # Parse tag: Enzyme:SEQ@POS[:contig_name]
            enzyme_seq_pos = tag_part.split(':')
            enzyme = enzyme_seq_pos[0]
            seq_pos = enzyme_seq_pos[1].split('@')
            seq = seq_pos[0]
            pos_str = seq_pos[1]
            # contig_name may be present after another ':'
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

    # Write multi-contig TGT
    with open(output_path, 'w') as f:
        # Header: keep original genome_id but add contig_names info in comment
        f.write(f">{genome_id}|length={total_length}\n")

        for i, contig_tags in enumerate(contigs):
            contig_name = f"contig_{i+1}"
            if not contig_tags:
                continue

            # Write contig tags
            for j, (enzyme, seq, pos, gap) in enumerate(contig_tags):
                # Position is global (relative to whole genome), contig_name is suffix
                f.write(f"{enzyme}:{seq}@{pos}:{contig_name}")
                if j < len(contig_tags) - 1:
                    f.write(f" -{gap}- ")
            f.write("\n")

    # Report
    total_tags = sum(len(c) for c in contigs)
    print(f"Split {total_tags} tags into {n_contigs} contigs")
    for i, c in enumerate(contigs):
        if c:
            print(f"  contig_{i+1}: {len(c)} tags, pos range {c[0][2]}-{c[-1][2]}")
        else:
            print(f"  contig_{i+1}: 0 tags")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: split_tgt.py <input.tgt> <output.tgt> <n_contigs>")
        sys.exit(1)

    split_tgt_into_contigs(sys.argv[1], sys.argv[2], int(sys.argv[3]))
