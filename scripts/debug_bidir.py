#!/usr/bin/env python3
"""Debug script to analyze scaffold orientation results for new TGT files."""

import re

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def parse_tgt(path):
    with open(path) as f:
        lines = f.readlines()
    all_tags = {}
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        for i in range(0, len(parts), 2):
            tag_part = parts[i]
            enzyme, rest = tag_part.split(':', 1)
            seq_pos = rest.split('@')
            seq = seq_pos[0]
            pos_name = seq_pos[1]
            if ':' in pos_name:
                pos, name = pos_name.split(':', 1)
            else:
                pos = pos_name
                name = "unknown"
            pos = int(pos)
            if name not in all_tags:
                all_tags[name] = []
            all_tags[name].append((seq, pos))
    return all_tags

def analyze(ref_tags, draft_tags):
    ref_map = {}
    for contig, tags in ref_tags.items():
        for seq, pos in tags:
            ref_map[seq] = pos
    
    for contig, tags in sorted(draft_tags.items()):
        fwd = sum(1 for seq, _ in tags if seq in ref_map)
        rev = sum(1 for seq, _ in tags if reverse_complement(seq) in ref_map)
        print(f"{contig}: fwd={fwd} rev={rev} ratio={rev/fwd if fwd > 0 else 'inf'}")

ref = parse_tgt("/tmp/syn2b_benchmark/tgts/E_coli_K-12_MG1655.tgt")
draft = parse_tgt("/tmp/syn2b_benchmark/K12_draft_10contigs_rev.tgt")

print("Reference tags:", sum(len(t) for t in ref.values()))
print("Draft tags:", sum(len(t) for t in draft.values()))
print()
analyze(ref, draft)
