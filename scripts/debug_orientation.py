#!/usr/bin/env python3
"""Debug script to analyze scaffold orientation results."""

import re

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

def parse_tgt(path):
    with open(path) as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    genome_id = header[1:].split('|')[0]
    
    # Parse all tags
    all_tags = {}
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        contig_name = None
        for i in range(0, len(parts), 2):
            tag_part = parts[i]
            # Parse: Enzyme:SEQ@POS:contig_name
            enzyme, rest = tag_part.split(':', 1)
            seq_pos = rest.split('@')
            seq = seq_pos[0]
            pos_name = seq_pos[1]
            if ':' in pos_name:
                pos, name = pos_name.split(':', 1)
                contig_name = name
            else:
                pos = pos_name
            pos = int(pos)
            
            if contig_name not in all_tags:
                all_tags[contig_name] = []
            all_tags[contig_name].append((seq, pos))
    
    return all_tags

def analyze_orientation(ref_tags, draft_tags):
    # Build ref map
    ref_map = {}
    for contig, tags in ref_tags.items():
        for seq, pos in tags:
            ref_map[seq] = pos
    
    # Build rc ref map
    ref_rc_map = {}
    for contig, tags in ref_tags.items():
        for seq, pos in tags:
            rc = reverse_complement(seq)
            ref_rc_map[rc] = pos
    
    for contig, tags in sorted(draft_tags.items()):
        # Forward orientation
        fwd_matches = []
        for seq, pos in tags:
            if seq in ref_map:
                fwd_matches.append((ref_map[seq], pos))
        
        # Reverse orientation
        rev_matches = []
        for seq, pos in tags:
            rc = reverse_complement(seq)
            if rc in ref_map:
                rev_matches.append((ref_map[rc], pos))
        
        # Compute concordance
        def concordance(matches):
            if len(matches) < 2:
                return 0
            matches.sort(key=lambda x: x[1])  # sort by draft pos
            conc = 0
            for i in range(1, len(matches)):
                draft_order = matches[i][1] > matches[i-1][1]
                ref_order = matches[i][0] > matches[i-1][0]
                if draft_order == ref_order:
                    conc += 1
            return conc
        
        fwd_conc = concordance(fwd_matches)
        rev_conc = concordance(rev_matches)
        
        chosen = "REV" if rev_conc > fwd_conc or (rev_conc == fwd_conc and len(rev_matches) >= len(fwd_matches)) else "FWD"
        
        print(f"{contig}: fwd={len(fwd_matches)}(c={fwd_conc}) rev={len(rev_matches)}(c={rev_conc}) -> {chosen}")

ref = parse_tgt("/tmp/syn2b_benchmark/tgts/E_coli_K-12_MG1655.tgt")
draft = parse_tgt("/tmp/syn2b_benchmark/K12_draft_10contigs_rev.tgt")

print("Reference contigs:", list(ref.keys()))
print("Draft contigs:", list(draft.keys()))
print()
analyze_orientation(ref, draft)
