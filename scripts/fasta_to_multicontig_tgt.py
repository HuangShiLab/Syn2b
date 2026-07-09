#!/usr/bin/env python3
"""Generate multi-contig TGT file by digesting a multi-contig FASTA with BcgI."""

import sys
import re


def digest_bcgI(sequence, contig_id, offset):
    """Find all BcgI sites (GAAGGCC) and extract 32bp tags."""
    site = b"GAAGGCC"
    site_len = len(site)
    tag_len = 32
    tags = []
    
    seq = sequence.upper()
    
    i = 0
    while i <= len(seq) - site_len:
        if seq[i:i+site_len] == site:
            # BcgI: 5' overhang, 3bp
            # tag_5p = 4bp before site + 28bp from site end = 32bp
            tag_5p_start = max(0, i - 4)
            tag_5p_end = i + site_len + 24
            if tag_5p_end <= len(seq) and (i - 4) >= 0:
                tag_5p = seq[tag_5p_start:tag_5p_end]
                pos_5p = offset + i - 4
                tags.append((tag_5p.decode(), pos_5p))
            
            # tag_3p = 28bp from site start + 4bp after = 32bp
            tag_3p_start = i + 3
            tag_3p_end = i + 3 + 32
            if tag_3p_end <= len(seq):
                tag_3p = seq[tag_3p_start:tag_3p_end]
                pos_3p = offset + i + 3
                tags.append((tag_3p.decode(), pos_3p))
            
            i += 1  # move forward 1 (sites can overlap)
        else:
            i += 1
    
    return tags


def generate_multi_contig_tgt(fasta_path, output_path):
    """Read multi-contig FASTA and generate TGT file."""
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    contigs = []
    current_name = None
    current_seq = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_name is not None:
                contigs.append((current_name, ''.join(current_seq)))
            current_name = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_name is not None:
        contigs.append((current_name, ''.join(current_seq)))
    
    total_length = sum(len(seq) for _, seq in contigs)
    genome_id = contigs[0][0].split('.')[0]  # use first contig prefix as genome ID
    
    # Digest all contigs
    all_tags = []
    contig_names = []
    offset = 0
    
    for idx, (name, seq) in enumerate(contigs):
        contig_id = idx + 1  # 1-based
        contig_names.append(name)
        
        tags = digest_bcgI(seq.encode(), contig_id, offset)
        for tag_seq, pos in tags:
            all_tags.append((tag_seq, pos, contig_id, name))
        
        offset += len(seq)
    
    # Sort by contig_id, then position
    all_tags.sort(key=lambda x: (x[2], x[1]))
    
    # Write TGT
    with open(output_path, 'w') as f:
        f.write(f">{genome_id}|length={total_length}\n")
        
        # Group by contig
        current_contig = None
        contig_tags = []
        
        for tag_seq, pos, contig_id, name in all_tags:
            if contig_id != current_contig:
                if contig_tags:
                    # Write previous contig tags
                    for j, (t, p) in enumerate(contig_tags):
                        f.write(f"BcgI:{t}@{p}:{name}")
                        if j < len(contig_tags) - 1:
                            gap = contig_tags[j+1][1] - p
                            f.write(f" -{gap}- ")
                    f.write("\n")
                current_contig = contig_id
                contig_tags = []
            contig_tags.append((tag_seq, pos))
        
        # Write last contig
        if contig_tags:
            for j, (t, p) in enumerate(contig_tags):
                f.write(f"BcgI:{t}@{p}:{name}")
                if j < len(contig_tags) - 1:
                    gap = contig_tags[j+1][1] - p
                    f.write(f" -{gap}- ")
            f.write("\n")
    
    print(f"Genome: {genome_id}, total length: {total_length}")
    print(f"Contigs: {len(contig_names)}")
    print(f"Total tags: {len(all_tags)}")
    for i, (name, seq) in enumerate(contigs[:5]):
        n_tags = sum(1 for t in all_tags if t[3] == name)
        print(f"  {name}: {len(seq)} bp, {n_tags} tags")
    print(f"  ...")
    for i, (name, seq) in enumerate(contigs[-5:]):
        n_tags = sum(1 for t in all_tags if t[3] == name)
        print(f"  {name}: {len(seq)} bp, {n_tags} tags")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: fasta_to_multicontig_tgt.py <input.fasta> <output.tgt>")
        sys.exit(1)
    
    generate_multi_contig_tgt(sys.argv[1], sys.argv[2])
