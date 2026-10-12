#!/usr/bin/env python3
"""Convert antiSMASH 8 JSON outputs into Syn2b feature tables."""
from __future__ import annotations
import argparse,csv,glob,json,re
from pathlib import Path
FIELDS=["gene_id","genome_id","contig_id","start","end","strand","product","ko","pfam","cog","eggnog","smash_tool","cluster_id"]
def parse_location(text):
 # BioJSON location format: [start:end](+/-)
 match=re.fullmatch(r"\[(\d+):(\d+)\](?:\(([+-])\))?",text.strip())
 if not match:raise ValueError(f"unsupported location {text}")
 return int(match.group(1))+1,int(match.group(2)),match.group(3) or "+"
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--antismash-dir',required=True);ap.add_argument('--output',required=True)
 args=ap.parse_args();features=[];regions_seen=0
 for path in Path(args.antismash_dir).glob('*/*.json'):
  data=json.load(open(path));genome=Path(data.get('input_file','genome')).stem
  for record in data.get('records',[]):
   regions=[];cds=[]
   for feature in record.get('features',[]):
    kind=feature.get('type');qualifiers=feature.get('qualifiers',{})
    start,end,strand=parse_location(feature['location'])
    if kind=='region':
     regions_seen+=1;number=(qualifiers.get('region_number') or [''])[0]
     product=';'.join(qualifiers.get('product',[])) or 'BGC'
     cluster=f"ANTISMASH:{record['id']}:region{number}";regions.append((start,end,cluster,product,qualifiers))
    elif kind=='CDS':
     locus=(qualifiers.get('locus_tag') or [f"{record['id']}_{start}_{end}"])[0]
     gene_functions=';'.join(qualifiers.get('gene_functions',[]))
     pfam=';'.join(sorted(set(re.findall(r'PF\d{5}',gene_functions))))
     cds.append((f"{record['id']}|{locus}",record['id'],start,end,strand,pfam))
   for gene_id,contig,start,end,strand,pfam in cds:
    hit=next(((rs,re_,cluster,product,qual) for rs,re_,cluster,product,qual in regions if start>=rs and end<=re_),None)
    if hit is None:continue
    region_start,region_end,cluster,product,qualifiers=hit
    category=';'.join(qualifiers.get('product',[]))
    features.append({"gene_id":gene_id,"genome_id":genome,"contig_id":contig,"start":start,"end":end,
      "strand":strand,"product":f"antiSMASH_gene;{product}","ko":"","pfam":pfam,"cog":"","eggnog":"",
      "smash_tool":f"ANTISMASH:{category}","cluster_id":cluster})
 with open(args.output,'w',newline='') as fh:
  writer=csv.DictWriter(fh,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');writer.writeheader()
  for feature in sorted(features,key=lambda r:(r['cluster_id'],r['contig_id'],r['start'],r['end'],r['gene_id'])):writer.writerow(feature)
 print(f"wrote {len(features)} CDS features from {regions_seen} antiSMASH regions")
if __name__=='__main__':main()
