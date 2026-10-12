#!/usr/bin/env python3
"""Select highest-quality OAPGC genome per SGB cluster."""
from __future__ import annotations
import argparse,csv,os
from collections import defaultdict
from pathlib import Path
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--metadata',required=True);ap.add_argument('--tgt-dir',required=True);ap.add_argument('--link-dir',required=True);ap.add_argument('--manifest',required=True)
 args=ap.parse_args();link=Path(args.link_dir);link.mkdir(parents=True,exist_ok=True);groups=defaultdict(list)
 with open(args.metadata) as fh:
  for row in csv.DictReader(fh,delimiter='\t'):
   path=Path(args.tgt_dir)/f"{row['genome_id']}.tgt"
   try:comp=float(row.get('completeness') or 0);cont=float(row.get('contamination') or 0);n50=int(row.get('n50') or 0)
   except ValueError:continue
   if path.exists():groups[row['cluster']].append((row,path,comp,cont,n50))
 selected=[]
 for cluster,rows in sorted(groups.items()):
  rows.sort(key=lambda x:(-x[2],x[3],-x[4],x[0]['genome_id']));row,path,_,_,_=rows[0]
  dst=link/f"{row['genome_id']}.tgt"
  if not dst.exists():os.symlink(path,dst)
  selected.append({'cluster':cluster,'genome_id':row['genome_id'],'completeness':row['completeness'],'contamination':row['contamination'],'n50':row['n50'],'path':str(path)})
 with open(args.manifest,'w',newline='') as fh:
  writer=csv.DictWriter(fh,fieldnames=list(selected[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(selected)
 print(f"selected {len(selected)} representatives")
if __name__=='__main__':main()
