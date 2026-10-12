#!/usr/bin/env python3
import argparse,csv
from collections import defaultdict
if __name__=="__main__":
 ap=argparse.ArgumentParser();ap.add_argument("--antismash",required=True);ap.add_argument("--eggnog",required=True);ap.add_argument("--representatives",required=True);ap.add_argument("--output",required=True)
 a=ap.parse_args()
 rep={}
 with open(a.representatives) as fh:
  for row in csv.DictReader(fh,delimiter="\t"):rep[row["genome_id"]]=row["cluster"]
 by=defaultdict(list)
 with open(a.eggnog) as fh:
  for row in csv.DictReader(fh,delimiter="\t"):
   cluster=rep.get(row["genome_id"],row["genome_id"])
   suffix=row["contig_id"].rsplit("_ctg_",1)[-1]
   by[(cluster,suffix)].append((int(row["start"]),int(row["end"]),row))
 n=sum(len(v) for v in by.values());hits=updated=0
 with open(a.antismash) as src,open(a.output,"w",newline="") as out:
  fields=None;writer=None
  for row in csv.DictReader(src,delimiter="\t"):
   if fields is None:
    fields=list(row);writer=csv.DictWriter(out,fieldnames=fields,delimiter="\t",lineterminator="\n");writer.writeheader()
   cluster=rep.get(row["genome_id"],row["genome_id"]);suffix=row["contig_id"].rsplit("_ctg_",1)[-1]
   start=int(row["start"]);end=int(row["end"]);candidates=[]
   for es,ee,egg in by.get((cluster,suffix),[]):
    if es<=end and ee>=start:candidates.append((min(ee,end)-max(es,start)+1,egg))
   if candidates:
    hits+=1;overlap,egg=max(candidates,key=lambda x:(x[0],x[2]["gene_id"]))
    for field in ["product","ko","pfam","eggnog"]:
     if not row[field] and egg[field]:row[field]=egg[field];updated+=1
   writer.writerow(row)
 print(f"eggnog_intervals={n} antismash_genes={hits} updated_fields={updated}")
