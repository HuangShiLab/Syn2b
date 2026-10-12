#!/usr/bin/env python3
import argparse,csv,glob,gzip,hashlib,json
from pathlib import Path
from collections import defaultdict
if __name__=="__main__":
 ap=argparse.ArgumentParser();ap.add_argument("--faa-dir",required=True);ap.add_argument("--eggnog-features",required=True);ap.add_argument("--antismash-dir",required=True);ap.add_argument("--output",required=True)
 a=ap.parse_args()
 seqmap={}
 for path in glob.glob(a.faa_dir+"/**/*.faa.gz",recursive=True):
  name=None;parts=[]
  with gzip.open(path,'rt') as fh:
   for line in fh:
    line=line.rstrip()
    if line.startswith(">"):
     if name is not None:seqmap.setdefault(hashlib.sha256("".join(parts).replace("*","").encode()).hexdigest(),name)
     name=line[1:].split()[0];parts=[]
    else:parts.append(line)
 print('seqmap',len(seqmap))
 egg={}
 with open(a.eggnog_features) as eh:
  for row in csv.DictReader(eh,delimiter='\t'):egg[row['gene_id']]=row
 print('egg',len(egg))
 fields=None;rows=[]
 njson=0
 for path in glob.glob(a.antismash_dir+"/*/*.json"):
  njson+=1
  if njson<=2:print('json',path)
  data=json.load(open(path));genome=Path(data.get("input_file","genome")).stem
  regions={};cdss=[]
  for record in data.get("records",[]):
   regions=[];cdss=[]
   for feature in record.get("features",[]):
    loc=feature.get("location","");kind=feature.get("type");qual=feature.get("qualifiers",{})
    import re
    m=re.fullmatch(r"\[(\d+):(\d+)\](?:\(([+-])\))?",loc)
    if not m:continue
    start=int(m.group(1))+1;end=int(m.group(2));strand=m.group(3) or "+"
    if kind=="region":
     cluster=f"ANTISMASH:{record['id']}:{Path(path).stem}:region{(qual.get('region_number') or [''])[0]}"
     product=";".join(qual.get("product",[])) or "BGC";regions.append((start,end,cluster,product,qual))
    elif kind=="CDS":
     protein="".join(qual.get("translation",[]));locus=(qual.get("locus_tag") or [f"{record['id']}_{start}_{end}"])[0]
     cdss.append((f"{genome}|{locus}",record["id"],start,end,strand,protein,qual))
   if cdss and len(rows)==0:print('sample cdss',len(cdss),[(x[0],len(x[5])) for x in cdss[:3]])
   for gene_id,record_id,start,end,strand,protein,qual in cdss:
    key=hashlib.sha256(protein.encode()).hexdigest();source=seqmap.get(key)
    hit=next(((rs,re_,cluster,product,qual) for rs,re_,cluster,product,qual in regions if start>=rs and end<=re_),None)
    if source is None or hit is None:continue
    rs,re_,cluster,product,qual=hit
    eggrow=egg.get(source)
    rows.append({"gene_id":gene_id,"genome_id":genome,"contig_id":record_id,"start":start,"end":end,
      "strand":strand,"product":"","ko":"","pfam":"","cog":"","eggnog":"","smash_tool":"","cluster_id":cluster,
      "sgb_gene_id":source,"translation":protein})
 # second pass load egg annotations and enrich
 print('seqmap',len(seqmap))
 egg={}
 with open(a.eggnog_features) as fh:
  for row in csv.DictReader(fh,delimiter="\t"):egg[row["gene_id"]]=row
 fields=["gene_id","genome_id","contig_id","start","end","strand","product","ko","pfam","cog","eggnog","smash_tool","cluster_id","sgb_gene_id","matched","translation"]
 with open(a.output,"w",newline="") as out:
  writer=csv.DictWriter(out,fieldnames=fields,delimiter="\t",lineterminator="\n");writer.writeheader();matched=0
  for row in rows:
   eggrow=egg.get(row["sgb_gene_id"])
   row["matched"]=bool(eggrow)
   if eggrow:
    matched+=1
    for field in ["product","ko","pfam","cog","eggnog"]:
     if not row[field]:row[field]=eggrow[field]
   writer.writerow(row)
 print(f"matched {matched}/{len(rows)}")
