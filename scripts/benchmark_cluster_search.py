#!/usr/bin/env python3
"""Synthetic rearrangement and held-out family retrieval benchmarks.

Mode `rearrangement` indexes each explicit source cluster and asks whether
perturbed copies recover it. Mode `holdout` removes selected explicit clusters,
searches them against remaining clusters, and evaluates family accuracy from an
optional cluster_id -> family mapping.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,random,subprocess,sys,tempfile

def stable(value):return int(hashlib.sha256(value.encode()).hexdigest()[:12],16)
from collections import defaultdict
from pathlib import Path
FIELDS=["gene_id","genome_id","contig_id","start","end","strand","product","ko","pfam","cog","eggnog","smash_tool","cluster_id"]
def run(cmd):
 print("+"," ".join(map(str,cmd)),file=sys.stderr);subprocess.run(list(map(str,cmd)),check=True)
def read_features(path):
 with open(path) as fh:return list(csv.DictReader(fh,delimiter="\t"))
def write_features(path,rows):
 with open(path,"w",newline="") as fh:
  w=csv.DictWriter(fh,fieldnames=FIELDS,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def rewrite_cluster(rows,new_id,genome="query",event="source"):
 out=[];positions=0
 for i,row in enumerate(rows):
  length=int(row["end"])-int(row["start"])+1;start=positions;end=start+length-1;positions=end+100
  out.append({**row,"gene_id":f"{event}_{new_id}_{i:03d}","genome_id":genome,"contig_id":"query",
              "start":start,"end":end,"strand":row.get("strand","+"),"cluster_id":new_id})
 if event=="inversion":
  middle=out[len(out)//2:(len(out)//2)+(max(2,len(out)//3))]
  middle.reverse()
  for j,row in enumerate(middle):row["strand"]="+" if row["strand"]!="-" else "-"
  out=out[:len(out)//2]+middle+out[(len(out)//2)+len(middle):]
 elif event=="loss":
  out=out[:-1] if len(out)>3 else out
 elif event=="duplication":
  duplicate=dict(out[0]);duplicate["gene_id"]+="_dup";out.insert(1,duplicate)
 elif event=="truncation":
  out=out[:max(3,int(len(out)*0.7))]
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--features",required=True);ap.add_argument("--output-dir",required=True)
 ap.add_argument("--syn2b",default="target/release/syn2b");ap.add_argument("--mode",choices=["rearrangement","holdout"],default="rearrangement")
 ap.add_argument("--family-map");ap.add_argument("--queries",type=int,default=100);ap.add_argument("--max-gap",type=int,default=5000);ap.add_argument("--seed",type=int,default=20261007)
 ap.add_argument("--min-functional-jaccard",type=float,default=0.25)
 ap.add_argument("--min-synteny",type=float,default=0.10)
 args=ap.parse_args();out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True);rng=random.Random(args.seed)
 rows=read_features(args.features);groups=defaultdict(list)
 for row in rows:
  if row["cluster_id"]:groups[row["cluster_id"]].append(row)
 clusters=sorted(groups)
 if not clusters:raise RuntimeError("no explicit cluster_id values found")
 family={}
 if args.family_map:
  with open(args.family_map) as fh:
   for row in csv.DictReader(fh,delimiter="\t"):family[row["cluster_id"]]=row.get("family") or row.get("cluster_id")
 if args.mode=="rearrangement":
  target=rows;queries=[];truth={}
  for source in clusters[:args.queries]:
   for event in ["loss","duplication","inversion","truncation"]:
    qid=f"{event}__{source}";variant=rewrite_cluster(groups[source],qid,event=event)
    queries.extend(variant);truth[qid]={"truth_cluster":source,"event":event}
  target_path=out/"target_features.tsv";query_path=out/"query_features.tsv";write_features(target_path,target);write_features(query_path,queries)
  s2cl=out/"target.s2cl";index=out/"target.s2cldb";hits=out/"hits.tsv"
  run([args.syn2b,"cluster-build","--features",target_path,"--output",s2cl,"--max-gap",args.max_gap])
  run([args.syn2b,"cluster-index","--cluster-files",s2cl,"--output",index])
  run([args.syn2b,"search-cluster","--query",query_path,"--index",index,"--output",hits,"--min-functional-jaccard",str(args.min_functional_jaccard),"--min-synteny",str(args.min_synteny)])
  hit_rows=defaultdict(list)
  for row in csv.DictReader(open(hits),delimiter="\t"):
   # build_clusters prefixes explicit IDs with genome_id; evaluate by suffix.
   key=row["query_cluster_id"].rsplit(":",1)[-1]
   hit_rows[key].append(row)
  detail=[];by=defaultdict(list)
  for qid,item in truth.items():
   ranked=hit_rows[qid];top=ranked[0] if ranked else None
   correct=bool(top and top["target_cluster_id"].rsplit(":",1)[-1]==item["truth_cluster"])
   row={**item,"query_id":qid,"hit":bool(ranked),"top1_correct":correct,
        "top5_correct":any(x["target_cluster_id"].rsplit(":",1)[-1]==item["truth_cluster"] for x in ranked[:5]),
        "score":float(top["cluster_score"]) if top else 0.0}
   detail.append(row);by[item["event"]].append(row)
  summary=[{"event":"ALL","queries":len(detail),"hit_fraction":sum(x["hit"] for x in detail)/len(detail),
            "top1_accuracy":sum(x["top1_correct"] for x in detail)/len(detail),"top5_accuracy":sum(x["top5_correct"] for x in detail)/len(detail)}]
  for event in sorted(by):
   values=by[event];summary.append({"event":event,"queries":len(values),"hit_fraction":sum(x["hit"] for x in values)/len(values),
                                   "top1_accuracy":sum(x["top1_correct"] for x in values)/len(values),"top5_accuracy":sum(x["top5_correct"] for x in values)/len(values)})
 else:
  if args.family_map:
   by_family=defaultdict(list)
   for cluster_id in clusters:
    family_name=family.get(cluster_id,cluster_id);by_family[family_name].append(cluster_id)
   held=[]
   for family_name in sorted(by_family):
    held.append(sorted(by_family[family_name],key=stable)[0])
    if len(held)>=args.queries:break
  else:held=clusters[:args.queries]
  held_set=set(held);target=[row for row in rows if row["cluster_id"] not in held_set]
  queries=[]
  truth={}
  for source in held:
   qid=f"query__{source}";queries.extend(rewrite_cluster(groups[source],qid));truth[qid]={"truth_cluster":source,"truth_family":family.get(source,source)}
  if not target:raise RuntimeError("holdout mode needs clusters outside the held-out set")
  target_path=out/"target_features.tsv";query_path=out/"query_features.tsv";write_features(target_path,target);write_features(query_path,queries)
  s2cl=out/"target.s2cl";index=out/"target.s2cldb";hits=out/"hits.tsv"
  run([args.syn2b,"cluster-build","--features",target_path,"--output",s2cl,"--max-gap",args.max_gap])
  run([args.syn2b,"cluster-index","--cluster-files",s2cl,"--output",index])
  run([args.syn2b,"search-cluster","--query",query_path,"--index",index,"--output",hits,"--min-functional-jaccard",str(args.min_functional_jaccard),"--min-synteny",str(args.min_synteny)])
  hit_rows=defaultdict(list)
  for row in csv.DictReader(open(hits),delimiter="\t"):
   # build_clusters prefixes explicit IDs with genome_id; evaluate by suffix.
   key=row["query_cluster_id"].rsplit(":",1)[-1]
   hit_rows[key].append(row)
  detail=[]
  for qid,item in truth.items():
   ranked=hit_rows[qid];top=ranked[0] if ranked else None
   predicted=None
   if top:predicted=family.get(top["target_cluster_id"].rsplit(":",1)[-1],top["target_cluster_id"].rsplit(":",1)[-1])
   detail.append({"query_id":qid,"truth_cluster":item["truth_cluster"],"truth_family":item["truth_family"],
                  "predicted_family":predicted or "*","hit":bool(ranked),"top1_family_correct":predicted==item["truth_family"],
                  "top5_family_correct":any((family.get(x["target_cluster_id"].rsplit(":",1)[-1],x["target_cluster_id"].rsplit(":",1)[-1])==item["truth_family"]) for x in ranked[:5]),
                  "score":float(top["cluster_score"]) if top else 0.0})
  summary=[{"event":"HOLDOUT","queries":len(detail),"hit_fraction":sum(x["hit"] for x in detail)/len(detail),
            "top1_accuracy":sum(x["top1_family_correct"] for x in detail)/len(detail),"top5_accuracy":sum(x["top5_family_correct"] for x in detail)/len(detail)}]
 with open(out/"detail.tsv","w",newline="") as fh:
  fields=list(detail[0]);w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(detail)
 with open(out/"summary.tsv","w",newline="") as fh:
  fields=list(summary[0]);w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(summary)
 print(summary)
if __name__=="__main__":main()
