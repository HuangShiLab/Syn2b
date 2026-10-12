#!/usr/bin/env python3
"""Leave-one-source-out, leave-one-contig-out placement benchmark.

Each selected MAG is fragmented at real TGT contig boundaries. All selected
source genomes are excluded from the index; cluster mates and a random OAPGC
background remain available. This is stricter than excluding only one source
at a time and evaluates whether individual contigs recover the correct species
cluster without their source genome.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,os,re,random,subprocess,sys,time
from collections import defaultdict
from pathlib import Path
TAG=re.compile(r"^([A-Za-z0-9_]+):([ACGTacgt]+)@([0-9]+)(?::(.+))?(/[-+])?$")
HEADER=re.compile(r"^>(.*?)\|length=([0-9]+)")

def stable(value):return int(hashlib.sha256(value.encode()).hexdigest()[:12],16)
def run(cmd):
 print("+"," ".join(map(str,cmd)),file=sys.stderr)
 subprocess.run(list(map(str,cmd)),check=True)
def parse_tgt(path):
 genome=None;length=0;contigs=defaultdict(list)
 with open(path) as fh:
  for line in fh:
   if line.startswith(">"):
    m=HEADER.match(line.strip())
    if m:genome,length=m.group(1),int(m.group(2))
   elif line.startswith("#") or not line.strip():continue
   m=TAG.match(line.strip()) if False else None
   for token in line.split():
    m=TAG.match(token)
    if m:
     raw=m.group(4) or "";suffix=m.group(5) or ""
     contig=raw.rsplit("/",1)[0] if raw else "contig"
     contigs[contig].append((int(m.group(3)),m.group(1),m.group(2),suffix))
 for values in contigs.values():values.sort()
 return genome,length,dict(contigs)
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--metadata",default="/Users/macstudio/Downloads/OAPGC/oapgc_genome_metadata.tsv")
 ap.add_argument("--tgt-dir",default="/Users/macstudio/Downloads/OAPGC/census/tgt")
 ap.add_argument("--syn2b",default="target/release/syn2b")
 ap.add_argument("--output-dir",required=True)
 ap.add_argument("--sources",type=int,default=50)
 ap.add_argument("--min-contigs",type=int,default=2)
 ap.add_argument("--background-genomes",type=int,default=5000)
 ap.add_argument("--max-contigs-per-source",type=int,default=5)
 ap.add_argument("--seed",type=int,default=20261007)
 args=ap.parse_args()
 out=Path(args.output_dir);(out/"ref_tgt").mkdir(parents=True,exist_ok=True)
 rng=random.Random(args.seed);rows=[];clusters=defaultdict(list)
 with open(args.metadata) as fh:
  for row in csv.DictReader(fh,delimiter="\t"):
   path=Path(args.tgt_dir)/f"{row['genome_id']}.tgt"
   if not path.exists():continue
   try:comp=float(row["completeness"] or 0);cont=float(row["contamination"] or 0);n50=int(row["n50"] or 0)
   except ValueError:continue
   row.update(path=str(path),comp=comp,cont=cont,n50_i=n50);rows.append(row);clusters[row["cluster"]].append(row)
 eligible=[]
 for cluster,mates in clusters.items():
  high=[r for r in mates if r["comp"]>=95 and r["cont"]<=2 and r["n50_i"]>=50_000]
  if len(high)>=2:eligible.append((cluster,high))
 eligible.sort(key=lambda x:stable(x[0]));selected=eligible[:args.sources]
 if len(selected)<args.sources:raise RuntimeError(f"only {len(selected)} eligible sources")
 source_rows=[];source_ids=set();mate_ids=set()
 for cluster,high in selected:
  # Prefer MAGs with several contigs so this is truly leave-one-contig-out.
  candidates=[]
  for row in high:
   try:_,_,contigs=parse_tgt(row["path"])
   except Exception:continue
   big=[c for c,v in contigs.items() if len(v)>=10]
   if len(big)>=args.min_contigs:candidates.append((len(big),row,contigs))
  candidates.sort(key=lambda x:(-x[0],stable(x[1]["genome_id"])))
  if not candidates:continue
  _,row,contigs=candidates[0];source_rows.append((row,contigs));source_ids.add(row["genome_id"])
  mate_ids.update(r["genome_id"] for r in high if r["genome_id"] not in source_ids)
 candidates=[r for r in rows if r["genome_id"] not in source_ids and r["genome_id"] not in mate_ids]
 candidates.sort(key=lambda r:stable("background:"+r["genome_id"]));background=candidates[:args.background_genomes]
 reference_ids=sorted(mate_ids|{r["genome_id"] for r in background})
 meta={r["genome_id"]:r for r in rows}
 for gid in reference_ids:
  dst=out/"ref_tgt"/f"{gid}.tgt"
  if not dst.exists():os.symlink(meta[gid]["path"],dst)
 with open(out/"reference_metadata.tsv","w",newline="") as fh:
  w=csv.writer(fh,delimiter="\t",lineterminator="\n");w.writerow(["genome_id","species_cluster","taxonomy","fasta_path","is_representative"])
  for gid in reference_ids:
   r=meta[gid];w.writerow([gid,r["cluster"],"",r["path"],"false"])
 truth={};queries=[]
 for source,contigs in source_rows:
  gid=source["genome_id"];truth[gid]=source["cluster"]
  chosen=sorted(contigs,key=lambda c:(-len(contigs[c]),stable(c)))[:args.max_contigs_per_source]
  for ci,contig in enumerate(chosen):
   positions=[x[0] for x in contigs[contig]];start=min(positions);end=max(positions);qlen=end-start+1
   qid=f"loc__{gid}__{ci:03d}";tokens=[]
   for pos,enzyme,seq,suffix in contigs[contig]:
    tokens.append(f"{enzyme}:{seq}@{pos-start}{suffix}")
   queries.append((qid,qlen,tokens))
 query=out/"queries.tgt"
 with query.open("w") as fh:
  for qid,qlen,tokens in queries:fh.write(f">{qid}|length={qlen}\n"+' '.join(tokens)+"\n")
 with open(out/"truth_manifest.tsv","w",newline="") as fh:
  w=csv.writer(fh,delimiter="\t",lineterminator="\n");w.writerow(["query_id","truth_genome","truth_cluster","tags"])
  for qid,qlen,tokens in queries:
   gid=qid.split("__")[1];w.writerow([qid,gid,truth[gid],len(tokens)])
 index=out/"reference.s2cidx"
 started=time.perf_counter();run([args.syn2b,"contig-index","--genome-tgt-dir",out/"ref_tgt","--output",index,"--genome-metadata",out/"reference_metadata.tsv"]);build=time.perf_counter()-started
 placements=out/"placements.tsv";details=out/"details.tsv";hits=placements.with_name("placements.hits.tsv")
 started=time.perf_counter();run([args.syn2b,"place-contig","--index",index,"--query",query,"--output",placements,"--details",details]);search=time.perf_counter()-started
 truth_by_q={r["query_id"]:r for r in csv.DictReader(open(out/"truth_manifest.tsv"),delimiter="\t")}
 hit_rows=defaultdict(list)
 for row in csv.DictReader(open(hits),delimiter="\t"):hit_rows[row["query_id"]].append(row)
 source_groups=defaultdict(list)
 for qid in truth_by_q:source_groups[qid.rsplit("__",1)[0]].append(qid)
 consensus={}
 for source,qids in source_groups.items():
  votes=defaultdict(float)
  for qid in qids:
   ranked=hit_rows[qid]
   if ranked:
    top=ranked[0]
    votes[top["species_cluster"]]+=float(top["score"])*(int(top["shared_adjacencies"])+1)
  consensus_cluster=max(votes, key=lambda cluster:(votes[cluster],cluster)) if votes else "*"
  for qid in qids:consensus[qid]=consensus_cluster
 summary=[];by=defaultdict(list)
 for qid,item in truth_by_q.items():
  ranked=hit_rows[qid];top=ranked[0] if ranked else None
  row={"query_id":qid,"truth_genome":item["truth_genome"],"truth_cluster":item["truth_cluster"],
       "predicted_cluster":top["species_cluster"] if top else "*","hit":bool(ranked),
       "top1_species_correct":bool(top and top["species_cluster"]==item["truth_cluster"]),
       "top5_species_correct":any(r["species_cluster"]==item["truth_cluster"] for r in ranked[:5]),
       "consensus_cluster":consensus[qid],
       "consensus_top1_species_correct":consensus[qid]==item["truth_cluster"]}
  by[item["truth_cluster"]].append(row)
 allrows=[r for values in by.values() for r in values]
 summary.append({"level":"ALL","queries":len(allrows),"hit_fraction":sum(r["hit"] for r in allrows)/len(allrows),
  "top1_species_accuracy":sum(r["top1_species_correct"] for r in allrows)/len(allrows),
  "top5_species_accuracy":sum(r["top5_species_correct"] for r in allrows)/len(allrows),
  "consensus_top1_species_accuracy":sum(r["consensus_top1_species_correct"] for r in allrows)/len(allrows),
  "cluster_level_top1_mean":sum(sum(r["top1_species_correct"] for r in values)/len(values) for values in by.values())/len(by)})
 with open(out/"summary.tsv","w",newline="") as fh:
  fields=list(summary[0]);w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(summary)
 with open(out/"detail.tsv","w",newline="") as fh:
  fields=list(allrows[0]);w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(allrows)
 stats={"reference_genomes":len(reference_ids),"sources":len(source_rows),"queries":len(queries),"index_build_seconds":build,"search_seconds":search,"queries_per_second":len(queries)/search,"index_bytes":index.stat().st_size}
 with open(out/"run_stats.json","w") as fh:json.dump(stats,fh,indent=2) if False else fh.write(json.dumps(stats,indent=2))
 print(stats);print(summary[0])
if __name__=="__main__":main()
