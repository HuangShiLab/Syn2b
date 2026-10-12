#!/usr/bin/env python3
"""Index-build and query runtime benchmark at a requested reference scale."""
from __future__ import annotations
import argparse,csv,hashlib,json,os,re,subprocess,time
from pathlib import Path
def stable(value):return int(hashlib.sha256(value.encode()).hexdigest()[:12],16)
def run(cmd):
 print("+"," ".join(map(str,cmd)),file=sys.stderr);subprocess.run(list(map(str,cmd)),check=True)
def timed(cmd,log):
 cmd=["/usr/bin/time","-l","-o",str(log)]+list(map(str,cmd));started=time.perf_counter()
 subprocess.run(list(map(str,cmd)),check=True);wall=time.perf_counter()-started
 text=Path(log).read_text();rss=None
 for pattern in [r"([0-9]+)\s+maximum resident set size",r"([0-9]+)\s+peak memory footprint"]:
  rss=re.search(pattern,text)
  if rss:break
 real=re.search(r"([0-9.]+) real",text)
 return {"wall_seconds":wall,"real_seconds":float(real.group(1)) if real else wall,"max_rss_bytes":int(rss.group(1)) if rss else 0}
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--metadata",default="/Users/macstudio/Downloads/OAPGC/oapgc_genome_metadata.tsv")
 ap.add_argument("--tgt-dir",default="/Users/macstudio/Downloads/OAPGC/census/tgt");ap.add_argument("--syn2b",default="target/release/syn2b")
 ap.add_argument("--query-tgt",required=True);ap.add_argument("--output-dir",required=True);ap.add_argument("--target-genomes",type=int,default=10_000)
 ap.add_argument("--heldout-manifest");ap.add_argument("--min-completeness",type=float,default=90);ap.add_argument("--max-contamination",type=float,default=5)
 args=ap.parse_args();out=Path(args.output_dir);ref=out/"ref_tgt";ref.mkdir(parents=True,exist_ok=True)
 rows=[];held=set()
 if args.heldout_manifest:
  with open(args.heldout_manifest) as fh:
   for row in csv.DictReader(fh,delimiter="\t"):
    if "truth_genome" in row:held.add(row["truth_genome"])
 with open(args.metadata) as fh:
  for row in csv.DictReader(fh,delimiter="\t"):
   gid=row["genome_id"];path=Path(args.tgt_dir)/f"{gid}.tgt"
   try:comp=float(row["completeness"] or 0);cont=float(row["contamination"] or 0)
   except ValueError:continue
   if path.exists() and gid not in held and comp>=args.min_completeness and cont<=args.max_contamination:rows.append((gid,path))
 rows.sort(key=lambda pair:stable(pair[0]));selected=rows[:args.target_genomes]
 if len(selected)<args.target_genomes:raise RuntimeError(f"only {len(selected)} genomes available")
 for gid,path in selected:
  link=ref/f"{gid}.tgt"
  if not link.exists():os.symlink(path,link)
 metadata=out/"reference_metadata.tsv"
 with metadata.open("w",newline="") as fh:
  writer=csv.writer(fh,delimiter="\t",lineterminator="\n");writer.writerow(["genome_id","species_cluster","taxonomy","fasta_path","is_representative"])
  for gid,path in selected:writer.writerow([gid,gid,"",str(path),"false"])
 index=out/"reference.s2cidx";build=timed([args.syn2b,"contig-index","--genome-tgt-dir",ref,"--output",index,"--genome-metadata",metadata],out/"time_index.log")
 query=timed([args.syn2b,"place-contig","--index",index,"--query",args.query_tgt,"--output",out/"placements.tsv","--details",out/"details.tsv"],out/"time_search.log")
 query_count=sum(1 for line in open(args.query_tgt) if line.startswith(">"))
 stats={"target_genomes_requested":args.target_genomes,"reference_genomes":len(selected),"queries":query_count,
        "index_bytes":index.stat().st_size,**{f"index_{k}":v for k,v in build.items()},**{f"search_{k}":v for k,v in query.items()},
        "queries_per_second":query_count/query["wall_seconds"]}
 with open(out/"run_stats.json","w") as fh:fh.write(json.dumps(stats,indent=2))
 print(json.dumps(stats,indent=2))
if __name__=="__main__":main()
