#!/usr/bin/env python3
"""Generate synthetic gene-cluster databases and measure index/query runtime."""
from __future__ import annotations
import argparse,csv,json,subprocess,time
from pathlib import Path
FIELDS=["gene_id","genome_id","contig_id","start","end","strand","product","ko","pfam","cog","eggnog","smash_tool","cluster_id"]
def run(cmd):
 print("+"," ".join(map(str,cmd)),file=__import__("sys").stderr);subprocess.run(list(map(str,cmd)),check=True)
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--output-dir",required=True);ap.add_argument("--clusters",type=int,default=100_000);ap.add_argument("--genes",type=int,default=5);ap.add_argument("--query-clusters",type=int,default=1_000);ap.add_argument("--syn2b",default="target/release/syn2b")
 ap.add_argument("--min-functional-jaccard",type=float,default=0.5)
 ap.add_argument("--min-synteny",type=float,default=0.25)
 args=ap.parse_args();out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
 target=out/"target_features.tsv";query=out/"query_features.tsv"
 with target.open("w",newline="") as fh:
  w=csv.DictWriter(fh,fieldnames=FIELDS,delimiter="\t",lineterminator="\n");w.writeheader()
  for cluster in range(args.clusters):
   genome=f"genome_{cluster:07d}";cid=f"C{cluster:07d}"
   for gene in range(args.genes):
    start=gene*1000
    w.writerow({"gene_id":f"{cid}_g{gene}","genome_id":genome,"contig_id":"ctg_1","start":start,"end":start+900,
      "strand":"+","product":f"role_{gene}","ko":f"K{cluster:07d}_{gene:02d}","pfam":f"PF{cluster:07d}_{gene:02d}","cog":f"COG{gene:05d}",
      "eggnog":f"ENOG{cluster:07d}_{gene:02d}","smash_tool":"PKS" if gene%2==0 else "T3SS","cluster_id":cid})
 with query.open("w",newline="") as fh:
  w=csv.DictWriter(fh,fieldnames=FIELDS,delimiter="\t",lineterminator="\n");w.writeheader()
  for cluster in range(args.query_clusters):
   genome="query";cid=f"Q{cluster:07d}"
   for gene in range(args.genes):
    start=gene*1000
    w.writerow({"gene_id":f"{cid}_g{gene}","genome_id":genome,"contig_id":"ctg_1","start":start,"end":start+900,
      "strand":"+","product":f"role_{gene}","ko":f"K{cluster:07d}_{gene:02d}","pfam":f"PF{cluster:07d}_{gene:02d}","cog":f"COG{gene:05d}",
      "eggnog":f"ENOG{cluster:07d}_{gene:02d}","smash_tool":"PKS" if gene%2==0 else "T3SS","cluster_id":cid})
 s2cl=out/"target.s2cl";index=out/"clusters.s2cldb";hits=out/"hits.tsv"
 started=time.perf_counter();run([args.syn2b,"cluster-build","--features",target,"--output",s2cl]);build_clusters=time.perf_counter()-started
 started=time.perf_counter();run([args.syn2b,"cluster-index","--cluster-files",s2cl,"--output",index]);build_index=time.perf_counter()-started
 started=time.perf_counter();run([args.syn2b,"search-cluster","--query",query,"--index",index,"--output",hits,"--min-functional-jaccard",str(args.min_functional_jaccard),"--min-synteny",str(args.min_synteny)]);search=time.perf_counter()-started
 stats={"clusters":args.clusters,"genes_per_cluster":args.genes,"feature_rows":args.clusters*args.genes,
        "query_clusters":args.query_clusters,"cluster_build_seconds":build_clusters,
        "index_build_seconds":build_index,"search_seconds":search,"queries_per_second":args.query_clusters/search,
        "index_bytes":index.stat().st_size,"hits":sum(1 for _ in open(hits))-1}
 with open(out/"run_stats.json","w") as fh:fh.write(json.dumps(stats,indent=2))
 print(json.dumps(stats,indent=2))
if __name__=="__main__":main()
