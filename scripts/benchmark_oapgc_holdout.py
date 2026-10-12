#!/usr/bin/env python3
"""Leave-one-genome-out placement validation using prebuilt OAPGC TGTs.

For each selected species cluster, one high-quality genome is held out and
fragmented at the landmark-coordinate level. Its cluster mates plus a random
background of unrelated OAPGC genomes are symlinked into a reference TGT
directory. The held-out genome is never indexed.

This measures species-cluster placement (the handoff's Benchmark B), not
top-1 recovery of the absent source genome.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, random, re, subprocess, sys, time
from collections import defaultdict
from pathlib import Path

TAG = re.compile(r"^([A-Za-z0-9_]+):([ACGTacgt]+)@([0-9]+)(?::(.+))?(/[-+])?$")
HEADER = re.compile(r"^>(.*?)\|length=([0-9]+)")
LEVELS = [500_000, 250_000, 100_000, 50_000, 25_000, 10_000]

def stable(s): return int(hashlib.sha256(s.encode()).hexdigest()[:12], 16)
def run_timed(cmd, log):
    cmd = ["/usr/bin/time", "-l", "-o", str(log)] + list(map(str, cmd))
    started = time.perf_counter()
    cp = subprocess.run(cmd, text=True, capture_output=True)
    elapsed = time.perf_counter() - started
    if cp.returncode: raise RuntimeError(f"command failed: {' '.join(map(str,cmd))}\n{cp.stdout}\n{cp.stderr}")
    return elapsed

def time_stats(log):
    text = Path(log).read_text()
    real = re.search(r"\s([0-9.]+) real", text)
    rss = re.search(r"\s([0-9]+) maximum resident set size", text)
    return float(real.group(1)) if real else 0.0, int(rss.group(1)) if rss else 0

def parse_tags(path):
    """Return genome, length, and per-contig landmark records."""
    genome=None; length=0; by_contig=defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                m=HEADER.match(line.strip())
                if m: genome,length=m.group(1),int(m.group(2))
            elif line.startswith("#") or not line.strip(): continue
            for token in line.split():
                m=TAG.match(token)
                if not m: continue
                enzyme,seq,pos=m.group(1),m.group(2),int(m.group(3))
                suffix=m.group(5) or ""
                raw=m.group(4) or ""
                contig=raw.rsplit("/",1)[0] if raw else "contig"
                by_contig[contig].append((pos,enzyme,seq,suffix))
    for contig in by_contig: by_contig[contig].sort()
    return genome,length,dict(by_contig)

def fragment(path, level, max_windows):
    genome,length,contigs=parse_tags(path)
    windows=[]
    for contig,tags in sorted(contigs.items()):
        if len(tags)<10: continue
        start0=tags[0][0]; end0=tags[-1][0]; span=end0-start0+1
        if span<min(10_000,level): continue
        step=min(level,span)
        these=[]
        for start in range(start0,start0+span-step+1,step):
            end=start+step-1
            selected=[(pos,enzyme,seq,suffix) for pos,enzyme,seq,suffix in tags if start<=pos<=end]
            if len(selected)>=5: these.append((contig,start,end,selected))
        if not these: continue
        if len(these)<=max_windows: selected_windows=these
        else:
            indexes=sorted(set(round(i*(len(these)-1)/(max_windows-1)) for i in range(max_windows)))
            selected_windows=[these[i] for i in indexes]
        windows.extend(selected_windows)
    return genome,length,windows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--metadata",default="/Users/macstudio/Downloads/OAPGC/oapgc_genome_metadata.tsv")
    ap.add_argument("--tgt-dir",default="/Users/macstudio/Downloads/OAPGC/census/tgt")
    ap.add_argument("--syn2b",default="target/release/syn2b")
    ap.add_argument("--output-dir",required=True)
    ap.add_argument("--target-clusters",type=int,default=30)
    ap.add_argument("--include-sources",action="store_true",help="index held-out source genomes for Benchmark-A fragmentation recovery")
    ap.add_argument("--background-genomes",type=int,default=970)
    ap.add_argument("--max-windows-per-level",type=int,default=2)
    ap.add_argument("--min-cluster-mates",type=int,default=2)
    ap.add_argument("--seed",type=int,default=20261006)
    ap.add_argument("--levels",type=int,nargs="*",default=LEVELS)
    args=ap.parse_args()
    out=Path(args.output_dir); (out/"ref_tgt").mkdir(parents=True,exist_ok=True)
    rng=random.Random(args.seed)

    rows=[]; by_cluster=defaultdict(list)
    with open(args.metadata) as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            gid=row["genome_id"]; path=Path(args.tgt_dir)/f"{gid}.tgt"
            try: comp=float(row["completeness"] or 0); cont=float(row["contamination"] or 0); n50=int(row["n50"] or 0)
            except ValueError: continue
            if not path.exists(): continue
            row.update(path=str(path),completeness_f=comp,contamination_f=cont,n50_i=n50)
            rows.append(row); by_cluster[row["cluster"]].append(row)

    eligible=[]
    for cluster,mates in by_cluster.items():
        high=[r for r in mates if r["completeness_f"]>=95 and r["contamination_f"]<=2 and r["n50_i"]>=50_000]
        if len(high)>=args.min_cluster_mates+1:
            eligible.append((cluster,high))
    eligible.sort(key=lambda x: stable(x[0]))
    selected=eligible[:args.target_clusters]
    if len(selected)<args.target_clusters: raise RuntimeError(f"only {len(selected)} eligible clusters")

    targets=[]; target_ids=set(); cluster_mate_ids=set()
    for cluster,high in selected:
        hold=sorted(high,key=lambda r:(-r["n50_i"],stable(r["genome_id"])))[0]
        targets.append(hold); target_ids.add(hold["genome_id"])
        if args.include_sources:
            cluster_mate_ids.add(hold["genome_id"])
        cluster_mate_ids.update(r["genome_id"] for r in high if r["genome_id"]!=hold["genome_id"])

    candidates=[r for r in rows if r["genome_id"] not in cluster_mate_ids]
    candidates.sort(key=lambda r: stable("background:"+r["genome_id"]))
    background=candidates[:args.background_genomes]
    reference_ids=sorted(cluster_mate_ids | {r["genome_id"] for r in background})
    for gid in reference_ids:
        src=Path(args.tgt_dir)/f"{gid}.tgt"; dst=out/"ref_tgt"/f"{gid}.tgt"
        if not dst.exists(): os.symlink(src,dst)

    with open(out/"reference_metadata.tsv","w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n")
        w.writerow(["genome_id","species_cluster","taxonomy","fasta_path","is_representative"])
        meta={r["genome_id"]:r for r in rows}
        for gid in reference_ids:
            r=meta[gid]; w.writerow([gid,r["cluster"],"",str(r["path"]),"false"])

    truth={r["genome_id"]:r["cluster"] for r in targets}
    query_records=[]; manifest=[]
    for target in targets:
        gid=target["genome_id"]
        for level in args.levels:
            packed=fragment(Path(target["path"]),level,args.max_windows_per_level)
            if not packed: continue
            _,length,windows=packed
            for qi,(contig,start,end,selected) in enumerate(windows):
                qid=f"frag__{gid}__{level}__{qi:04d}"
                qlen=end-start+1
                tokens=[]
                for pos,enzyme,seq,suffix in selected:
                    tokens.append(f"{enzyme}:{seq}@{pos-start}{suffix}")
                query_records.append((qid,qlen,tokens))
                manifest.append({"query_id":qid,"truth_genome":gid,"truth_cluster":truth[gid],
                                 "level":level,"contig":contig,"start":start,"end":end,"tags":len(tokens)})
    query=out/"queries.tgt"
    with query.open("w") as fh:
        for qid,qlen,tokens in query_records:
            fh.write(f">{qid}|length={qlen}\n")
            fh.write(" ".join(tokens)+"\n")
    with open(out/"truth_manifest.tsv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(manifest[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(manifest)

    index=out/"reference.s2cidx"
    if index.exists():
        # Reuse an index when iterating quickly on query construction.
        index_elapsed, index_rss = time_stats(out/"time_index.log")
    else:
        index_elapsed=run_timed([args.syn2b,"contig-index","--genome-tgt-dir",out/"ref_tgt",
                                 "--output",index,"--genome-metadata",out/"reference_metadata.tsv"],out/"time_index.log")
        _,index_rss=time_stats(out/"time_index.log")
    placements=out/"placements.tsv"; details=out/"details.tsv"
    hits=placements.with_name(placements.stem + ".hits.tsv")
    search_elapsed=run_timed([args.syn2b,"place-contig","--index",index,"--query",query,
                              "--output",placements,"--details",details,
                              "--min-tags",5,"--min-adjacencies",3,"--max-targets",25],out/"time_search.log")
    _,search_rss=time_stats(out/"time_search.log")

    with open(hits) as fh: hit_rows=list(csv.DictReader(fh,delimiter="\t"))
    hits_by_q=defaultdict(list)
    for row in hit_rows: hits_by_q[row["query_id"]].append(row)
    with open(placements) as fh: placement_rows={r["query_id"]:r for r in csv.DictReader(fh,delimiter="\t")}
    detail=[]
    for item in manifest:
        qid=item["query_id"]; ranked=hits_by_q[qid]
        top=ranked[0] if ranked else None
        top_cluster=top["species_cluster"] if top else "*"
        top5=[r["species_cluster"] for r in ranked[:5]]
        placement=placement_rows[qid]
        detail.append({**item,
            "predicted_genome":top["genome"] if top else "*", "predicted_cluster":top_cluster,
            "top1_genome_correct":top is not None and top["genome"]==item["truth_genome"],
            "top1_species_correct":top_cluster==item["truth_cluster"],
            "top5_species_correct":item["truth_cluster"] in top5,
            "orientation":top["orientation"] if top else "*", "forward_orientation":(top or {}).get("orientation")=="+" ,
            "rank_true_cluster":next((i+1 for i,c in enumerate([r["species_cluster"] for r in ranked]) if c==item["truth_cluster"]),"NA"),
            "top_score":float(top["score"]) if top else 0.0,
            "placement_score":float(placement["placement_score"]), "confidence":placement["confidence"],
            "query_coverage":float(placement["query_coverage"]),
            "confident_and_correct":placement["confidence"] in ("high","medium") and top_cluster==item["truth_cluster"]})
    by=defaultdict(list)
    for row in detail: by[row["level"]].append(row)
    summary=[]
    all_true=[r["top1_species_correct"] for r in detail]; all_top5=[r["top5_species_correct"] for r in detail]
    all_genome=[r["top1_genome_correct"] for r in detail]
    all_conf=[r["confident_and_correct"] for r in detail]; all_cov=[r["query_coverage"] for r in detail]
    summary.append({"level_bp":"ALL","queries":len(detail),"top1_species_cluster_accuracy":sum(all_true)/len(detail),
        "top5_species_cluster_accuracy":sum(all_top5)/len(detail),"top1_genome_accuracy":sum(all_genome)/len(detail),"median_true_cluster_rank":sorted(float(r["rank_true_cluster"]) for r in detail if r["rank_true_cluster"]!="NA")[len([r for r in detail if r["rank_true_cluster"]!="NA"])//2],
        "forward_orientation_fraction":sum(r["forward_orientation"] for r in detail)/len(detail),
        "mean_query_coverage":sum(all_cov)/len(detail),"confident_correct_fraction":sum(all_conf)/len(detail)})
    for level in sorted(by):
        rs=by[level]; ranks=[float(r["rank_true_cluster"]) for r in rs if r["rank_true_cluster"]!="NA"]
        summary.append({"level_bp":level,"queries":len(rs),
            "top1_species_cluster_accuracy":sum(r["top1_species_correct"] for r in rs)/len(rs),
            "top5_species_cluster_accuracy":sum(r["top5_species_correct"] for r in rs)/len(rs),
            "top1_genome_accuracy":sum(r["top1_genome_correct"] for r in rs)/len(rs),
            "median_true_cluster_rank":sorted(ranks)[len(ranks)//2] if ranks else "NA",
            "forward_orientation_fraction":sum(r["forward_orientation"] for r in rs)/len(rs),
            "mean_query_coverage":sum(r["query_coverage"] for r in rs)/len(rs),
            "confident_correct_fraction":sum(r["confident_and_correct"] for r in rs)/len(rs)})
    fields=list(summary[0])
    with open(out/"summary.tsv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(summary)
    fields=list(detail[0])
    with open(out/"detail.tsv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(detail)
    stats={"reference_genomes":len(reference_ids),"target_clusters":len(targets),"queries":len(detail),
           "index_build_seconds_wall":index_elapsed,"index_build_real_seconds":time_stats(out/"time_index.log")[0],
           "index_build_max_rss_bytes":index_rss,"search_seconds_wall":search_elapsed,
           "search_real_seconds":time_stats(out/"time_search.log")[0],"search_max_rss_bytes":search_rss,
           "queries_per_second":len(detail)/search_elapsed,"index_bytes":index.stat().st_size,
           "query_tgt_bytes":query.stat().st_size}
    with open(out/"run_stats.json","w") as fh: json.dump(stats,fh,indent=2)
    print(json.dumps(stats,indent=2)); print(summary[0])
if __name__=="__main__": main()
