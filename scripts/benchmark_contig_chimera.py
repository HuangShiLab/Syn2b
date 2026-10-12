#!/usr/bin/env python3
"""Simulated chimera benchmark using existing OAPGC TGT records."""
from __future__ import annotations
import argparse,csv,hashlib,json,os,re,random,subprocess,sys,time
from collections import defaultdict
from pathlib import Path
TAG=re.compile(r"^([A-Za-z0-9_]+):([ACGTacgt]+)@([0-9]+)(?::(.+))?(/[-+])?$")
HEADER=re.compile(r"^>(.*?)\|length=([0-9]+)")
def stable(s):return int(hashlib.sha256(s.encode()).hexdigest()[:12],16)
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
   for token in line.split():
    m=TAG.match(token)
    if m:
     raw=m.group(4) or "";suffix=m.group(5) or ""
     contig=raw.rsplit("/",1)[0] if raw else "contig"
     contigs[contig].append((int(m.group(3)),m.group(1),m.group(2),suffix))
 for values in contigs.values():values.sort()
 return genome,length,dict(contigs)
def rewrite(records):
 records=sorted(records);start=records[0][0];end=records[-1][0]
 tokens=[f"{enzyme}:{seq}@{pos-start}{suffix}" for pos,enzyme,seq,suffix in records]
 return end-start+1,tokens,start,end
def link(rows,gids,dst):
 for gid in gids:
  path=dst/f"{gid}.tgt"
  if not path.exists():os.symlink(rows[gid]["path"],path)
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--metadata",default="/Users/macstudio/Downloads/OAPGC/oapgc_genome_metadata.tsv")
 ap.add_argument("--tgt-dir",default="/Users/macstudio/Downloads/OAPGC/census/tgt")
 ap.add_argument("--syn2b",default="target/release/syn2b")
 ap.add_argument("--output-dir",required=True)
 ap.add_argument("--positives",type=int,default=50)
 ap.add_argument("--controls",type=int,default=50)
 ap.add_argument("--left-min",type=int,default=250_000)
 ap.add_argument("--left-max",type=int,default=500_000)
 ap.add_argument("--right-min",type=int,default=250_000)
 ap.add_argument("--right-max",type=int,default=500_000)
 ap.add_argument("--background-genomes",type=int,default=2_000)
 ap.add_argument("--window",default="50kb")
 ap.add_argument("--seed",type=int,default=20261007)
 args=ap.parse_args();out=Path(args.output_dir);ref=out/"ref_tgt";ref.mkdir(parents=True,exist_ok=True)
 rows={};clusters=defaultdict(list)
 with open(args.metadata) as fh:
  for raw in csv.DictReader(fh,delimiter="\t"):
   gid=raw["genome_id"];path=Path(args.tgt_dir)/f"{gid}.tgt"
   try:comp=float(raw["completeness"] or 0);cont=float(raw["contamination"] or 0);n50=int(raw["n50"] or 0)
   except ValueError:continue
   if not path.exists():continue
   row=dict(raw);row.update(path=str(path),comp=comp,cont=cont,n50=n50);rows[gid]=row;clusters[raw["cluster"]].append(row)
 eligible=[r for r in rows.values() if r["comp"]>=95 and r["cont"]<=2 and r["n50"]>=100_000]
 eligible.sort(key=lambda r:stable(r["genome_id"]))
 need=2*args.positives+args.controls
 selected=eligible[:need]
 if len(selected)<need:raise RuntimeError(f"only {len(selected)} eligible genomes; need {need}")
 truth=[];ref_ids=set()
 for i in range(args.positives):
  a=selected[2*i];b=selected[2*i+1]
  _,_,ca=parse_tgt(a["path"]);_,_,cb=parse_tgt(b["path"])
  ca=sorted([values for values in ca.values() if len(values)>=10],key=len)
  cb=sorted([values for values in cb.values() if len(values)>=10],key=len)
  if not ca or not cb:continue
  ca=ca[-1];cb=cb[-1]
  ll=random.Random(args.seed+i*2).randint(args.left_min,args.left_max)
  rl=random.Random(args.seed+i*2+1).randint(args.right_min,args.right_max)
  left=[record for record in ca if record[0]<ca[0][0]+ll]
  right=[record for record in cb if record[0]>=cb[-1][0]-rl]
  if len(left)<5 or len(right)<5:continue
  left_length,left_tokens,left_start,left_end=rewrite(left)
  right_length,right_tokens,right_start,right_end=rewrite(right)
  breakpoint=left_end-left_start+1
  qid=f"chim__{a['genome_id']}__{b['genome_id']}__{i:04d}"
  right_offset=breakpoint
  shifted_right=[f"{token.split('@')[0]}@{int(token.split('@')[1].split('/')[0])+right_offset}{('/'+token.split('/')[-1]) if '/' in token else ''}" for token in right_tokens]
  query_tokens=left_tokens+shifted_right
  truth.append({"query_id":qid,"truth":"chimera","left_genome":a["genome_id"],"right_genome":b["genome_id"],
                "left_cluster":a["cluster"],"right_cluster":b["cluster"],"breakpoint":breakpoint,
                "query_length":left_length+right_end-right_start+1,"tags":len(query_tokens),"tokens":query_tokens})
  ref_ids.update((a["genome_id"],b["genome_id"]))
 for i in range(args.controls):
  row=selected[-1-i];_,_,contigs=parse_tgt(row["path"])
  large=sorted([values for values in contigs.values() if len(values)>=10],key=len)
  if not large:continue
  qlen,tokens,start,end=rewrite(large[-1])
  qid=f"control__{row['genome_id']}__{i:04d}"
  truth.append({"query_id":qid,"truth":"control","cluster":row["cluster"],"query_length":qlen,"tags":len(tokens),"tokens":tokens})
  ref_ids.add(row["genome_id"])
 if not truth:raise RuntimeError("no chimera/control queries generated")
 link(rows,ref_ids,ref)
 candidates=[row for gid,row in rows.items() if gid not in ref_ids and row["comp"]>=90]
 candidates.sort(key=lambda row:stable("background:"+row["genome_id"]))
 link(rows,[row["genome_id"] for row in candidates[:args.background_genomes]],ref)
 query=out/"queries.tgt"
 with query.open("w") as fh:
  for item in truth:
   fh.write(f">{item['query_id']}|length={item['query_length']}\n")
   fh.write(' '.join(item["tokens"])+"\n")
 with open(out/"truth_manifest.tsv","w",newline="") as fh:
  fields=sorted({key for row in truth for key in row});writer=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n")
  writer.writeheader();writer.writerows(truth)
 index=out/"reference.s2cidx";started=time.perf_counter()
 run([args.syn2b,"contig-index","--genome-tgt-dir",ref,"--output",index]);build_seconds=time.perf_counter()-started
 chimerism=out/"chimerism.tsv";started=time.perf_counter()
 run([args.syn2b,"scan-chimera","--index",index,"--query",query,"--output",chimerism,"--window",args.window])
 scan_seconds=time.perf_counter()-started
 segments=defaultdict(list)
 for row in csv.DictReader(open(chimerism),delimiter="\t"):
  # scan output appends :start-end to the query contig name.
  segments[row["query_id"].rsplit(":",1)[0]].append(row)
 rows=[]
 for item in truth:
  window_segments=segments.get(item["query_id"],[]);transitions=[]
  for segment in window_segments:
   if segment["cluster_transition"].lower()=="true":
    start=int(segment["start"]);end=int(segment["end"]);transitions.append((start+end)//2)
  detected=bool(transitions)
  error=min((abs(position-item.get("breakpoint",0)) for position in transitions),default=None)
  rows.append({**item,"detected_chimera":detected,"observed_breakpoint":transitions[0] if transitions else "NA",
               "breakpoint_error":error if error is not None else "NA"})
 positives=[row for row in rows if row["truth"]=="chimera"];controls=[row for row in rows if row["truth"]=="control"]
 errors=[row["breakpoint_error"] for row in positives if row["breakpoint_error"]!="NA"]
 stats={"reference_genomes":len(list(ref.glob('*.tgt'))),"positives":len(positives),"controls":len(controls),
        "sensitivity":sum(row["detected_chimera"] for row in positives)/max(1,len(positives)),
        "specificity":1-sum(row["detected_chimera"] for row in controls)/max(1,len(controls)),
        "false_positive_rate":sum(row["detected_chimera"] for row in controls)/max(1,len(controls)),
        "median_breakpoint_error_bp":sorted(errors)[len(errors)//2] if errors else "NA",
        "index_build_seconds":build_seconds,"scan_seconds":scan_seconds,"index_bytes":index.stat().st_size}
 with open(out/"detail.tsv","w",newline="") as fh:
  fields=sorted({key for row in rows for key in row});writer=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n");writer.writeheader();writer.writerows(rows)
 with open(out/"run_stats.json","w") as fh:fh.write(json.dumps(stats,indent=2))
 print(json.dumps(stats,indent=2))
if __name__=="__main__":main()
