#!/usr/bin/env python3
"""Generate clean FASTA queries/references and run skani/minimap2 baselines.

OAPGC FASTAs contain '-' padding/gap characters. Syn2b TGT coordinates refer to
these original records. For a valid nucleotide baseline this script strips '-'
from references and from each query's original contig interval before running
skani and minimap2.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,os,re,shutil,subprocess,sys,time
from collections import defaultdict
from pathlib import Path
def stable(s):return int(hashlib.sha256(s.encode()).hexdigest()[:12],16)
def run(cmd,stdout=None):
 print("+"," ".join(map(str,cmd)),file=sys.stderr)
 with open(stdout,"wb") if stdout else open(os.devnull,"wb") as out:
  subprocess.run(list(map(str,cmd)),stdout=out,check=True)
def fasta(path):
 name=None;parts=[]
 with open(path) as fh:
  for line in fh:
   if line.startswith(">"):
    if name is not None:yield name,"".join(parts)
    name=line[1:].strip().split()[0];parts=[]
   else:parts.append(line.strip())
 if name is not None:yield name,"".join(parts)
def clean_sequence(seq):return re.sub(r"[^ACGTacgt]","N",seq,flags=re.I)
def contig_offsets(sequences):
 offsets={};offset=0
 for name,sequence in sequences.items():
  offsets[name]=offset;offset+=len(sequence)
 return offsets
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--reference-tgt-dir",required=True,help="directory of symlinks selecting the reference genomes")
 ap.add_argument("--query-truth-manifest",required=True)
 ap.add_argument("--original-fasta-dir",default="/Users/macstudio/Downloads/OAPGC/genomes")
 ap.add_argument("--metadata",default="/Users/macstudio/Downloads/OAPGC/oapgc_genome_metadata.tsv")
 ap.add_argument("--output-dir",required=True)
 ap.add_argument("--clean-fasta-dir",help="directory containing pre-extracted {genome}.clean.fna files")
 ap.add_argument("--threads",type=int,default=8)
 ap.add_argument("--tools",nargs="*",default=["skani","minimap2"])
 args=ap.parse_args();out=Path(args.output_dir);ref=out/"clean_refs";query=out/"queries.clean.fna";ref.mkdir(parents=True,exist_ok=True)
 source_contigs=defaultdict(dict)
 manifest=list(csv.DictReader(open(args.query_truth_manifest),delimiter="\t"))
 for row in manifest:
  gid=row["truth_genome"]
  if gid not in source_contigs:
   d={}
   source_path=Path(args.clean_fasta_dir)/f"{gid}.clean.fna" if args.clean_fasta_dir else Path(args.original_fasta_dir)/f"{gid}.fna"
   for name,seq in fasta(source_path):d[name]=seq
   source_contigs[gid]=(d,contig_offsets(d))
 with query.open("w") as outfh:
  for row in manifest:
   raw_name=row["contig"].split("|",1)[-1];sequences,offsets=source_contigs[row["truth_genome"]]
   seq=sequences.get(raw_name)
   if seq is None:raise RuntimeError(f"missing query contig {raw_name}")
   offset=offsets.get(raw_name,0);start=max(0,int(row["start"])-offset);end=min(len(seq)-1,int(row["end"])-offset)
   if start>end:continue
   seq=seq[start:end+1].upper()
   if args.clean_fasta_dir:seq=re.sub(r"[^ACGTacgt]","N",seq)
   else:seq=clean_sequence(seq)
   if not seq:continue
   outfh.write(f">{row['query_id']}\n{seq}\n")
 selected=[Path(p).resolve().stem for p in Path(args.reference_tgt_dir).glob("*.tgt")];selected.sort()
 meta={}
 with open(args.metadata) as fh:
  for row in csv.DictReader(fh,delimiter="\t"):meta[row["genome_id"]]=row
 ref_list=out/"clean_ref_list.txt"
 if args.clean_fasta_dir:
  written=[str(path) for path in sorted(Path(args.clean_fasta_dir).glob("*.clean.fna"))]
 else:
  written=[]
  for gid in selected:
   clean=ref/f"{gid}.clean.fna"
   if not clean.exists():
    with clean.open("w") as outfh:
     for name,seq in fasta(Path(args.original_fasta_dir)/f"{gid}.fna"):
      seq=clean_sequence(seq)
      if seq:outfh.write(f">{name}\n{seq}\n")
   written.append(str(clean))
 ref_list.write_text("\n".join(written)+"\n")
 results=[]
 if "skani" in args.tools and (out/"skani.raw.tsv").exists():
  raw=out/"skani.raw.tsv";sketch_seconds=0.0;search_seconds=0.0
  hits=defaultdict(list)
  for row in csv.DictReader(open(raw),delimiter="\t"):
   query_id=Path(row["Query_name"]).stem;reference=Path(row["Ref_file"]).stem.replace(".clean","")
   hits[query_id].append((reference,float(row["ANI"])))
  results.append({"tool":"skani","sketch_seconds":sketch_seconds,"search_seconds":search_seconds,"hits":hits})
 elif "skani" in args.tools:
  db=out/"skani_db"
  if db.exists():shutil.rmtree(db)
  started=time.perf_counter()
  run(["skani","sketch","-l",ref_list,"-o",db,"-t",args.threads],stdout=out/"skani_sketch.log")
  sketch_seconds=time.perf_counter()-started
  raw=out/"skani.raw.tsv";started=time.perf_counter()
  run(["skani","search","-d",db,"-q",query,"--qi","-n","25","-t",args.threads],stdout=raw)
  search_seconds=time.perf_counter()-started
  hits=defaultdict(list)
  for row in csv.DictReader(open(raw),delimiter="\t"):
   query_id=Path(row["Query_name"]).stem;reference=Path(row["Ref_file"]).stem.replace(".clean","")
   hits[query_id].append((reference,float(row["ANI"])))
  results.append({"tool":"skani","sketch_seconds":sketch_seconds,"search_seconds":search_seconds,"hits":hits})
 if "minimap2" in args.tools and (out/"minimap2.paf").exists():
  raw=out/"minimap2.paf";concat_seconds=0.0;search_seconds=0.0
  scores=defaultdict(lambda:defaultdict(float))
  for line in open(raw):
   fields=line.rstrip().split("\t")
   if len(fields)<10:continue
   scores[fields[0]][re.sub(r"_ctg_\d+$","",fields[5])]+=int(fields[10])
  hits={q:[(ref,score) for ref,score in sorted(items.items(),key=lambda pair:(-pair[1],pair[0]))] for q,items in scores.items()}
  results.append({"tool":"minimap2","concat_seconds":concat_seconds,"search_seconds":search_seconds,"hits":hits})
 elif "minimap2" in args.tools:
  concat=out/"reference.clean.concat.fna";started=time.perf_counter()
  with concat.open("wb") as outfh:
   for path in written:
    with open(path,"rb") as src:outfh.write(src.read())
  concat_seconds=time.perf_counter()-started
  raw=out/"minimap2.paf";started=time.perf_counter()
  run(["minimap2","-t",str(args.threads),"-x","asm20",concat,query],stdout=raw)
  search_seconds=time.perf_counter()-started
  scores=defaultdict(lambda:defaultdict(float))
  for line in open(raw):
   fields=line.rstrip().split("\t")
   if len(fields)<10:continue
   query_id=fields[0];reference=re.sub(r"_ctg_\d+$","",fields[5]);scores[query_id][reference]+=int(fields[10])
  hits={q:[(ref,score) for ref,score in sorted(items.items(),key=lambda pair:(-pair[1],pair[0]))] for q,items in scores.items()}
  results.append({"tool":"minimap2","concat_seconds":concat_seconds,"search_seconds":search_seconds,"hits":hits})
 truth={row["query_id"]:row for row in manifest}
 clusters={}
 with open(args.metadata) as fh:
  for row in csv.DictReader(fh,delimiter="\t"):clusters[row["genome_id"]]=row.get("cluster") or row["genome_id"]
 summary=[]
 for result in results:
  ranked={query_id:[genome for genome,score in values] for query_id,values in result.pop("hits").items()}
  result["queries"]=len(truth);result["hit_fraction"]=sum(bool(ranked.get(q)) for q in truth)/len(truth)
  result["top1_genome_accuracy"]=sum(next(iter(ranked.get(q,[])),None)==truth[q]["truth_genome"] for q in truth)/len(truth)
  result["top1_species_accuracy"]=sum(
   clusters.get(next(iter(ranked.get(q,[])),None),"*")==truth[q].get("truth_cluster",truth[q].get("cluster","")) for q in truth
  )/len(truth)
  result["top5_species_accuracy"]=sum(
   any(clusters.get(hit[0],hit[0])==truth[q].get("truth_cluster",truth[q].get("cluster","")) for hit in ranked.get(q,[])[:5]) for q in truth
  )/len(truth)
  summary.append(result)
 with open(out/"baseline_summary.json","w") as fh:fh.write(json.dumps(summary,indent=2))
 print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
