#!/usr/bin/env python3
"""Extract clean OAPGC FASTAs from the original Zenodo tar archives."""
from __future__ import annotations
import argparse,gzip,io,json,tarfile
from pathlib import Path
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--raw-dir",default="/Users/macstudio/Downloads/OAPGC/raw")
 ap.add_argument("--reference-tgt-dir",required=True,help="symlink dir selecting genomes")
 ap.add_argument("--output-dir",required=True)
 args=ap.parse_args();out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
 wanted={path.resolve().stem for path in Path(args.reference_tgt_dir).glob("*.tgt")}
 print(f"requested {len(wanted)} genomes",flush=True)
 remaining=set(wanted);archives=sorted(Path(args.raw_dir).glob("OAPGC_genomes_part*.tar.gz"))
 written=0
 for archive in archives:
  if not remaining:break
  print(f"scanning {archive}",flush=True)
  with tarfile.open(archive,"r:gz") as tf:
   members=[]
   for member in tf:
    if not member.isfile():continue
    name=Path(member.name).name
    if not name.endswith(".fna.gz"):continue
    genome_id=name[:-7]
    if genome_id in remaining:members.append((member,genome_id))
   for member,genome_id in members:
    destination=out/f"{genome_id}.clean.fna"
    if destination.exists():
     remaining.remove(genome_id);written+=1;continue
    raw=tf.extractfile(member)
    if raw is None:continue
    with gzip.open(io.BytesIO(raw.read()),"rt") as source,destination.open("w") as output:
     for line in source:output.write(line)
    remaining.remove(genome_id);written+=1
   print(f"  extracted cumulative={written} remaining={len(remaining)}",flush=True)
 stats={"requested":len(wanted),"written":written,"missing":sorted(remaining)}
 (out/"extract_stats.json").write_text(json.dumps(stats,indent=2))
 print(json.dumps(stats,indent=2))
 if remaining:raise SystemExit(2)
if __name__=="__main__":main()
