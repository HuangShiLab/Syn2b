#!/usr/bin/env python3
"""Import MIBiG, CARD, AMRFinderPlus, and ResFinder into feature tables.

These are reference/locus imports, not gene predictions. Coordinates are used
when provided; reference-set imports that do not represent contiguous loci use
stable synthetic coordinates and explicit `*_reference` contig ids.
"""
from __future__ import annotations
import argparse,csv,glob,json,re,sys
from collections import defaultdict
from pathlib import Path
FIELDS=["gene_id","genome_id","contig_id","start","end","strand","product","ko","pfam","cog","eggnog","smash_tool","cluster_id"]
def slug(value):return re.sub(r"[^A-Za-z0-9_.-]+","_",value.strip()).strip("_")
def writer(path):
 fh=open(path,"w",newline="");w=csv.DictWriter(fh,fieldnames=FIELDS,delimiter="\t",lineterminator="\n");w.writeheader();return fh,w
def row(gene_id,genome,contig,start,end,strand,product,cluster,smash="",ko="",pfam="",cog="",eggnog=""):
 return {"gene_id":gene_id,"genome_id":genome,"contig_id":contig,"start":start,"end":end,
         "strand":strand,"product":product,"ko":ko,"pfam":pfam,"cog":cog,"eggnog":eggnog,
         "smash_tool":smash,"cluster_id":cluster}
def fasta(path):
 name=None;seq=[]
 with open(path) as fh:
  for line in fh:
   line=line.rstrip("\n\r")
   if line.startswith(">"):
    if name is not None:yield name,"".join(seq)
    name=line[1:].strip();seq=[]
   else:seq.append(line.strip())
 if name is not None:yield name,"".join(seq)
def import_mibig(root,out):
 json_dir=root/"mibig_json_4.0/mibig_json_4.0";classes={}
 for path in json_dir.glob("*.json"):
  data=json.load(open(path));accession=data.get("accession");terms=[]
  for item in data.get("biosynthesis",{}).get("classes",[]):
   for key in ["class","subclass","ripp_type"]:
    if item.get(key):terms.append(str(item[key]))
  if accession:classes[accession]=";".join(dict.fromkeys(terms))
 fasta_path=root/"mibig_prot_seqs_4.0.fasta";features=[];seen=set()
 rx=re.compile(r"^(?P<accession>[^|]+)\|(?P<contig>[^|]+)\|(?P<start>\d+)-(?P<end>\d+)\|(?P<strand>[-+])\|(?P<protein>[^|]*)(?:\|(?P<product>.*))?$")
 for header,sequence in fasta(fasta_path):
  m=rx.match(header)
  if not m:continue
  accession=m.group("accession");contig=m.group("contig");start=int(m.group("start"));end=int(m.group("end"));strand=m.group("strand")
  protein=m.group("protein").strip();product=(m.group("product") or protein or "MIBiG_gene").replace("_"," ").strip()
  gene_id=protein or f"{accession}:{contig}:{start}-{end}"
  key=(accession,contig,gene_id,start,end)
  if key in seen:continue
  seen.add(key)
  features.append(row(gene_id,"MIBiG",f"{accession}|{contig}",start,end,strand,product,accession,
                      smash=classes.get(accession,""),))
 fh,w=writer(out)
 for feature in sorted(features,key=lambda r:(r["cluster_id"],r["contig_id"],r["start"],r["end"],r["gene_id"])):
  w.writerow(feature)
 fh.close()
 with open(out.parent/"mibig_family_map.tsv","w",newline="") as mapfh:
  mw=csv.writer(mapfh,delimiter="\t",lineterminator="\n");mw.writerow(["cluster_id","family"])
  for accession in sorted(classes):mw.writerow([accession,classes[accession]])
 return len(features)
def import_card(root,out):
 paths=list(root.glob("card-data-*/aro_index.tsv"))
 if not paths:raise RuntimeError("CARD aro_index.tsv not found")
 groups=defaultdict(list)
 with open(paths[0]) as fh:
  for record in csv.DictReader(fh,delimiter="\t"):
   family=(record.get("AMR Gene Family") or record.get("Drug Class") or "Unclassified").split(";")[0].strip() or "Unclassified"
   gene_id=record.get("Protein Accession") or record.get("Model Sequence ID") or record.get("ARO Accession")
   product=record.get("ARO Name") or record.get("Model Name") or gene_id
   groups[family].append((gene_id,product,record.get("Drug Class",""),record.get("Resistance Mechanism",""),record.get("ARO Accession","")))
 features=[]
 for family in sorted(groups):
  contig=f"CARD_reference:{slug(family)}"
  for index,(gene_id,product,drug,mechanism,aro) in enumerate(sorted(groups[family])):
   start=1000*(index+1);features.append(row(gene_id,"CARD",contig,start,start+900,"+",
      cluster=f"CARD:{slug(family)}",product=f"CARD:{slug(family)}",smash=f"CARD:{slug(family)};{slug(mechanism)}",ko=family))
 fh,w=writer(out)
 for feature in features:w.writerow(feature)
 fh.close()
 return len(features),len(groups)
def import_amrfinder(path,out):
 features=[]
 with open(path) as fh:
  for record in csv.DictReader(fh,delimiter="\t"):
   allele=record.get("allele","").strip();family=record.get("gene_family","").strip()
   product=record.get("product_name","").strip() or family or allele
   gene_id=record.get("refseq_protein_accession","").strip() or record.get("genbank_protein_accession","").strip() or allele or slug(product)
   if not family:continue
   class_name=record.get("class","").strip() or "Unclassified"
   cluster=f"AMRFINDERPLUS:{slug(class_name)}"
   try:start=int(record.get("genbank_start") or record.get("refseq_start") or 0);end=int(record.get("genbank_stop") or record.get("refseq_stop") or 0)
   except ValueError:start=end=0
   if start<=0 or end<start:start=1;end=max(1,len(product))
   strand=record.get("genbank_strand") or record.get("refseq_strand") or "+"
   features.append(row(gene_id or slug(product),"AMRFinderPlus","ReferenceGeneCatalog",start,max(end,start),strand,
      product,cluster,ko=family,smash=f"AMRFINDERPLUS:{slug(class_name)}"))
 fh,w=writer(out)
 for feature in sorted(features,key=lambda r:(r["cluster_id"],r["gene_id"])):w.writerow(feature)
 fh.close()
 return len(features)
def import_resfinder(path,out):
 root=next((p for p in path.glob("*/all.fsa")),None)
 if root is None:raise RuntimeError("ResFinder all.fsa not found")
 phenotypes={}
 phenotype_path=root.parent/"phenotypes.txt"
 with open(phenotype_path) as fh:
  for record in csv.DictReader(fh,delimiter="\t"):
   accession=record.get("Gene_accession no.","").strip()
   if accession:phenotypes[accession]=record.get("Class","").strip() or "Unclassified"
 features=[]
 for header,sequence in fasta(root):
  gene_id=header.split()[0]
  base=gene_id.rsplit("_",2)[0] if gene_id.count("_")>=2 else gene_id
  product=base.replace("_"," ")
  class_name=phenotypes.get(gene_id,"Unclassified")
  cluster=f"RESFINDER:{slug(class_name)}"
  features.append(row(gene_id,"ResFinder",f"ResFinder_reference:{slug(class_name)}",1,max(1,len(sequence)),"+",
      product,cluster,ko=base,smash=f"RESFINDER:{slug(class_name)}"))
 fh,w=writer(out)
 for feature in sorted(features,key=lambda r:(r["cluster_id"],r["gene_id"])):w.writerow(feature)
 fh.close()
 return len(features)
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--external-root",default="external_data");ap.add_argument("--output-dir",default="external_data/imported")
 args=ap.parse_args();root=Path(args.external_root);out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
 mibig=import_mibig(root/"mibig",out/"mibig_features.tsv")
 card,card_clusters=import_card(root/"card",out/"card_features.tsv")
 amr=import_amrfinder(root/"amrfinder/ReferenceGeneCatalog.txt",out/"amrfinder_features.tsv")
 res=import_resfinder(root/"resfinder/db",out/"resfinder_features.tsv")
 stats={"mibig_genes":mibig,"card_genes":card,"card_cluster_families":card_clusters,
        "amrfinder_records":amr,"resfinder_records":res}
 with open(out/"import_stats.json","w") as fh:json.dump(stats,fh,indent=2)
 print(json.dumps(stats,indent=2))
if __name__=="__main__":main()
