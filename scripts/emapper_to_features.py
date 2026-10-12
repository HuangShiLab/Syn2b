#!/usr/bin/env python3
"""Convert emapper annotations + Prodigal GFFs into Syn2b feature tables."""
from __future__ import annotations
import argparse,csv,gzip
from pathlib import Path

FIELDS=["gene_id","genome_id","contig_id","start","end","strand","product","ko","pfam","cog","eggnog","smash_tool","cluster_id"]

def attributes(text):
    out={}
    for item in text.rstrip().rstrip(";").split(";"):
        if "=" in item:
            key,value=item.split("=",1)
            out[key]=value
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--annotation-dir",required=True)
    parser.add_argument("--emapper",required=True)
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    root=Path(args.annotation_dir)
    coords={}
    for path in root.rglob("*.gff.gz"):
        genome=path.name.removesuffix(".gff.gz")
        with gzip.open(path,"rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields=line.rstrip().split("\t")
                if len(fields)<9 or fields[2]!="CDS":
                    continue
                attr=attributes(fields[8])
                protein=attr.get("ID")
                if protein:
                    coords.setdefault(genome,{}).setdefault(fields[0],{})[protein]=(fields[0],int(fields[3]),int(fields[4]),fields[6])
    features=[]
    with open(args.emapper) as handle:
        columns=None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                columns=line[1:].rstrip("\n\r").split("\t")
                continue
            if line.startswith("#"):
                continue
            values=line.rstrip("\n\r").split("\t")
            row=dict(zip(columns,values))
            query=row["query"]
            if query.count("_")<2:
                continue
            genome,remainder=query.split("_",1)
            contig_parts=remainder.split("_")
            contig=genome+"_"+contig_parts[0]
            protein=remainder
            coordinate=coords.get(genome,{}).get(contig,{}).get(protein)
            if coordinate is None:
                continue
            coord_contig,start,end,strand=coordinate
            ko=";".join(item for item in row.get("KEGG_ko","").split(",") if item)
            pfam=";".join(item for item in row.get("PFAMs","").split(",") if item)
            features.append({"gene_id":query,"genome_id":genome,"contig_id":coord_contig,
              "start":start,"end":end,"strand":strand,"product":row.get("Description") or "",
              "ko":ko,"pfam":pfam,"cog":row.get("COG_category",""),
              "eggnog":row.get("eggNOG_OGs",""),"smash_tool":"","cluster_id":""})
    with open(args.output,"w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=FIELDS,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        for feature in sorted(features,key=lambda item:(item["genome_id"],item["contig_id"],item["start"],item["end"],item["gene_id"])):
            writer.writerow(feature)
    print(f"wrote {len(features)} features to {args.output}")

if __name__=="__main__":
    main()
