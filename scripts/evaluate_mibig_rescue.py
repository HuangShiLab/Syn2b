#!/usr/bin/env python3
import argparse,csv,collections
ap=argparse.ArgumentParser()
ap.add_argument('--mibig-map',required=True)
ap.add_argument('--mmseqs',required=True)
ap.add_argument('--token-hits',required=True)
ap.add_argument('--detail',required=True)
ap.add_argument('--output',required=True)
a=ap.parse_args()
fam={r['cluster_id']:r['family'] for r in csv.DictReader(open(a.mibig_map),delimiter='\t')}
def family(x):return fam.get(x,x)
scores=collections.defaultdict(lambda:collections.defaultdict(float))
for line in open(a.mmseqs):
 f=line.rstrip().split('\t')
 if len(f)<12:continue
 q=f[0].split('|',1)[0];t=f[1].split('|',1)[0]
 scores[q][t]=max(scores[q][t],float(f[11]))
protein={q:[x for x,_ in sorted(v.items(),key=lambda z:(-z[1],z[0]))] for q,v in scores.items()}
token=collections.defaultdict(list)
for r in csv.DictReader(open(a.token_hits),delimiter='\t'):
 token[r['query_cluster_id'].rsplit(':',1)[-1]].append(r['target_cluster_id'].rsplit('MIBiG:',1)[-1])
rows=list(csv.DictReader(open(a.detail),delimiter='\t'))
def metrics(method):
 hit=top1=top5=0;combined=[]
 for row in rows:
  q=row['truth_cluster'];ta=token.get(q,[]);pa=protein.get(q,[]);truth=set(family(q).split(';'))
  if method=='token':final=ta
  elif method=='protein':final=pa
  else:
   token_ok=any(set(family(x).split(';'))&truth for x in ta[:5]);final=ta if token_ok else pa
  hit+=bool(final);top1+=bool(final and family(final[0])==family(q))
  top5+=any(family(x)==family(q) for x in final[:5]);combined.append((q,final))
 n=len(rows);return n,hit/n,top1/n,top5/n,combined
out=[]
for method in ['token','protein','token+protein_rescue']:
 n,h,t1,t5,combined=metrics(method)
 result={'method':method,'queries':n,'hit_fraction':h,'broad_family_top1_accuracy':t1,'broad_family_top5_accuracy':t5}
 out.append(result);print(result)
with open(a.output,'w',newline='') as fh:
 fields=['query_cluster','truth_family','token_family','protein_family','final_family','token_top5_broad_correct','protein_top1_broad_correct','final_top5_broad_correct']
 writer=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
 for q,final in combined:
  token_hits=token.get(q,[]);protein_hits=protein.get(q,[])
  writer.writerow({'query_cluster':q,'truth_family':family(q),'token_family':family(token_hits[0]) if token_hits else '*',
   'protein_family':family(protein_hits[0]) if protein_hits else '*','final_family':family(final[0]) if final else '*',
   'token_top5_broad_correct':any(family(x)==family(q) for x in token_hits[:5]),
   'protein_top1_broad_correct':bool(protein_hits and family(protein_hits[0])==family(q)),
   'final_top5_broad_correct':any(family(x)==family(q) for x in final[:5])})
