#!/usr/bin/env python3
"""Aggregate MMseqs2 gene hits into MIBiG family scores and combine with Syn2b token scores."""
import argparse,csv,collections,math
ap=argparse.ArgumentParser()
ap.add_argument('--mmseqs',required=True);ap.add_argument('--token-hits',required=True);ap.add_argument('--detail',required=True)
ap.add_argument('--family-map',required=True);ap.add_argument('--output',required=True)
a=ap.parse_args()
fam={}
for r in csv.DictReader(open(a.family_map),delimiter='\t'):fam[r['cluster_id']]=r['family']
def F(x):return fam.get(x,x)
def P(x):return x.split('|',1)[0]
def primary(x):return F(x).split(';')[0]
def broad(x):return ';'.join(F(x).split(';')[:2])
# best gene-level bits per query protein/target cluster
gene_best=collections.defaultdict(lambda:collections.defaultdict(float))
for line in open(a.mmseqs):
 f=line.rstrip().split('\t')
 if len(f)<12:continue
 q=P(f[0]);t=P(f[1]);bits=float(f[11])
 key=(q,f[0],t,f[1])
 if bits>gene_best[q][(t,f[0],f[1])]:gene_best[q][(t,f[0],f[1])]=bits
# aggregate: sum best bits across query genes
pscores=collections.defaultdict(lambda:collections.defaultdict(float));pgenes=collections.defaultdict(lambda:collections.defaultdict(int))
for q,hits in gene_best.items():
 for (t,qp,tp),bits in hits.items():
  pscores[q][t]+=bits;pgenes[q][t]+=1
token=collections.defaultdict(list)
for r in csv.DictReader(open(a.token_hits),delimiter='\t'):
 token[r['query_cluster_id'].rsplit(':',1)[-1]].append((r['target_cluster_id'].rsplit('MIBiG:',1)[-1],float(r['cluster_score']),float(r['functional_jaccard']),float(r['adjacency_score'])))
rows=list(csv.DictReader(open(a.detail),delimiter='\t'))
methods={}
methods['token']=lambda q,t,ps,c:(1.0 if t in [x[0] for x in token.get(q,[])] else 0.0)
methods['protein_sumbits']=lambda q,t,ps,c:ps.get(t,0.0)
methods['protein_meanbits']=lambda q,t,ps,c:(ps.get(t,0.0)/max(1,pgenes[q].get(t,0)))
# normalized combination of token score and log sumbits
def combined(q,t,ps,c):
 token_score=next((score for target,score,jac,adj in token.get(q,[]) if target==t),0.0)
 return token_score+math.log1p(ps.get(t,0.0))/10.0
methods['token+protein_logbits']=combined
methods['token+protein_sumbits']=lambda q,t,ps,c:0.5*(token_score:=next((score for target,score,jac,adj in token.get(q,[]) if target==t),0.0))+0.5*min(1.0,ps.get(t,0.0)/3000.0)
results=[]
for method,score_fn in methods.items():
 top1=top5=primary1=primary5=broad1=broad5=hit=0
 for row in rows:
  q=row['truth_cluster'];targets=collections.defaultdict(float)
  for target,score,jac,adj in token.get(q,[]):targets[target]=max(targets[target],score)
  for t in pscores.get(q,{}):targets[t]=max(targets[t],score_fn(q,t,pscores[q],None))
  ranked=sorted(targets,key=lambda t:(-targets[t],t));top_matches=ranked[:5]
  hit+=bool(top_matches);top1+=bool(top_matches and F(top_matches[0])==F(q));top5+=any(F(x)==F(q) for x in top_matches);primary1+=bool(top_matches and primary(F(top_matches[0]))==primary(F(q)));primary5+=any(primary(x)==primary(F(q)) for x in top_matches);broad1+=bool(top_matches and broad(F(top_matches[0]))==broad(F(q)));broad5+=any(broad(x)==broad(F(q)) for x in top_matches)
 n=len(rows);results.append({'method':method,'queries':n,'hit_fraction':hit/n,'exact_top1':top1/n,'exact_top5':top5/n,'primary_top1':primary1/n,'primary_top5':primary5/n,'broad_top1':broad1/n,'broad_top5':broad5/n})
with open(a.output,'w',newline='') as fh:
 fields=['method','queries','hit_fraction','exact_top1','exact_top5','primary_top1','primary_top5','broad_top1','broad_top5'];w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(results)
for row in results:print(row)
