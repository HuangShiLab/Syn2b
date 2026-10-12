#!/usr/bin/env python3
import csv,collections,glob,json
from pathlib import Path
root=Path('/Users/macstudio/Downloads/Syn2b/external_data')
carddir=next(d for d in (root/'card').glob('card-data-*') if d.is_dir())
aro={}
with open(carddir/'aro_index.tsv') as fh:
 for row in csv.DictReader(fh,delimiter='\t'):
  protein=row.get('Protein Accession','').strip()
  if protein:aro[protein]=row
amr={}
with open(root/'amrfinder/ReferenceGeneCatalog.txt') as fh:
 for row in csv.DictReader(fh,delimiter='\t'):
  protein=row.get('refseq_protein_accession','').strip() or row.get('genbank_protein_accession','').strip()
  if protein:amr[protein]=row
resclass={}
resdir=next(d for d in (root/'resfinder/db').glob('*') if d.is_dir())
with open(resdir/'phenotypes.txt') as fh:
 for row in csv.DictReader(fh,delimiter='\t'):
  accession=row.get('Gene_accession no.','').strip()
  if accession:resclass[accession]=row.get('Class','').strip() or 'Unclassified'
def best_hits(path):
 best={}
 for line in open(path):
  f=line.rstrip().split('\t')
  if len(f)<12:continue
  query,target,bits=f[0],f[1],float(f[11])
  if query not in best or bits>best[query][0]:best[query]=(bits,target,f[2],f[3])
 return best
prot_hits=best_hits(root/'amr_validation/card_vs_amrfinder.m8')
res_hits=best_hits(root/'amr_validation/card_vs_resfinder.m8')
rows=[]
for query,(bits,target,fident,alnlen) in sorted(prot_hits.items()):
 aro_id=next((x for x in query.split('|') if x.startswith('ARO:')),'')
 card=aro.get(aro_id,{})
 targetrow=amr.get(target,{})
 rows.append({'search':'CARD_vs_AMRFinderPlus','query_header':query,'aro_id':aro_id,
  'card_family':card.get('AMR Gene Family',''),'card_drug_class':card.get('Drug Class',''),
  'target_header':target,'target_gene_family':targetrow.get('gene_family',''),
  'target_class':targetrow.get('class',''),'target_subtype':targetrow.get('subtype',''),
  'bitscore':bits,'identity':float(fident),'alnlen':int(alnlen)})
for query,(bits,target,fident,alnlen) in sorted(res_hits.items()):
 aro_id=next((x for x in query.split('|') if x.startswith('ARO:')),'')
 card=aro.get(aro_id,{})
 targetbase=target.split()[0]
 rows.append({'search':'CARD_vs_ResFinder','query_header':query,'aro_id':aro_id,
  'card_family':card.get('AMR Gene Family',''),'card_drug_class':card.get('Drug Class',''),
  'target_header':target,'target_gene_family':'','target_class':resclass.get(targetbase,'Unclassified'),
  'target_subtype':'','bitscore':bits,'identity':float(fident),'alnlen':int(alnlen)})
out=root/'oapgc_functional/amr_crosswalk_best_hits.tsv'
fields=list(rows[0])
with open(out,'w',newline='') as fh:
 writer=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
bysearch=collections.defaultdict(collections.Counter)
for row in rows:
 bysearch[row['search']]['queries']+=1
 bysearch[row['search']]['identity_ge_90']+=float(row['identity'])>=0.90
 bysearch[row['search']]['identity_ge_80']+=float(row['identity'])>=0.80
summary=[]
for search,c in bysearch.items():
 summary.append({'search':search,'queries_with_hits':c['queries'],'identity_ge_80':c['identity_ge_80'],
                 'identity_ge_90':c['identity_ge_90'],
                 'fraction_identity_ge_80':c['identity_ge_80']/c['queries'],
                 'fraction_identity_ge_90':c['identity_ge_90']/c['queries']})
with open(root/'oapgc_functional/amr_crosswalk_summary.json','w') as fh:json.dump(summary,fh,indent=2)
print(json.dumps(summary,indent=2))
