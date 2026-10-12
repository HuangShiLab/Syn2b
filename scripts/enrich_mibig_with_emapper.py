#!/usr/bin/env python3
import argparse,csv
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--features',required=True);ap.add_argument('--emapper',required=True);ap.add_argument('--output',required=True)
 a=ap.parse_args()
 ann={}
 with open(a.emapper) as fh:
  columns=None
  for line in fh:
   if line.startswith('##'):continue
   if line.startswith('#query'):
    columns=line[1:].rstrip('\n\r').split('\t');continue
   if line.startswith('#'):continue
   row=dict(zip(columns,line.rstrip('\n\r').split('\t')))
   query=row['query'];protein=query.split('|')[4] if query.count('|')>=4 else query
   ann[protein]={
    'product':row.get('Description') or '',
    'ko':';'.join(x.replace('ko:','') for x in row.get('KEGG_ko','').split(',') if x and x!='-'),
    'pfam':';'.join(x for x in row.get('PFAMs','').split(',') if x and x!='-'),
    'eggnog':';'.join(x.split('|')[0] for x in row.get('eggNOG_OGs','').split(',') if x and x!='-'),
   }
 rows=[];matched=0
 with open(a.features) as fh:
  reader=csv.DictReader(fh,delimiter='\t');fields=reader.fieldnames
  for row in reader:
   old=None;e=ann.get(row['gene_id'])
   if e:
    matched+=1;old={k:row[k] for k in ['product','ko','pfam','eggnog']}
    for field in ['product','ko','pfam','eggnog']:
     if e[field]:row[field]=e[field]
   rows.append(row)
 with open(a.output,'w',newline='') as out:
  writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
 print(f'matched {matched}/{len(rows)}')
