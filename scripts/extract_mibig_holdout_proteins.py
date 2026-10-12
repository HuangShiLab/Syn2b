#!/usr/bin/env python3
import csv
from pathlib import Path
root=Path('/Users/macstudio/Downloads/Syn2b')
detail=root/'validation/mibig_family_holdout_versioned/detail.tsv'
rows=list(csv.DictReader(open(detail),delimiter='\t'))
queries=sorted({r['truth_cluster'] for r in rows})
out=root/'external_data/mibig/functional/heldout_proteins.faa'
target=root/'external_data/mibig/functional/target_proteins.faa'
qfh=open(out,'w')
tfh=open(target,'w')
nq=nt=0
name=None
seq=[]
def flush():
 global name,seq
 if name is not None and seq:
  cluster_version=name.split('|',1)[0]
  sequence=''.join(seq).replace('*','')
  handle=qfh if cluster_version in queries else tfh
  handle.write('>'+name+'\n'+sequence+'\n')
 name=None;seq=[]
with open(root/'external_data/mibig/functional/mibig_proteins.nostop.faa') as fh:
 for line in fh:
  line=line.rstrip()
  if line.startswith('>'):
   flush()
   name=line[1:].strip();seq=[]
  else:
   seq.append(line)
flush()
qfh.close();tfh.close()
print('query clusters',len(queries),'query proteins',sum(line.startswith('>') for line in open(out)),'target proteins',sum(line.startswith('>') for line in open(target)))
