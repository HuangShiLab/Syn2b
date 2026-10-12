#!/usr/bin/env python3
import argparse,json
ap=argparse.ArgumentParser();ap.add_argument('--antismash-dir',required=True);ap.add_argument('--output',required=True)
a=ap.parse_args();n=0
with open(a.output,'w') as out:
 for path in sorted((a.antismash_dir+'/') .__class__ and __import__('pathlib').Path(a.antismash_dir).glob('*/*.json')):
  data=json.load(open(path));genome=__import__('pathlib').Path(data.get('input_file','genome')).stem
  for record in data.get('records',[]):
   regions=[]
   for feature in record.get('features',[]):
    if feature.get('type')=='region':
     q=feature.get('qualifiers',{});regions.append((feature['location'],(q.get('region_number') or [''])[0]))
   for feature in record.get('features',[]):
    if feature.get('type')!='CDS':continue
    q=feature.get('qualifiers',{});locus=(q.get('locus_tag') or [''])[0];protein=''.join(q.get('translation',[]))
    if not protein:continue
    hit=next(((loc,num) for loc,num in regions),None)
    if not hit:continue
    region_loc,region_number=hit
    header=f"{genome}|{record['id']}|{region_number}|{locus}"
    out.write(f'>{header}\n{protein}\n');n+=1
print(n)
