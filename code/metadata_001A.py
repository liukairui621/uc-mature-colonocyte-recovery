from pathlib import Path
import csv, io, json, hashlib, datetime, re
from collections import Counter, defaultdict
ROOT=Path('/root/projects/UC_Treatment_Recovery')
OUT=ROOT/'runs/001A_metadata'
def parse(p):
    rows=[r for r in csv.reader(io.StringIO(p.read_text()),delimiter='\t') if r and r[0].startswith('!Sample_')]
    ids=next(r[1:] for r in rows if r[0]=='!Sample_geo_accession')
    assert len(ids)==len(set(ids))
    samples=[{'gsm':x} for x in ids]
    for row in rows:
        assert len(row)==len(ids)+1,(p.name,row[0],len(row),len(ids))
        key=row[0][8:]
        for s,v in zip(samples,row[1:]):
            if key.startswith('characteristics') and ': ' in v:
                k,val=v.split(': ',1)
                if k in s and s[k]!=val:raise ValueError((s['gsm'],k))
                s[k]=val
            elif key in ['title','source_name_ch1','platform_id','supplementary_file','data_processing']:
                s.setdefault(key,[]);s[key].append(v)
    return samples
def one(s,k):return s.get(k,[''])[0] if isinstance(s.get(k),list) else s.get(k,'')
def write_tsv(name,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(rows)
registry=[];raw={};issues=[]
for gs in ['GSE92415','GSE23597','GSE73661']:
    samples=parse(ROOT/'inputs/metadata'/f'{gs}.matrix_header.txt');raw[gs]=samples
    for s in samples:
        row=dict(series=gs,gsm=s['gsm'],title=one(s,'title'),platform=one(s,'platform_id'),role={'GSE92415':'discovery','GSE23597':'held_out_bulk','GSE73661':'held_out_context'}[gs])
        if gs=='GSE92415':
            row.update(subject=s.get('subject',''),disease=s['disease'],tissue=s['tissue'],treatment=s['treatment'],time=s['visit'],week={'Week 0':'0','Week 6':'6','Healthy':''}[s['visit']],response=s['wk6response'],endpoint='clinical_response_week6',mayo_total=s['mayo score'],age=s.get('age',''))
        elif gs=='GSE23597':
            row.update(subject=s['subject'],disease='UC',tissue=one(s,'source_name_ch1'),treatment='placebo' if s['dose']=='placebo' else 'IFX',dose=s['dose'],time=s['time'],week=s['time'].replace('W',''),response=s['wk8 response'],response_w30=s['wk30 response'],endpoint='clinical_response_week8_and30')
        else:
            arm=s['induction therapy_maintenance therapy'];week=s['week (w)']
            row.update(subject=s['study individual number'],disease='non_IBD_control' if arm=='CO' else 'UC',tissue=s['tissue'],treatment=arm,time=week,week=week.replace('W','') if '_' not in week else '4-6',mayo_endoscopic=s['mayo endoscopic subscore'],endpoint='mayo_endoscopic_at_visit',exposure_status='requires_protocol_timeline_reconstruction')
        row['subject_linkage']='source_ID' if row['subject'] else 'not_provided'
        registry.append(row)
write_tsv('sample_registry_v1.tsv',registry)
(ROOT/'runs/001A_metadata/raw_sample_fields.json').write_text(json.dumps(raw,indent=2))
d=[r for r in registry if r['series']=='GSE92415']; byid=defaultdict(list)
for s in d:
    if s['subject']:byid[s['subject']].append(s)
pairs=[];unpaired=[]
for subject,ss in sorted(byid.items()):
    bytime=defaultdict(list)
    for s in ss:bytime[s['time']].append(s)
    if len(bytime['Week 0'])==1 and len(bytime['Week 6'])==1:
        a,b=bytime['Week 0'][0],bytime['Week 6'][0]
        for k in ['treatment','response','age','disease','tissue']:
            assert a[k]==b[k],(subject,k,a[k],b[k])
        pairs.append(dict(subject=subject,arm=b['treatment'],response=b['response'],baseline_gsm=a['gsm'],post_gsm=b['gsm'],mayo_baseline=a['mayo_total'],mayo_post=b['mayo_total'],age=a['age']))
    else:
        unpaired.extend(ss)
write_tsv('discovery_pairs_v1.tsv',pairs);write_tsv('discovery_unpaired_v1.tsv',unpaired)
summary={'run_id':'001A_metadata','phase':'technical_asset_audit','status':'COMPLETE','not_final_sample_counts':True,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'arrays_by_series':dict(Counter(r['series'] for r in registry)),'discovery_pairs':dict(Counter(r['arm'] for r in pairs)),'paired_response':dict(Counter(r['arm']+'_'+r['response'] for r in pairs)),'discovery_unpaired':dict(Counter(r['time'] for r in unpaired)),'unresolved':['GSE92415 overview count versus unpaired subject IDs; do not infer matching','GSE23597 P13 repeated baseline; hold duplicate decision until source check','GSE73661 repeated control subject 7; actual induction/maintenance exposure reconstruction','Clinical response is distinct from mucosal or histological healing'],'validation_expression_accessed':False}
(OUT/'metadata_summary.json').write_text(json.dumps(summary,indent=2))
manifest={'run_id':'001A_metadata','inputs':[],'outputs':[],'code':str(ROOT/'code/metadata_001A.py'),'command':'python3 code/metadata_001A.py','status':'COMPLETE'}
for key,files in [('inputs',list((ROOT/'inputs/metadata').glob('*.matrix_header.txt'))),('outputs',list(OUT.glob('*')))]:
    for f in files:manifest[key].append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
manifest['code_sha256']=hashlib.sha256((ROOT/'code/metadata_001A.py').read_bytes()).hexdigest()
(ROOT/'provenance/manifests/metadata_001A.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(summary,indent=2))
