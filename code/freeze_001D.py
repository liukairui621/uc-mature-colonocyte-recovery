from pathlib import Path
import json,hashlib,datetime,subprocess
P=Path('/root/projects/UC_Treatment_Recovery');O=P/'runs/001D_annotation'
def rec(p):
 return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'mtime_utc':datetime.datetime.fromtimestamp(p.stat().st_mtime,datetime.timezone.utc).isoformat()}
review=P/'provenance/reviews/001D_annotation_protocol_independent_v1.json'
d=json.loads(review.read_text());assert d['status']=='REVIEWED',d['status']
draft=P/'planning/signature_lock_v1.1.draft.json';lock=json.loads(draft.read_text())
for r in lock['amended_frozen_inputs']:
 assert rec(Path(r['path']))['sha256']==r['sha256'],r['path']
# Every primary analytical input, script, config, output and log gets a checksum.
ins=['runs/001A_qc/discovery_preprocessed_v1.rds','runs/001C_programs/program_scores.rds','inputs/annotation/program_membership_001C.tsv','provenance/reviews/metadata_independent_v1_gsm_barcode.tsv','inputs/annotation/hgnc_complete_set_20260910.txt','inputs/annotation/GPL13158.annot.gz','inputs/annotation/GPL570.annot.gz','inputs/annotation/GPL6244.relay.annot.gz','renv.lock']
paths=[P/f for f in ins]+list(O.glob('*'))+[P/f for f in ['code/annotation_001D.py','code/discovery_001D.R','code/annotation_bridge_001D.R','code/validation_002A.R','code/download_annotations_001D.py','code/amend_lock_001D.py','planning/Analysis_Amendment_001D.md']]+[review]
manifest={'run_id':'001D_annotation_amendment','status':'REVIEWED','validation_expression_accessed':False,'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'commands':['python3 code/download_annotations_001D.py','python3 code/annotation_001D.py','Rscript code/discovery_001D.R','Rscript code/annotation_bridge_001D.R','python3 code/amend_lock_001D.py'],'working_directory':str(P),'seed':20260910,'software_environment':'renv.lock; runs/001D_annotation/sessionInfo.txt','branch':'Version1.1 HGNC/probe annotation plus pre-validation common194 amendment; old001C retained','review':rec(review),'artifacts':[rec(f) for f in sorted(set(paths)) if f.is_file()]}
mf=P/'provenance/manifests/stage_001D_reviewed.json';mf.write_text(json.dumps(manifest,indent=2))
for f in ['code','planning','provenance/reviews','provenance/manifests/stage_001D_reviewed.json']:
 subprocess.run(['git','add',f],cwd=P,check=True)
subprocess.run(['git','commit','-m','Fix cross-platform gene annotation and freeze prevalidation safeguards'],cwd=P,check=True)
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=P,text=True).strip()
lock['status']='FROZEN';lock['frozen_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();lock['source_code_commit']=commit
lock['adoption_review']=rec(review);lock['amendment_manifest']=rec(mf)
(P/'signature_lock.json').write_text(json.dumps(lock,ensure_ascii=False,indent=2))
(P/'planning/signature_lock_v1.1.frozen.json').write_bytes((P/'signature_lock.json').read_bytes())
sha=rec(P/'signature_lock.json')['sha256'];(P/'signature_lock.sha256').write_text(sha+'  signature_lock.json\n')
subprocess.run(['git','add','signature_lock.json','signature_lock.sha256','planning/signature_lock_v1.1.frozen.json'],cwd=P,check=True)
subprocess.run(['git','commit','-m','Record version1.1 freeze before held-out expression access'],cwd=P,check=True)
print(json.dumps({'version':'1.1','sha256':sha,'frozen_utc':lock['frozen_utc'],'source_code_commit':commit,'validation_expression_accessed':False},indent=2))
