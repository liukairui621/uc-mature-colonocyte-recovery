from pathlib import Path
import json,hashlib,datetime,subprocess,os
P=Path('/root/projects/UC_Treatment_Recovery');V=P/'runs/002A_GSE23597';R=P/'provenance/reviews'
def rec(p):
 return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
reviews=[R/'002A_validation_results_independent_v1.json',R/'002A_report_independent_v1.json',R/'002A_visual_QA.json']
for p in reviews:
 d=json.loads(p.read_text());assert d['status']=='REVIEWED',(p,d['status'])
lock=json.loads((P/'signature_lock.json').read_text())
checked=[]
for e in lock['amended_frozen_inputs']:
 actual=rec(Path(e['path']));assert actual['sha256']==e['sha256'],e['path'];checked.append(actual)
old=json.loads((P/'provenance/manifests/stage_001D_reviewed.json').read_text())
for e in old['artifacts']:
 assert rec(Path(e['path']))['sha256']==e['sha256'],e['path']
integrity={'status':'PASS','checked_v1_1_frozen_inputs':len(checked),'checked_001D_manifest_artifacts':len(old['artifacts']),'lock_sha256':rec(P/'signature_lock.json')['sha256'],'old_lock_v1_0_sha256':rec(P/'planning/signature_lock_v1.0.json')['sha256'],'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'002A_final_hash_integrity.json').write_text(json.dumps(integrity,indent=2))
report=P/'reports/Annotation_amendment_and_GSE23597_validation_20260910.md'
info=json.loads((V/'validation_summary.json').read_text())
manifest={'run_id':'002A_GSE23597','status':'REVIEWED','scope':'One frozen primary validation; no alternative primary','command':'Rscript code/validation_002A.R','seed':20260910,'cwd':str(P),'source_code_commit':lock['source_code_commit'],'lock':rec(P/'signature_lock.json'),'expression_access_event':rec(V/'expression_access_event.json'),'result_summary':info,'reviews':[rec(f) for f in reviews],'inputs':checked+[rec(P/'inputs/validation/GSE23597_series_matrix.txt.gz')],'outputs':[rec(f) for f in sorted(V.rglob('*')) if f.is_file()],'review_evidence':[rec(f) for f in sorted(R.glob('002A*')) if f.is_file()],'code':[rec(P/f) for f in ['code/validation_002A.R','code/plots_002A.R','code/report_002A.py','code/finalize_002A.py','code/figure_theme.R']],'report':rec(report),'environment':rec(P/'renv.lock'),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
ledger=P/'provenance/manifests/stage_002A_completion.json';ledger.write_text(json.dumps(manifest,indent=2))
status=P/'project_status.json';oldstatus=json.loads(status.read_text())
archive=P/'planning/project_status_before_002A.json'
if not archive.exists():archive.write_bytes(status.read_bytes())
oldstatus.update(current_phase='002A GSE23597 completed and independently reviewed; prespecified primary support criterion not met',candidate_frozen=True,signature_lock_version='1.1',independent_validation_started=True,independent_validation_completed=True,primary_validation_evaluable=True,primary_validation_support=False,validation_expression_accessed=True,next_phase='002B preassigned supportive/contextual work only; GSE73661 cannot replace primary; GSE282122 cell-source interpretation pending',review_pending=False,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage_report=str(report),independent_review=str(reviews[0]),completion_ledger=str(ledger))
oldstatus['analysis_input_files']=list(dict.fromkeys(oldstatus.get('analysis_input_files',[])+[str(P/'inputs/validation/GSE23597_series_matrix.txt.gz')]))
status.write_text(json.dumps(oldstatus,indent=2))
sup={'status':'REVIEWED','current_report':rec(report),'previous_report':'reports/Review_response_and_001C_report_20260910.md','old_artifacts':'Retained unchanged; current program estimates use001D canonical annotation and common194. Partial IBrD proxy retired. No old GO/CAMERA run has been relabelled as new-annotation output.','primary_validation':'GSE23597 technically evaluable, support rule false. GSE73661 never substitutes.'}
(P/'provenance/manifests/report_supersession_002A.json').write_text(json.dumps(sup,indent=2))
env=os.environ.copy();env.update(GIT_AUTHOR_NAME='Codex Research Agent',GIT_AUTHOR_EMAIL='codex@localhost',GIT_COMMITTER_NAME='Codex Research Agent',GIT_COMMITTER_EMAIL='codex@localhost')
subprocess.run(['git','add','code','reports','planning','provenance','project_status.json'],cwd=P,check=True)
subprocess.run(['git','commit','-m','Report frozen GSE23597 validation and independent numerical audit'],cwd=P,check=True,env=env)
print(json.dumps({'status':'REVIEWED','primary_support':False,'report':rec(report),'ledger':rec(ledger),'integrity':integrity,'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=P,text=True).strip(),'git_status':subprocess.check_output(['git','status','--short'],cwd=P,text=True)},indent=2))
