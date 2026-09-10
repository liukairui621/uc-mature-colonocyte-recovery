from pathlib import Path
import json, hashlib, datetime, subprocess, shutil
p=Path(__file__).resolve().parents[1]
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''): h.update(chunk)
 return h.hexdigest()
def rec(q): return {'path':str(q),'bytes':q.stat().st_size,'sha256':digest(q)}
def write(q,d): q.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
reviews=[p/'provenance/reviews'/f for f in [
 '002B_precision_independent_v1.json','002B_directional_power_independent_v1.json',
 '002B_GSE73661_IFX_metadata_independent_v1.json','002B_GSE73661_results_independent_v1.json',
 '002B_GSE282122_asset_audit_v1.json','002B_report_independent_v1.json','002B_visual_QA.json']]
for q in reviews: assert json.loads(q.read_text())['status']=='REVIEWED',str(q)
manifest=p/'provenance/manifests/stage_002B_completion.json'
assert not manifest.exists(),'Completion manifest already exists; do not overwrite.'
integrity=json.loads((p/'provenance/reviews/002B_prior_and_freeze_integrity.json').read_text())
assert integrity['status']=='PASS'
for item in integrity['checks']: assert item['match']
assert digest(p/'signature_lock.json')=='62462dfbc25bb857c94c872424173e215c6d1b8b6724a0bbc724a5c2683bfb36'
assert digest(p/'planning/support_plan_002B.frozen.json')=='d8e5549c0e9e0a60696488653b45933fda7e8470da174351ab4377fe7de43399'
artifactchecks=[]
def check_records(v):
 if isinstance(v,dict):
  if isinstance(v.get('path'),str) and isinstance(v.get('sha256'),str):
   q=Path(v['path'])
   if not q.is_absolute(): q=p/q
   if str(q).startswith(str(p)+'/'):
    assert q.is_file(),str(q)
    match=digest(q)==v['sha256']
    artifactchecks.append({'path':str(q),'match':match})
    assert match,str(q)
  for v1 in v.values():check_records(v1)
 elif isinstance(v,list):
  for v1 in v:check_records(v1)
for q in reviews: check_records(json.loads(q.read_text()))
reviewcheck=p/'provenance/reviews/002B_final_review_hash_integrity.json'
assert not reviewcheck.exists()
write(reviewcheck,{'status':'PASS','checked_utc':now,'checks':artifactchecks})
report=p/'reports/Precision_review_and_002B_support_20260910.md'
oldreport=p/'reports/Annotation_amendment_and_GSE23597_validation_20260910.md'
sup=p/'provenance/manifests/report_supersession_002B.json'
assert not sup.exists()
write(sup,{'updated_utc':now,'current_stage_report':rec(report),'prior_report_retained':rec(oldreport),'relationship':'Supplement with precision and preassigned supportive endpoint results. Does not reverse or replace the GSE23597 primary finding.'})
state=json.loads((p/'project_status.json').read_text())
state.update(current_phase='002B COMPLETE AND REVIEWED: precision audit, GSE73661 supportive endoscopic analysis, GSE282122 metadata/feature asset audit',
 next_phase='003 GSE282122 cell-source, composition and within-state exploration; expression X not yet accessed',
 updated_utc=now,review_pending=False,stage_report=str(report),
 independent_review=str(p/'provenance/reviews/002B_report_independent_v1.json'),
 independent_result_review=str(p/'provenance/reviews/002B_GSE73661_results_independent_v1.json'),
 completion_ledger=str(manifest),primary_validation_support=False,
 GSE73661_expression_accessed=True,GSE73661_analysis_completed=True,
 GSE73661_evidence_role='descriptive_supportive_only; cannot replace GSE23597',
 GSE73661_endoscopic_positive_association=True,
 GSE282122_metadata_and_features_accessed=True,GSE282122_expression_X_accessed=False,
 GSE282122_cell_analyses_completed=False)
src=str(p/'inputs/support/GSE73661_series_matrix.txt.gz')
state['analysis_input_files']=list(dict.fromkeys(state.get('analysis_input_files',[])+[src]))
write(p/'project_status.json',state)
snapshot=p/'planning/project_status_at_002B_completion.json'
assert not snapshot.exists()
write(snapshot,state)
inputs=set()
for f in [
 'signature_lock.json','planning/support_plan_002B.frozen.json','renv.lock',
 'inputs/support/GSE73661_series_matrix.txt.gz','inputs/metadata/GSE73661.matrix_header.txt',
 'inputs/annotation/hgnc_complete_set_20260910.txt',
 'runs/001D_annotation/discovery_preprocessed_canonical.rds','runs/001D_annotation/program_scores.rds',
 'runs/001D_annotation/discovery_model_fits.rds','runs/001D_annotation/GPL6244_canonical_probe_mapping.tsv',
 'runs/001D_annotation/program_membership_canonical.tsv','runs/001D_annotation/frozen_common_inflammation_members.tsv',
 'runs/002A_GSE23597/validation_scores.rds','runs/002A_GSE23597/candidate_model_fits.rds',
 'runs/002A_GSE23597/primary_patient_data.tsv','runs/002A_GSE23597/validation_summary.json',
 'provenance/manifests/stage_001D_reviewed.json','provenance/manifests/stage_002A_completion.json']:
 inputs.add(p/f)
for folder in ['inputs/reviews_002B','inputs/002B_methods','inputs/002B_asset_audit']:
 inputs.update(q for q in (p/folder).rglob('*') if q.is_file())
outputs=set([report,snapshot,sup,p/'planning/Next_stage_003_cell_context.md'])
for folder in ['runs/002B_precision','runs/002B_GSE73661_support','runs/002B_figures']:
 outputs.update(q for q in (p/folder).rglob('*') if q.is_file())
code=set((p/'code').glob('*002B*'))
code.update([p/'code/figure_theme.R'])
evidence=set((p/'provenance/reviews').glob('002B*'))
logs=set((p/'logs').glob('*002B*'))
head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=p,text=True).strip()
d={'run_id':'002B','status':'REVIEWED','completed_utc':now,
 'scope':'Conditional precision audit; internally prespecified descriptive GSE73661 endoscopic support; GSE282122 metadata and feature audit only',
 'current_source_commit':head,'support_execution_source_commit':'8153b329a5ddddbe4f47f8ff487fb29878478702',
 'support_freeze_commit':'67e765d07cd93227b5f454db455d7de4af1c8aa6',
 'commands':['Rscript code/precision_002B.R','Rscript code/design_73661_002B.R','Rscript code/directional_power_002B.R','Rscript code/support_73661_002B.R','Rscript code/plots_002B.R','python3 code/report_002B.py'],
 'seed':20260910,'cwd':str(p),
 'primary_validation_support_remains_false':True,'additional_confirmatory_primary_tests':0,
 'GSE282122_expression_X_accessed':False,
 'precision_summary':json.loads((p/'runs/002B_precision/precision_summary.json').read_text()),
 'GSE73661_summary':json.loads((p/'runs/002B_GSE73661_support/support_summary.json').read_text()),
 'reviews':[rec(q) for q in reviews],
 'inputs':[rec(q) for q in sorted(inputs)],
 'outputs':[rec(q) for q in sorted(outputs)],
 'code':[rec(q) for q in sorted(code)],
 'review_evidence':[rec(q) for q in sorted(evidence)],
 'logs':[rec(q) for q in sorted(logs)],
 'mutable_project_status_policy':'Completion snapshot is hashed; mutable root project_status.json is not hashed so later stages can update it.',
 'next_phase_status':'PLANNED; no cell-level expression or composition results yet'}
write(manifest,d)
print(json.dumps({'manifest':str(manifest),'sha256':digest(manifest),'inputs':len(inputs),'outputs':len(outputs),'reviews':len(reviews),'checked_review_artifact_records':len(artifactchecks)},ensure_ascii=False))
