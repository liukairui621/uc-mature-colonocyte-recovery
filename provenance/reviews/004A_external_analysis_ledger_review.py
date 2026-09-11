from pathlib import Path
import json,csv,hashlib,math,collections,datetime,platform
p=Path('/root/projects/UC_Treatment_Recovery')
out=p/'provenance/reviews/004A_external_analysis_ledger_review.json'
if out.exists():raise SystemExit('Refuse overwrite')
O=p/'runs/004A_evidence_ledger'
def table(rel):return list(csv.DictReader((p/rel).open(),delimiter='\t'))
def one(rows,predicate=lambda r:True):
 selected=[r for r in rows if predicate(r)]
 assert len(selected)==1
 return selected[0]
def fm(rel):
 f=p/rel
 return {'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size}
overview=table('runs/004A_evidence_ledger/external_analysis_overview.tsv')
a=one(table('runs/002A_GSE23597/primary_test.tsv'))
b=one(table('runs/002B_GSE73661_support/candidate_support_models.tsv'),lambda r:r['model']=='adjusted_support' and r['method']=='OLS')
vd=table('runs/002D_VDZ_support/candidate_models.tsv')
w6=one(vd,lambda r:r['visit']=='W6' and r['model']=='adjusted_support' and r['method']=='OLS')
w12=one(vd,lambda r:r['visit']=='W12' and r['model']=='adjusted_support' and r['method']=='OLS')
sc=table('runs/003_cell_context/patient_group_models.tsv')
ct=one(sc,lambda r:r['branch']=='exclude_predicted_doublets::all_fixed::threshold20' and r['metric']=='delta_ct_state_score')
orig={'002A':(a,'p'),'002B':(b,'p'),'002D_W6':(w6,'p'),'002D_W12':(w12,'p'),'003_CT':(ct,'p_exact_permutation')}
assert set(r['analysis_id'] for r in overview)==set(orig)
for row in overview:
 source,pc=orig[row['analysis_id']]
 for col in ['n','estimate','lower','upper']:
  assert math.isclose(float(row[col]),float(source[col]),rel_tol=1e-12,abs_tol=1e-14)
 assert math.isclose(float(row['p']),float(source[pc]),rel_tol=1e-12,abs_tol=1e-14)
ap=table('runs/002A_GSE23597/primary_patient_data.tsv')
bp=table('runs/002B_GSE73661_support/support_patient_data.tsv')
vp=table('runs/002D_VDZ_support/patient_level_scores.tsv')
cp=table('runs/003_cell_context/joint_primary_patient_metrics.tsv')
counts={'002A':(len(ap),sum(r['response']=='Yes' for r in ap)),
'002B':(len(bp),sum(r['endoscopic_healing']=='Yes' for r in bp)),
'002D_W6':(sum(r['visit']=='W6' for r in vp),sum(r['visit']=='W6' and r['endoscopic_healing']=='Yes' for r in vp)),
'002D_W12':(sum(r['visit']=='W12' for r in vp),sum(r['visit']=='W12' and r['endoscopic_healing']=='Yes' for r in vp)),
'003_CT':(len(cp),sum(r['remission']=='Remission' for r in cp))}
for r in overview:
 n,pos=counts[r['analysis_id']]
 assert int(r['n'])==n and int(r['n_outcome_positive'])==pos and int(r['n_outcome_negative'])==n-pos
raw=table('runs/003_cell_context/fixed_epithelial_state_models.tsv')
dc=table('runs/003D_decontx_state_family/decontx_fixed15_state_models.tsv')
rawct=one(raw,lambda r:r['branch']=='Non ileal CT colonocyte')
dcct=one(dc,lambda r:r['final_analysis']=='Non ileal CT colonocyte')
line=one(overview,lambda r:r['analysis_id']=='003_CT')
for k,v in [('BH_two_axes',ct['BH_two_main']),('BH_raw_fixed15',rawct['BH_fixed_state_family_n15']),('BH_DecontX_fixed15',dcct['BH_fixed_state_family_n15'])]:
 assert math.isclose(float(line[k]),float(v),rel_tol=1e-12,abs_tol=1e-14)
assert len(raw)==len(dc)==15
assert sum(r['p_exact_permutation']!='NA' for r in raw)==sum(r['p_exact_permutation']!='NA' for r in dc)==10
lock=json.loads((p/'signature_lock.json').read_text())
pb=json.loads((p/'planning/support_plan_002B.frozen.json').read_text())
pd=json.loads((p/'planning/support_plan_002D_VDZ.frozen.json').read_text())
ps=json.loads((p/'planning/analysis_plan_003.frozen.json').read_text())
assert lock['primary_validation']['n_primary_tests']==1
assert pb['n_confirmatory_primary_tests']==pd['n_confirmatory_primary_tests']==ps['n_confirmatory_primary_tests']==0
assert 'effect estimation only' in pb['reporting'] and 'no rescue success criterion' in pb['reporting']
index=table('runs/004A_evidence_ledger/external_statistic_files.tsv')
mirror=[json.loads(x) for x in (O/'external_statistic_rows.jsonl').read_text().splitlines()]
icount=collections.Counter(x['source'] for x in mirror);cache={};maxerr=0.;fields=0
for z in index:
 f=p/z['path'];assert hashlib.sha256(f.read_bytes()).hexdigest()==z['sha256']
 cache[z['path']]=table(z['path'])
 assert len(cache[z['path']])==int(z['rows'])==icount[z['path']]
for m in mirror:
 old=cache[m['source']][m['source_row_1based']-1]
 assert set(old)==set(m['statistics'])
 for k,v in m['statistics'].items():
  rawvalue=old[k];fields+=1
  if v is None:assert rawvalue in ['','NA','NaN','nan','N/A']
  elif isinstance(v,bool):assert rawvalue.lower()==str(v).lower()
  elif isinstance(v,(float,int)):
   maxerr=max(maxerr,abs(float(rawvalue)-v))
   assert math.isclose(float(rawvalue),v,rel_tol=1e-12,abs_tol=5e-15)
  else:assert rawvalue==v
man=json.loads((p/'provenance/manifests/stage_004A_evidence_ledger.json').read_text())
hashchecks=[]
def walk(z):
 if isinstance(z,dict):
  if 'path' in z and 'sha256' in z:
   f=Path(z['path']);hashchecks.append({'path':str(f),'sha256_match':hashlib.sha256(f.read_bytes()).hexdigest()==z['sha256'],'size_match':f.stat().st_size==z['bytes']})
  for v in z.values():walk(v)
 elif isinstance(z,list):
  for v in z:walk(v)
walk(man)
assert all(z['sha256_match'] and z['size_match'] for z in hashchecks)
compat=json.loads((O/'interval_compatibility.json').read_text())
assert float(w6['lower'])<=float(b['estimate'])<=float(w6['upper'])
assert float(w6['lower'])<=0<=float(w6['upper'])
assert compat['VDZ_W6_OLS_CI_contains_IFX_point'] and compat['VDZ_W6_OLS_CI_contains_zero']
result={
'artifact_id':'004A_external_analysis_ledger_independent_review_v1',
'status':'REVIEWED',
'decision':'Adopt the compiled ledger and report with their explicit original evidence roles. No unresolved blocker.',
'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
'reviewer':'independent agent /root/vdz_metadata_review',
'scope':'Read-only verification of original frozen plans and executed statistical exports followed by numeric compilation checks. No regression, influence diagnostic or new scientific test run.',
'outputs_reviewed':[fm('runs/004A_evidence_ledger/'+f) for f in ['external_analysis_overview.tsv','external_statistic_files.tsv','external_statistic_rows.jsonl','inventory_scope.json','interval_compatibility.json']]+[fm('reports/External_analysis_ledger_and_interpretation_20260911.md')],
'producer_code':fm('code/build_004A_external_analysis_ledger.py'),
'producer_manifest_snapshot':fm('provenance/manifests/stage_004A_evidence_ledger.json'),
'reviewer_code':fm('provenance/reviews/004A_external_analysis_ledger_review.py'),
'counts_verified_from_patient_sources':counts,
'count_source_files':[fm(f) for f in ['runs/002A_GSE23597/primary_patient_data.tsv','runs/002B_GSE73661_support/support_patient_data.tsv','runs/002D_VDZ_support/patient_level_scores.tsv','runs/003_cell_context/joint_primary_patient_metrics.tsv']],
'numeric_checks':{'all_five_headline_estimates_intervals_P_and_n_match':True,'CT_two_axis_BH_separate_from_raw_fixed15_and_DecontX_fixed15':True,'fixed_state_family':15,'evaluable_states':10,'included_tables':len(index),'indexed_rows':len(mirror),'mirrored_fields_compared':fields,'maximum_serialization_error':maxerr,'mirrored_rows_are_not_independent_test_count':True},
'evidence_role_checks':{'original_confirmatory_primary_test_count':1,'002B_effect_estimation_without_success_criterion':True,'002D_W6_principal_supportive_not_confirmatory':True,'002D_W12_secondary_not_replacement':True,'003_exploratory':True,'DecontX_same_patient_sensitivity_not_sixth_validation':True,'explicit_post_result_BH15_implementation_amendment_disclosed':True},
'retained_context':{'composition_main19_and_joint16_in_inventory':True,'threshold10_exact_nonsignificant_and_HC3_nominal_positive_disclosed':True,'all15_states_including5_non_evaluable_retained':True,'TA_LGR5_raw_positive_and_DecontX_BH_nonsignificant_retained':True,'cross_sectional_sample_state_contamination_descriptive_only':True,'repeated_source_P_values_not_new_tests':True,'original_primary_negative_retained':True},
'interpretation_checks':{'VDZ_W6_OLS_CI_contains_IFX_point':True,'VDZ_W6_OLS_CI_contains_zero':True,'containment_not_equality_or_equivalence_test':True,'HC3_not_assumed_more_conservative_than_OLS':True,'HC3_OLS_difference_not_used_to_infer_heteroskedasticity_or_patient_dominance':True,'no_effect_size_ranking_across_bulk_and_singlecell_units':True,'no_global_any_positive_success_claim':True},
'resolved_comment':'Producer replaced hardcoded outcome counts with source-table calculations, added count inputs to manifest, and retained threshold10 HC3 and joint-patient comparison limitations.',
'manifest_checks':hashchecks,
'unresolved_blockers':[],
'new_scientific_tests':False,
'old_reports_or_frozen_plans_modified_by_reviewer':False,
'environment':{'python':platform.python_version()}
}
out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'review_path':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'tables':len(index),'rows':len(mirror),'manifest_entries':len(hashchecks),'unresolved_blockers':0}))
