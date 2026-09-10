
from pathlib import Path
import json,datetime,hashlib,subprocess
r=Path('/root/projects/UC_Treatment_Recovery')
def info(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
review=r/'provenance/reviews/final_stage_adoption_001C.json'
rev=json.loads(review.read_text())
print('ADOPTION_REVIEW',json.dumps({k:v for k,v in rev.items() if k in ['status','decision','blocking_items','blocking_errors','can_freeze_candidate']},ensure_ascii=False))
assert hashlib.sha256(review.read_bytes()).hexdigest()=='eb7c4b6ac092d67e962f1ffe363a3e1d59f9f60de60b3f51a76b69d98d26d5b8'
d=json.loads((r/'planning/signature_lock.draft.json').read_text())
d['status']='FROZEN'
d['frozen_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
d['source_code_commit']=subprocess.run(['git','rev-parse','HEAD'],cwd=r,capture_output=True,text=True).stdout.strip()
d['adoption_review']=info(review)
d['freeze_scope']='Genes,weights,score,validation endpoint/model/sample rules; not a frozen claim of discovery, mechanism or validation success.'
dest=r/'signature_lock.json';assert not dest.exists()
dest.write_text(json.dumps(d,ensure_ascii=False,indent=2))
(r/'signature_lock.sha256').write_text(hashlib.sha256(dest.read_bytes()).hexdigest()+'  signature_lock.json\n')
# Append only the final adoption status; numerical text and scientific claims remain unchanged from the reviewed report.
report=r/'reports/Review_response_and_001C_report_20260910.md'
pre=info(report)
report.write_text(report.read_text()+'\n## 最终采用状态\n\n本轮已通过独立数值与文字审查，状态为REVIEWED。六基因候选及GSE23597 W8主要验证规则已写入根目录signature_lock.json并冻结；冻结范围为分析对象和验证规则，不是生物学结论。独立验证表达尚未分析。\n',encoding='utf-8')
closure={'stage':'001B2 corrections plus001C/001C2 patient-program refinement','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'REVIEWED','candidate_status':'FROZEN','signature_lock':info(dest),'report':info(report),'report_version_reviewed_before_status_append':pre,'report_append_scope':'Administrative adoption status only; numerical content unchanged','reviews':[info(r/'provenance/reviews'/x) for x in ['001B2_calibration_independent_v1.json','001C_identity_independent_v1.json','001C_001C2_programs_independent_v1.json','prior_programs_001C.json','final_stage_adoption_001C.json','stage_001C_hash_integrity.json','visual_QA_001C.json']],'run_manifests':[info(r/'provenance/manifests'/x) for x in ['001B2_calibration.json','001C_identity_qc.json','001C_programs.json','001C2_refinement.json','001C2_HC3.json','001C_figures.json']],'provenance_note':'As-generated manifests retain their original COMPLETE/PENDING fields. This adoption ledger supplies their final reviewed status without changing reviewed hashes.','independent_validation_started':False,'validation_expression_used':False,'next_phase':'002A GSE23597 technical audit followed by locked W8 test','allowed_claims':['A fixed exploratory mature-epithelial-associated candidate has conditional clinical-response association in discovery','No novel cell state or original MARS/IBrD superiority established']}
(r/'provenance/manifests/stage_001C_completion.json').write_text(json.dumps(closure,ensure_ascii=False,indent=2))
report_manifest={'artifact_id':'Review_response_and_001C_report_20260910','status':'REVIEWED','phase':'exploratory report','generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cwd':str(r),'command':'python3 code/build_report_001C.py followed by administrative adoption-status append','seed':None,'output':info(report),'code':[info(r/'code/build_report_001C.py')],'inputs':[info(r/p) for p in ['runs/001B2_calibration/Hallmark_patient_effects_CI.tsv','runs/001B2_calibration/review_numeric_audit.json','runs/001C_programs/program_overlap_and_delta_correlations.tsv','runs/001C2_refinement/refinement_response_effects_CI.tsv','runs/001C2_refinement/HC3_common_model_sensitivity.tsv','runs/001C2_refinement/refinement_leave_one_patient_out.tsv','runs/001C_programs/leave_one_candidate_gene_out.tsv','runs/001B_discovery/discovery_model_summary.json']],'environment':info(r/'renv.lock'),'review':info(review),'adoption_ledger':str(r/'provenance/manifests/stage_001C_completion.json'),'limitations':['Exploratory conditional associations; validation not executed','MARS/IBrD proxy limitations explicitly retained']}
(r/'provenance/manifests/report_001C.json').write_text(json.dumps(report_manifest,indent=2))
sup={'superseded_report':info(r/'reports/Stage1_report_20260910.md'),'current_report':info(report),'reason':'Replace significance-count framing with patient effects and CI; disclose pathway sensitivity; report new exploratory response-axis analysis','unchanged_valid_scope':'Prior source-ID pairing, expression preprocessing and under-declared-model numerical reproducibility remain valid.'}
(r/'provenance/manifests/report_supersession_001C.json').write_text(json.dumps(sup,indent=2))
status=json.loads((r/'project_status.json').read_text())
status.update({'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'current_phase':'001C completed and independently reviewed; candidate and validation rules frozen','candidate_frozen':True,'signature_lock':str(dest),'independent_validation_started':False,'review_pending':False,'stage_report':str(report),'completion_ledger':str(r/'provenance/manifests/stage_001C_completion.json'),'next_phase':'002A GSE23597 technical audit then prespecified W8 validation'})
(r/'project_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
readme=r/'README.md'
readme.write_text('# Current status: review corrections and001C complete\n\nCandidate and independent-test rules: signature_lock.json.\nCurrent report: reports/Review_response_and_001C_report_20260910.md.\nState ledger: provenance/manifests/stage_001C_completion.json.\nThe earlier report is retained for history; use the revised report for interpretation.\nIndependent validation has not been analyzed.\n\n'+readme.read_text())
print('FROZEN',info(dest));print('REPORT',info(report))
