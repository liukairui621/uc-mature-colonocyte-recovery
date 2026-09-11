from pathlib import Path
import csv,json,collections,subprocess,hashlib,datetime,sys,platform
root=Path('/root/projects/UC_Treatment_Recovery')
out=root/'provenance/reviews/002D_VDZ_metadata_review.json'
if out.exists(): raise SystemExit('Refusing to overwrite existing review')
header=root/'inputs/metadata/GSE73661.matrix_header.txt'
rows=list(csv.reader(header.read_text().splitlines(),delimiter='\t'))
cs=[r for r in rows if r and r[0]=='!Sample_characteristics_ch1']
gsms=next(r[1:] for r in rows if r and r[0]=='!Sample_geo_accession')
samples=[]
for i,gsm in enumerate(gsms):
 d={'gsm':gsm}
 for r in cs:
  k,v=r[i+1].split(': ',1);d[k]=v
 samples.append(d)
T='induction therapy_maintenance therapy';I='study individual number';W='week (w)';M='mayo endoscopic subscore'
vdz={d[I] for d in samples if 'vdz' in d[T]}
ifx={d[I] for d in samples if d[T]=='IFX'}
pairsets={}
for week in ['W6','W12']:
 pairs=[]
 for pid in sorted(vdz,key=int):
  b=[d for d in samples if d[I]==pid and d[W]=='W0']
  f=[d for d in samples if d[I]==pid and d[W]==week]
  if b and f:
   assert len(b)==len(f)==1
   pairs.append({'subject':pid,'path':f[0][T],'baseline_gsm':b[0]['gsm'],'followup_gsm':f[0]['gsm'],'baseline_mayo':int(b[0][M]),'followup_mayo':int(f[0][M]),'healing':int(f[0][M])<=1})
 pairsets[week]={'n_pairs':len(pairs),'n_healing':sum(d['healing'] for d in pairs),'n_nonhealing':sum(not d['healing'] for d in pairs),'pairs':pairs}
r='suppressPackageStartupMessages(library(jsonlite)); obj<-readRDS("/root/projects/UC_Treatment_Recovery/runs/002B_GSE73661_support/support_scores.rds"); cat(toJSON(list(object_names=names(obj),score_dimensions=dim(obj$sample_scores),gene_expression_dimensions=dim(obj$gene_expression),score_gsm_names=colnames(obj$sample_scores),IFX_pair_count=nrow(obj$pairs),IFX_pair_subjects=as.character(obj$pairs$subject),R_version=R.version.string),auto_unbox=TRUE))'
proc=subprocess.run(['Rscript','--vanilla','-e',r],capture_output=True,text=True,check=True)
structure=json.loads(proc.stdout)
assert structure['score_dimensions']==[5,178]
assert set(structure['score_gsm_names'])==set(gsms)
assert set(structure['IFX_pair_subjects'])==ifx
inputs=['inputs/metadata/GSE73661.matrix_header.txt','runs/002B_GSE73661_support/downloaded_matrix_header.txt','code/support_73661_002B.R','runs/002B_GSE73661_support/support_scores.rds','planning/Study_Design_v1.md','planning/support_plan_002B.frozen.json','renv.lock']
def filemeta(rel):
 p=root/rel
 return {'path':rel,'size':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
result={
'artifact_id':'002D_VDZ_metadata_independent_review_v1',
'run_id':'002D_VDZ_extension_metadata_only',
'reviewer':'independent agent /root/vdz_metadata_review',
'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
'scope':'Metadata counts, prior scoring scope, and interpretation of expansion design only; no new VDZ expression values, scores, deltas or outcome models computed or displayed.',
'status':'REVIEWED',
'plan_execution_approval':'Not assessed: frozen 002D plan not yet provided to this reviewer',
'inputs':[filemeta(f) for f in inputs],
'methods':{'metadata':'Python csv parser of published GEO matrix header; exact subject IDs and scheduled week labels; healing defined as released Mayo endoscopic subscore <=1','RDS':'Read existing object; inspect only names, dimensions, sample IDs and IFX pair IDs; no score values extracted or printed','code_evidence':{'whole_matrix_scoring_line':37,'IFX_pair_subset_line':38,'save_full_object_line':39},'random_seed':'not applicable; deterministic metadata review','command':'python3 provenance/reviews/002D_VDZ_metadata_review.py','python':platform.python_version(),'r_structure_command':['Rscript','--vanilla','-e',r]},
'counts':{'all_samples':len(samples),'active_VDZ_unique_subjects':len(vdz),'IFX_unique_subjects':len(ifx),'all_trial_subjects_including_placebo':len({d[I] for d in samples if d[T] not in ['IFX','CO']}),'VDZ_IFX_subject_overlap':sorted(vdz&ifx),'week6_week12_subject_overlap':sorted({p['subject'] for p in pairsets['W6']['pairs']}&{p['subject'] for p in pairsets['W12']['pairs']}),'paths_by_week':{t:dict(collections.Counter(d[W] for d in samples if d[T]==t)) for t in sorted({d[T] for d in samples})},'pairs':pairsets},
'prior_object_structure':structure,
'findings':[
{'severity':'BLOCKING_if_unfixed','id':'prior_score_access','finding':'The claim that VDZ sample scores were never computed is false. Existing 002B code computed all 178 sample scores and saved them in support_scores.rds. Only IFX paired differences and outcome association were the target of the inspected 002B code.','required_resolution':'Disclose matrix and sample scores already computed; describe plan as prospective to new VDZ paired outcome-association analysis, not prospective to expression or scoring access. Do not claim that source files prove no human ever inspected the VDZ values.'},
{'severity':'MAJOR','id':'small_healing_group','finding':'W6 has 27 pairs but only 6 healing; W12 13 pairs but only 3 healing. Total pairs alone overstates precision and model support.'},
{'severity':'MAJOR','id':'future_path_conditioning','finding':'All 3 vdz_plac W6 samples are healing, accounting for half the W6 healing group. Excluding them by future maintenance path would substantially select the endpoint group. Future maintenance treatment is not a baseline covariate.'},
{'severity':'MAJOR','id':'W12_distinct_group','finding':'W12 and W6 subjects do not overlap. Do not describe W12 as an extended follow-up of W6 participants or freely pool them; verify LTS exposure histories with original trial sources.'},
{'severity':'MAJOR','id':'power_not_hard_threshold','finding':'The 72-128 pair calculation is a conditional planning approximation under fixed beta, residual variance, response balance and design collinearity. It is not a minimum validity criterion and does not establish impossibility of clinical-response validation. The fixed beta was selected in discovery and remains uncertain.'},
{'severity':'MAJOR','id':'negative_VDZ_not_specificity','finding':'A nonsignificant VDZ association does not demonstrate anti-TNF specificity; cross-drug effect differences require compatible estimands and an interaction/heterogeneity assessment, which this cross-study design may not identify causally.'},
{'severity':'MAJOR','id':'positive_VDZ_not_common_mechanism','finding':'A positive VDZ association would support cross-drug endoscopic association in sampled patients, not establish a universal recovery mechanism or drug efficacy ranking.'},
{'severity':'MAJOR','id':'same_dataset_evidence_role','finding':'Nonoverlapping subject IDs support distinct patients within the same GEO study. This is cross-drug patient evidence sharing study, platform and processing context; it is not a new external dataset.'}
],
'boundaries':['Keep GSE23597 failed clinical-response primary test unchanged.','New endoscopic analysis cannot replace the clinical-response primary or be pooled with its P value.','Keep candidate genes and scoring rules fixed.','Retain vdz_plac W6 as induction exposure after checking original trial timing.','Do not classify all W52 samples as continuous VDZ.','List off-window visits separately without relabeling to W6 or W12.','Placebo W6 has only 3 pairs; avoid credible drug-placebo efficacy claims.','Retain endpoint-score construct overlap and bulk tissue composition limitations.','GSE16879 independence remains unverified until patient-level evidence resolves potential overlap; shared institution is not proof of duplicate patients.'],
'new_biological_models_run':False,
'new_VDZ_scores_computed':False,
'expression_values_displayed':False
}
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'review_path':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'W6_n':pairsets['W6']['n_pairs'],'W6_healing':pairsets['W6']['n_healing'],'W12_n':pairsets['W12']['n_pairs'],'W12_healing':pairsets['W12']['n_healing'],'overlap_VDZ_IFX':len(vdz&ifx)}))
