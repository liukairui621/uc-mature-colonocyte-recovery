import pathlib,json,hashlib,datetime,sys
import pandas as pd
from scipy.stats import norm,nct,t
P=pathlib.Path('/root/projects/UC_Treatment_Recovery')
out=P/'runs/002D_VDZ_metadata';out.mkdir(exist_ok=True)
m=pd.read_csv(P/'runs/001A_metadata/sample_registry_v1.tsv',sep='\t',dtype=str).fillna('')
m=m[m.series=='GSE73661'].copy()
v=m[m.treatment.str.startswith('vdz')].copy()
vis=v.groupby(['treatment','time']).size().reset_index(name='arrays')
vis.to_csv(out/'visit_inventory.tsv',sep='\t',index=False)
rows=[]
for endpoint,paths in [('W6',['vdz_vdz4w','vdz_vdz8w','vdz_plac']),('W12',['vdz4w'])]:
 for (path,subject),z in v[v.treatment.isin(paths)].groupby(['treatment','subject']):
  b=z[z.time=='W0'];f=z[z.time==endpoint]
  assert len(b)==1 and len(f)<=1
  if len(f)==0: continue
  b=b.iloc[0];f=f.iloc[0]
  assert b.mayo_endoscopic in ['0','1','2','3'] and f.mayo_endoscopic in ['0','1','2','3']
  rows.append(dict(subject=subject,treatment_path=path,visit=endpoint,baseline_gsm=b.gsm,followup_gsm=f.gsm,baseline_mayo_endoscopic=int(b.mayo_endoscopic),followup_mayo_endoscopic=int(f.mayo_endoscopic),endoscopic_healing='Yes' if int(f.mayo_endoscopic)<=1 else 'No'))
pairs=pd.DataFrame(rows)
assert len(pairs)==40 and not pairs.subject.duplicated().any()
assert not set(pairs.subject)&set(m[m.treatment=='IFX'].subject)
pairs.to_csv(out/'VDZ_pairs.tsv',sep='\t',index=False)
use=set(pairs.baseline_gsm)|set(pairs.followup_gsm)
v['disposition']=v.gsm.map(lambda x:'in_W6_or_W12_pair' if x in use else 'not_in_target_pair')
v.to_csv(out/'VDZ_all_samples_disposition.tsv',sep='\t',index=False)
ct=pairs.groupby(['visit','endoscopic_healing']).size().unstack(fill_value=0)
ct.to_csv(out/'healing_counts.tsv',sep='\t')
pp=pd.read_csv(P/'runs/002B_precision/fitted_precision_components.tsv',sep='\t')
effect=float(pp[pp.cohort=='GSE92415'].estimate.iloc[0])
se_target=effect/(norm.ppf(.975)+norm.ppf(.8))
sc=[]
for _,r in pp.iterrows():
 n_norm=float(r.n*(r.SE/se_target)**2)
 n_t=next(n for n in range(int(r.n),1001) if nct.cdf(-t.ppf(.975,n-(int(r.n)-int(r.df))),n-(int(r.n)-int(r.df)),effect/(r.SE*(r.n/n)**.5))+nct.sf(t.ppf(.975,n-(int(r.n)-int(r.df))),n-(int(r.n)-int(r.df)),effect/(r.SE*(r.n/n)**.5))>=.8)
 sc.append(dict(source=r.cohort,assumed_effect=effect,target_SE_normal=se_target,normal_approx_required_n=n_norm,noncentral_t_required_n=n_t,assumptions='Residual variance, response proportion, and response VIF held at source values; not a transport guarantee or minimum eligible sample size'))
pd.DataFrame(sc).to_csv(out/'conditional_sample_size_scenarios.tsv',sep='\t',index=False)
r=pp[pp.cohort=='GSE23597'].iloc[0]
se62=float(r.SE*(32/62)**.5)
pow62=float(norm.cdf(-norm.ppf(.975)-effect/se62)+norm.sf(norm.ppf(.975)-effect/se62))
summary=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_W6=27,healing_W6=int(ct.loc['W6','Yes']),n_W12=13,healing_W12=int(ct.loc['W12','Yes']),distinct_VDZ_subjects=40,IFX_subject_ID_overlap=[],W6_missing_subjects=['30'],placebo_W6_n=3,conditional_62_pair_normal_power=pow62,conditional_62_pair_SE=se62,raw_response_information_IFX=23*(8/23)*(15/23),raw_response_information_VDZ_W6=27*(6/27)*(21/27),information_caveat='Before accounting for residual variance and covariates; neither a power calculation nor a drug comparison',expression_association_computed_in_this_audit=False)
(out/'audit_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print(pd.DataFrame(sc).to_string(index=False))
