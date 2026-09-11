from pathlib import Path
import json, hashlib, itertools, math, datetime
import numpy as np
import pandas as pd
from scipy import stats
from functools import lru_cache
P=Path('/root/projects/UC_Treatment_Recovery')
O=P/'provenance/reviews'
genes=['AQP8','HMGCS2','GUCA2A','CA2','SLC26A3','MS4A12']
counts=[g+'_counts' for g in genes]; norms=[g+'_norm10k_sum' for g in genes]
metrics=['delta_ct_logit_epi','delta_ct_logit_all','delta_ct_state_score','delta_epi_score','composition_contribution','within_state_contribution','total_epi_linear_change']
ag=pd.read_csv(P/'runs/003_extraction/sample_state_candidate_aggregates.tsv.gz',sep='\t')
pairs=pd.read_csv(P/'runs/003_preexpression/fixed_UC_colonic_site_pairs.tsv',sep='\t')
fullkeys=['cell_set','sample_id','Patient','Site','Treatment','Remission_status','LibraryType','Batch','final_analysis','major','minor']
assert ag.duplicated(fullkeys).sum()==0
assert pairs.duplicated(['patient','site']).sum()==0
assert len(pairs.pre_sample_id.unique())==len(pairs)==len(pairs.post_sample_id.unique())
assert set(pairs.pre_sample_id).isdisjoint(set(pairs.post_sample_id))
@lru_cache(None)
def comb(n,k): return np.array(list(itertools.combinations(range(n),k)),dtype=np.int16)
def fit(d,metric,branch):
 z=d[d[metric].notna()].copy(); y=z[metric].to_numpy(float); g=(z.remission=='Remission').to_numpy(int); n=len(y); k=int(g.sum())
 out=dict(branch=branch,metric=metric,n=n,remission=k,nonremission=n-k)
 if n<4 or k in (0,n): return out
 y1=y[g==1];y0=y[g==0];mean1=float(np.mean(y1));mean0=float(np.mean(y0));b=mean1-mean0
 df=n-2; sse=np.sum((y1-mean1)**2)+np.sum((y0-mean0)**2)
 se=np.sqrt(sse/df*(1/k+1/(n-k)))
 # Independent HC3 using closed-form two-group covariance, rather than production matrix sandwich.
 seh=np.sqrt(np.sum((y1-mean1)**2)/(k-1)**2+np.sum((y0-mean0)**2)/(n-k-1)**2)
 cr=stats.t.ppf(.975,df)
 sums=y[comb(n,k)].sum(axis=1); diff=sums/k-(y.sum()-sums)/(n-k)
 pe=float(np.mean(np.abs(diff)>=abs(b)-1e-12))
 out.update(mean_remission=mean1,mean_nonremission=mean0,estimate=b,se=float(se),lower=float(b-cr*se),upper=float(b+cr*se),p_OLS=float(2*stats.t.sf(abs(b/se),df)),se_HC3=float(seh),lower_HC3=float(b-cr*seh),upper_HC3=float(b+cr*seh),p_HC3=float(2*stats.t.sf(abs(b/seh),df)),p_exact_permutation=pe,df=df,permutation_assignments=math.comb(n,k))
 return out
samples={}
for (cs,sid),a in ag.groupby(['cell_set','sample_id'],sort=True):
 meta=a[['Patient','Site','Treatment','Remission_status','LibraryType','Batch']].drop_duplicates(); assert len(meta)==1
 e=a[a.major=='Non_ileal_epithelium'];ct=a[a.final_analysis=='Non ileal CT colonocyte'];assert (ct.major=='Non_ileal_epithelium').all()
 def sums(q): return dict(n=int(q.n_cells.sum()),umi=float(q.total_UMI.sum()),count=q[counts].sum().to_numpy(float),norm=q[norms].sum().to_numpy(float))
 s=dict(meta.iloc[0]); s.update(all=sums(a),epi=sums(e),ct=sums(ct));s['other_n']=s['epi']['n']-s['ct']['n']
 for den,name in [('epi','ct_logit_epi'),('all','ct_logit_all')]:
  n=s[den]['n']; s[name]=float(np.log((s['ct']['n']+.5)/(n-s['ct']['n']+.5))) if n>0 else np.nan
 for den,name in [('ct','ct_state_score'),('epi','epi_score')]: s[name]=float(np.log2((s[den]['count']+.5)/(s[den]['umi']+1)*1e6).mean())
 s['ct_linear']=float(s['ct']['norm'].mean()/s['ct']['n']) if s['ct']['n']>0 else np.nan
 s['other_linear']=float((s['epi']['norm']-s['ct']['norm']).mean()/s['other_n']) if s['other_n']>0 else np.nan
 s['epi_linear']=float(s['epi']['norm'].mean()/s['epi']['n']) if s['epi']['n']>0 else np.nan
 samples[(cs,sid)]=s
sites=[]
for cs in ['exclude_predicted_doublets','author_annotated']:
 for th in [10,20,50]:
  for scope in ['all_fixed','v31_only','author_both']:
   pp=pairs.copy()
   if scope=='v31_only': pp=pp[pp.library_type=='3 prime v3.1']
   if scope=='author_both': pp=pp[pp.both_in_author_paired_list.astype(str).str.lower()=='true']
   for r in pp.to_dict('records'):
    a=samples[(cs,r['pre_sample_id'])];b=samples[(cs,r['post_sample_id'])]
    for s,visit in [(a,'Pre'),(b,'Post')]:
     assert s['Patient']==r['patient'] and s['Site']==r['site'] and s['Remission_status']==r['remission'] and s['Treatment']==visit and s['LibraryType']==r['library_type'] and s['Batch']==r['batch']
    d={**r,'cell_set':cs,'scope':scope,'threshold':th}
    for met in ['ct_logit_epi','ct_logit_all','epi_score']:d['delta_'+met]=b[met]-a[met]
    ok=a['ct']['n']>=th and b['ct']['n']>=th
    d.update(pre_ct_n=a['ct']['n'],post_ct_n=b['ct']['n'],pre_epi_n=a['epi']['n'],post_epi_n=b['epi']['n'],state_evaluable=ok,delta_ct_state_score=b['ct_state_score']-a['ct_state_score'] if ok else np.nan)
    for nm in ['composition_contribution','within_state_contribution','ct_within_contribution','other_within_contribution','total_epi_linear_change','decomposition_error']:d[nm]=np.nan
    if ok and a['other_n']>0 and b['other_n']>0:
     w0=a['ct']['n']/a['epi']['n'];w1=b['ct']['n']/b['epi']['n'];wb=(w0+w1)/2
     d['composition_contribution']=(w1-w0)*((a['ct_linear']+b['ct_linear'])/2-(a['other_linear']+b['other_linear'])/2)
     d['ct_within_contribution']=wb*(b['ct_linear']-a['ct_linear'])
     d['other_within_contribution']=(1-wb)*(b['other_linear']-a['other_linear'])
     d['within_state_contribution']=d['ct_within_contribution']+d['other_within_contribution']
     d['total_epi_linear_change']=b['epi_linear']-a['epi_linear']
     d['decomposition_error']=d['total_epi_linear_change']-d['composition_contribution']-d['within_state_contribution']
    sites.append(d)
S=pd.DataFrame(sites)
keys=['cell_set','scope','threshold','patient','remission','library_type','batch']
extra=['ct_within_contribution','other_within_contribution','decomposition_error']
D=S.groupby(keys,dropna=False)[metrics+extra].mean().reset_index()
D['n_site_pairs']=S.groupby(keys,dropna=False).size().values
D['n_joint_axis_site_pairs']=S.assign(ok=S.delta_ct_logit_epi.notna()&S.delta_ct_state_score.notna()).groupby(keys,dropna=False).ok.sum().values
models=[]
for (cs,sc,th),d in D.groupby(['cell_set','scope','threshold']):
 for m in metrics: models.append(fit(d,m,f'{cs}::{sc}::threshold{th}'))
M=pd.DataFrame(models);M['BH_two_main']=np.nan
ix=M.index[(M.branch=='exclude_predicted_doublets::all_fixed::threshold20')&M.metric.isin(['delta_ct_logit_epi','delta_ct_state_score'])]
p=M.loc[ix,'p_exact_permutation'].to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.minimum(q,1);M.loc[ix[order],'BH_two_main']=q
# Correct joint analysis: select jointly evaluable SITE PAIRS before patient equal-site means.
jointS=S[(S.cell_set=='exclude_predicted_doublets')&(S.scope=='all_fixed')&(S.threshold==20)&S.delta_ct_logit_epi.notna()&S.delta_ct_state_score.notna()]
J=jointS.groupby(keys,dropna=False)[metrics+extra].mean().reset_index();J['n_site_pairs']=jointS.groupby(keys,dropna=False).size().values
JM=pd.DataFrame([fit(J,m,'joint_primary_patient_set') for m in ['delta_ct_logit_epi','delta_ct_state_score','composition_contribution','within_state_contribution']])
checks={}
def compare(ind,path,keycols,cols):
 ref=pd.read_csv(P/path,sep='\t');u=ind.merge(ref,on=keycols,suffixes=('_ind','_ref'),validate='one_to_one');assert len(u)==len(ref)==len(ind),(path,len(u),len(ref),len(ind))
 errs={}
 for c in cols:
  x=u[c+'_ind'].to_numpy(float);y=u[c+'_ref'].to_numpy(float);assert np.array_equal(np.isnan(x),np.isnan(y)),(path,c,'NA mismatch');errs[c]=float(np.nanmax(np.abs(x-y))) if np.isfinite(x).any() else 0.
 return dict(rows=len(u),max_abs_difference=max(errs.values()),by_column=errs)
checks['patient_metrics']=compare(D,'runs/003_cell_context/patient_metrics.tsv',['cell_set','scope','threshold','patient'],metrics+['decomposition_error','n_site_pairs','n_joint_axis_site_pairs'])
checks['site_metrics']=compare(S,'runs/003_cell_context/site_pair_metrics.tsv',['cell_set','scope','threshold','patient','site'],metrics+['decomposition_error','pre_ct_n','post_ct_n','pre_epi_n','post_epi_n'])
checks['models']=compare(M,'runs/003_cell_context/patient_group_models.tsv',['branch','metric'],['n','remission','nonremission','mean_remission','mean_nonremission','estimate','se','lower','upper','p_OLS','se_HC3','lower_HC3','upper_HC3','p_HC3','p_exact_permutation','BH_two_main'])
# Chunk repair audit, independently summing old aggregates under the exact full released key.
old=pd.read_csv(P/'runs/003_extraction_pre_chunk_fix_invalid/sample_state_candidate_aggregates.tsv.gz',sep='\t')
sumcols=['n_cells','total_UMI']+counts+norms
old['_ng']=old.mean_n_genes*old.n_cells;old['_mt']=old.mean_pct_mt*old.n_cells
c=old.groupby(fullkeys,dropna=False)[sumcols+['_ng','_mt']].sum().reset_index();c['mean_n_genes']=c._ng/c.n_cells;c['mean_pct_mt']=c._mt/c.n_cells
u=c.merge(ag,on=fullkeys,suffixes=('_old','_new'),validate='one_to_one');assert len(u)==len(ag)==len(c)
chunkdiff={z:float(np.max(np.abs(u[z+'_old']-u[z+'_new']))) for z in sumcols+['mean_n_genes','mean_pct_mt']}
oldm=pd.read_csv(P/'runs/003_cell_context_pre_chunk_fix_invalid/patient_group_models.tsv',sep='\t')
checks['pre_post_chunk_main_models']=compare(oldm,'runs/003_cell_context/patient_group_models.tsv',['branch','metric'],['n','estimate','lower','upper','p_OLS','p_HC3','p_exact_permutation','BH_two_main'])
checks['correct_common_site_joint_vs_current']=compare(JM,'runs/003_cell_context/joint_primary_patient_models.tsv',['branch','metric'],['n','mean_remission','mean_nonremission','estimate','se','lower','upper','p_OLS','p_HC3','p_exact_permutation'])
# Family state CT keys, strict read-level invariants.
epi=ag[ag.major=='Non_ileal_epithelium'];meta=pairs[['patient','remission','library_type','batch']].drop_duplicates();assert meta.patient.nunique()==len(meta)
cell_equal=ag[ag.cell_set=='author_annotated'].drop(columns='cell_set').reset_index(drop=True).equals(ag[ag.cell_set=='exclude_predicted_doublets'].drop(columns='cell_set').reset_index(drop=True))
mainS=S[(S.cell_set=='exclude_predicted_doublets')&(S.scope=='all_fixed')&(S.threshold==20)]
attrition=[]
for th in [10,20,50]:
 ss=S[(S.cell_set=='exclude_predicted_doublets')&(S.scope=='all_fixed')&(S.threshold==th)]
 for met in ['delta_ct_logit_epi','delta_ct_state_score']:
  z=ss[ss[met].notna()];attrition.append(dict(threshold=th,metric=met,site_pairs=len(z),patients=z.patient.nunique(),remission_patients=z[z.remission=='Remission'].patient.nunique(),excluded_patients=sorted(set(ss.patient)-set(z.patient))))
for name,table in [('patient_metrics',D),('site_metrics',S),('models',M),('common_site_patient_metrics',J),('common_site_models',JM)]: table.to_csv(O/f'003_recomputed_{name}_v1.tsv',sep='\t',index=False)
num=dict(status='NUMERICALLY_RECOMPUTED',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),algorithm='Python raw aggregate reconstruction; closed-form two-group OLS and HC3; all exact group-label combinations; manual BH',software=dict(numpy=np.__version__,pandas=pd.__version__,scipy=__import__('scipy').__version__),checks=checks,aggregate=dict(rows=len(ag),full_key_columns=fullkeys,full_key_duplicates=int(ag.duplicated(fullkeys).sum()),simple_cellset_sample_finalstate_duplicates=int(ag.duplicated(['cell_set','sample_id','final_analysis']).sum()),epithelial_simple_key_duplicates=int(epi.duplicated(['cell_set','sample_id','final_analysis']).sum()),cells_per_set=ag.groupby('cell_set').n_cells.sum().to_dict(),cell_sets_identical=cell_equal,old_rows=len(old),old_full_key_duplicates=int(old.duplicated(fullkeys).sum()),reconsolidation_differences=chunkdiff),metadata=dict(site_pairs=len(pairs),patients=len(meta),crosswalk=meta.groupby(['remission','library_type','batch']).size().reset_index(name='patients').to_dict('records'),sample_metadata_pairs_all_match=True),attrition=attrition,primary=M.loc[ix].to_dict('records'),common_site_joint=JM.to_dict('records'),common_site_pairs=len(jointS),max_abs_Kitagawa_error=float(S.decomposition_error.abs().max()),common_site_ct_vs_other_terms=J.groupby('remission')[['composition_contribution','ct_within_contribution','other_within_contribution','within_state_contribution','total_epi_linear_change']].mean().reset_index().to_dict('records'))
def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,list):return [clean(v) for v in x]
 if isinstance(x,float) and not np.isfinite(x):return None
 if isinstance(x,np.generic):return x.item()
 return x
(O/'003_cell_context_numeric_independent_v1.json').write_text(json.dumps(clean(num),indent=2,allow_nan=False)+'\n')
print(json.dumps(clean(dict(primary=num['primary'],joint=num['common_site_joint'],attrition=attrition,checks=checks,aggregate=num['aggregate'],decomp=num['common_site_ct_vs_other_terms'])),indent=2,allow_nan=False))
