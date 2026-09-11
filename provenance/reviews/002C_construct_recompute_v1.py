import pathlib,csv,json,math,collections
import numpy as np
from scipy.stats import rankdata,t,mannwhitneyu
P=pathlib.Path('/root/projects/UC_Treatment_Recovery');S=P/'runs/002B_GSE73661_support';O=P/'runs/002C_construct_overlap';R=P/'provenance/reviews'
def read(f):return list(csv.DictReader(f.open(),delimiter='\t'))
data=read(S/'support_patient_data.tsv');diagnostics=read(S/'patient_diagnostics.tsv');lopo=read(S/'leave_one_patient_out.tsv')
assert len(data)==23 and len({r['subject'] for r in data})==23
ids=[r['subject'] for r in data];yes=np.array([r['endoscopic_healing']=='Yes' for r in data]);delta=np.array([float(r['delta_candidate']) for r in data]);baseline=np.array([float(r['baseline_candidate']) for r in data]);post=baseline+delta
bendo=np.array([int(r['baseline_mayo_endoscopic']) for r in data]);pendo=np.array([int(r['followup_mayo_endoscopic']) for r in data]);dendo=pendo-bendo
assert np.array_equal(yes,pendo<=1) and yes.sum()==8
nr=int(yes.sum());nn=len(data)-nr;ties=len(np.unique(delta))!=len(delta);assert not ties
wins=float(np.sum(delta[yes,None]>delta[None,~yes])+0.5*np.sum(delta[yes,None]==delta[None,~yes]));auc=wins/(nr*nn)
ranks=rankdata(delta);rank_sum=int(ranks[yes].sum());U=rank_sum-nr*(nr+1)//2;assert U==wins
# Exact permutation null distribution via counts of rank subsets, not asymptotic normal approximation.
dp=[collections.defaultdict(int) for _ in range(nr+1)];dp[0][0]=1
for rank in range(1,len(data)+1):
 for k in range(min(rank,nr),0,-1):
  for total,count in list(dp[k-1].items()):dp[k][total+rank]+=count
nullU={ss-nr*(nr+1)//2:c for ss,c in dp[nr].items()};den=sum(nullU.values());assert den==math.comb(23,8)
lower=sum(c for u,c in nullU.items() if u<=U)/den;upper=sum(c for u,c in nullU.items() if u>=U)/den;exactp=min(1,2*min(lower,upper))
scipyMW=mannwhitneyu(delta[yes],delta[~yes],method='exact',alternative='two-sided');assert abs(scipyMW.statistic-U)<1e-12 and abs(scipyMW.pvalue-exactp)<1e-12
def corr(x,y,rank=False):
 if rank:x=rankdata(x,method='average');y=rankdata(y,method='average')
 xc=x-x.mean();yc=y-y.mean();r=float(xc@yc/math.sqrt((xc@xc)*(yc@yc)));p=float(2*t.sf(abs(r)*math.sqrt((len(x)-2)/(1-r*r)),len(x)-2));return r,p
spdd=corr(delta,dendo,True);sppost=corr(post,pendo,True);spbase=corr(baseline,bendo,True);ppost=corr(post,pendo);pbase=corr(baseline,bendo)
expected=[{'analysis':'AUC_healing_from_delta_candidate','estimate':auc,'p':None},{'analysis':'Mann_Whitney_delta_by_healing','estimate':U,'p':exactp}]
for label,stat in [('Spearman_delta_candidate_vs_delta_endoscopic',spdd),('Spearman_post_candidate_vs_post_endoscopic',sppost),('Spearman_baseline_candidate_vs_baseline_endoscopic',spbase),('Pearson_post_candidate_vs_post_endoscopic',ppost),('Pearson_baseline_candidate_vs_baseline_endoscopic',pbase)]:expected.append({'analysis':label,'estimate':stat[0],'p':stat[1]})
observed={r['analysis']:r for r in read(O/'construct_overlap_statistics.tsv')};assert set(observed)=={r['analysis'] for r in expected};mxstat=0
for row in expected:
 b=observed[row['analysis']];assert int(b['n'])==23;mxstat=max(mxstat,abs(float(b['estimate'])-row['estimate']))
 if row['p'] is not None:mxstat=max(mxstat,abs(float(b['p'])-row['p']))
assert mxstat<1e-10
# Reconstruct fixed descriptive support model only to check influence and leave-one-out signs.
infl=np.array([float(r['delta_common194_inflammation']) for r in data]);X=np.column_stack([np.ones(23),yes.astype(float),baseline,infl])
def fit(X,y):
 u,s,v=np.linalg.svd(X,full_matrices=False);rank=np.linalg.matrix_rank(X);assert rank==X.shape[1];beta=v.T@((u.T@y)/s);e=y-X@beta;df=len(y)-rank;mse=e@e/df;h=np.sum(u*u,axis=1);bread=(v.T/(s*s))@v
 se=np.sqrt(np.diag(bread)*mse);pv=float(2*t.sf(abs(beta[1]/se[1]),df));cook=e*e/(rank*mse)*h/(1-h)**2;ext=np.sqrt((e@e-e*e/(1-h))/(df-1));student=e/(ext*np.sqrt(1-h))
 return {'beta':beta,'p':pv,'resid':e,'h':h,'cook':cook,'student':student}
f=fit(X,delta);diagby={r['subject']:r for r in diagnostics};lpby={r['omitted_subject']:r for r in lopo};mxdiag=0;mxlo=0;ownlo={}
for i,s in enumerate(ids):
 for field,vec in [('leverage',f['h']),('cooks_distance',f['cook']),('residual',f['resid']),('studentized_residual',f['student'])]:mxdiag=max(mxdiag,abs(float(diagby[s][field])-vec[i]))
 mask=np.arange(23)!=i;z=fit(X[mask],delta[mask]);ownlo[s]={'estimate':float(z['beta'][1]),'p':z['p']};mxlo=max(mxlo,abs(float(lpby[s]['estimate'])-z['beta'][1]),abs(float(lpby[s]['p'])-z['p']))
assert mxdiag<1e-10 and mxlo<1e-10
sourceby={r['subject']:r for r in read(O/'patient_level_source.tsv')};assert set(sourceby)==set(ids);mxmerged=0
for i,s in enumerate(ids):
 m=sourceby[s]
 for key in data[i]:
  if key in ['baseline_mayo_endoscopic','followup_mayo_endoscopic','delta_candidate','baseline_candidate','delta_common194_inflammation']:mxmerged=max(mxmerged,abs(float(m[key])-float(data[i][key])))
  else:assert m[key]==data[i][key]
 for key,val in [('post_candidate',post[i]),('delta_endoscopic',dendo[i])]:mxmerged=max(mxmerged,abs(float(m[key])-val))
 for key in ['leverage','cooks_distance','residual','studentized_residual']:mxmerged=max(mxmerged,abs(float(m[key])-float(diagby[s][key])))
 for key in ['estimate','p']:mxmerged=max(mxmerged,abs(float(m[key])-float(lpby[s][key])))
assert mxmerged<1e-10
groups=[]
for label,mask in [('healing_yes',yes),('healing_no',~yes)]:groups.append({'group':label,'n':int(mask.sum()),'delta_min':float(delta[mask].min()),'delta_max':float(delta[mask].max()),'n_above_healing_min':int(np.sum(delta[mask]>=delta[yes].min()))})
obsgrp={r['group']:r for r in read(O/'group_delta_ranges.tsv')}
for g in groups:
 for k in g:
  if k!='group':assert abs(float(obsgrp[g['group']][k])-g[k])<1e-10
summary=json.loads((O/'construct_overlap_summary.json').read_text());assert summary['baseline_endoscopic_distribution']=={str(k):v for k,v in collections.Counter(bendo.tolist()).items()}
above=[ids[i] for i in range(23) if not yes[i] and delta[i]>=delta[yes].min()];assert above==[str(x) for x in summary['nonhealing_above_healing_min']]
maxCook=ids[int(np.argmax(f['cook']))];i43=ids.index('43');assert maxCook=='43';assert abs(summary['max_Cooks_D']-f['cook'][i43])<1e-10;assert abs(summary['LOPO_after_max_Cooks_D_omission']-ownlo['43']['estimate'])<1e-10
result={'status':'RECOMPUTATION_PASSED','statistics':expected,'AUC_pairwise_wins':U,'AUC_pair_comparisons':nr*nn,'MW_rank_sum':rank_sum,'MW_exact_null_assignments':den,'MW_upper_tail_assignments':sum(c for u,c in nullU.items() if u>=U),'MW_no_score_ties':True,'endoscopic_change_definition':'followup minus baseline; negative values are endoscopic improvement','Spearman_method':'Pearson correlation of average ranks, asymptotic t reference with n-2 df because ordinal endpoint ties; not exact Spearman P','baseline_Mayo_distribution':dict(collections.Counter(map(str,bendo))),'baseline_Mayo_by_healing':{label:dict(collections.Counter(map(str,bendo[mask]))) for label,mask in [('Yes',yes),('No',~yes)]},'group_ranges':groups,'nonhealing_at_or_above_healing_min':above,'patient43':{'healing':data[i43]['endoscopic_healing'],'endoscopic_baseline':int(bendo[i43]),'endoscopic_followup':int(pendo[i43]),'delta_candidate':float(delta[i43]),'CookD':float(f['cook'][i43]),'leverage':float(f['h'][i43]),'all_patient_adjusted_beta':float(f['beta'][1]),'beta_without_patient43':ownlo['43']['estimate'],'p_without_patient43':ownlo['43']['p'],'coefficient_change_on_omission':ownlo['43']['estimate']-float(f['beta'][1]),'direction':'Omitting patient43 increases the positive coefficient; patient43 does not create or solely drive the positive association.'},'checks':{'statistics7_rows_max_diff':mxstat,'patient_diagnostics23_rows_max_diff':mxdiag,'LOPO23_rows_max_diff':mxlo,'merged_source23_rows_max_diff':mxmerged},'limitations':['AUC is empirical contemporaneous separation using post-minus-baseline expression, not prospective baseline prediction.','These construct diagnostics are post-result descriptive and do not constitute independently prespecified confirmation.','Baseline endoscopy has only2/3 with19of23 at3; lack of baseline correlation does not establish absence of biological association.','High AUC and correlations are compatible with shared construct content, not measurement equivalence or causal identity.']}
with (R/'002C_construct_numeric_independent_v1.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
