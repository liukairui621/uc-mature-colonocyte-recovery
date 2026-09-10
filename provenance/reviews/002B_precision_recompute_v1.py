import pathlib,json,csv,math,hashlib
import numpy as np
import scipy
from scipy.stats import t,nct,chi2
from scipy.optimize import brentq
P=pathlib.Path('/root/projects/UC_Treatment_Recovery');O=P/'runs/002B_precision';R=P/'provenance/reviews'
raw=json.loads((R/'002B_precision_recompute_inputs_v1.json').read_text())
def power(effect,se,df):
 crit=t.ppf(.975,df);return float(nct.sf(crit,df,-effect/se)+nct.sf(crit,df,effect/se))
def positive_support(effect,se,df):return float(nct.sf(t.ppf(.975,df),df,effect/se))
def mde(se,df,positive=False):
 fun=positive_support if positive else power
 return float(brentq(lambda e:fun(e,se,df)-.8,0,20*se,xtol=1e-12))
models={};infos={}
for name,src in raw['models'].items():
 X=np.asarray(src['X'],float);y=np.asarray(src['y'],float);cols=src['columns'];n,p=X.shape
 U,s,Vh=np.linalg.svd(X,full_matrices=False);rank=int(np.linalg.matrix_rank(X));assert rank==p
 beta=Vh.T@((U.T@y)/s);res=y-X@beta;df=n-rank;sigma=float(np.linalg.norm(res)/math.sqrt(df));j=cols.index('responseYes')
 covariance=(Vh.T/(s*s))@Vh;se=float(sigma*math.sqrt(covariance[j,j]))
 response=X[:,j];Z=np.delete(X,j,axis=1);rr=response-Z@np.linalg.lstsq(Z,response,rcond=None)[0];sst=float(np.sum((response-response.mean())**2));rss=float(rr@rr)
 assert abs(se-sigma/math.sqrt(rss))<1e-10
 models[name]={'X':X,'y':y,'beta':beta,'columns':cols,'response_index':j}
 lo=float(sigma*math.sqrt(df/chi2.ppf(.975,df)));hi=float(sigma*math.sqrt(df/chi2.ppf(.025,df)))
 info={'cohort':name,'n':n,'p_response':float(response.mean()),'df':df,'delta_SD':float(np.std(y,ddof=1)),'residual_sigma':sigma,'residual_sigma_lower':lo,'residual_sigma_upper':hi,'response_SST':sst,'residualized_response_SS':rss,'response_VIF':sst/rss,'SE':se,'estimate':float(beta[j]),'MDE80_plugin':mde(se,df),'MDE80_lower_sigma':mde(se*lo/sigma,df),'MDE80_upper_sigma':mde(se*hi/sigma,df)}
 infos[name]=info
D=infos['GSE92415'];V=infos['GSE23597']
def rows(name):return list(csv.DictReader((O/name).open(),delimiter='\t'))
checks=[]
def compare(name,expected,keys,columns,tol=1e-7):
 observed=rows(name)
 def key(r):return tuple(round(float(r[k]),10) if isinstance(r[k],(int,float,np.number)) or k in ['assumed_effect','assumed_true_effect','n','assumed_response_VIF'] else r[k] for k in keys)
 obs={key(r):r for r in observed};exp={key(r):r for r in expected};assert len(obs)==len(observed) and len(exp)==len(expected) and set(obs)==set(exp),(name,set(obs)^set(exp))
 mx=0.0
 for k,r in exp.items():
  for c in columns:mx=max(mx,abs(float(r[c])-float(obs[k][c])))
 assert mx<tol,(name,mx)
 checks.append({'file':name,'rows':len(expected),'max_absolute_numeric_difference':mx})
compare('fitted_precision_components.tsv',[D,V],['cohort'],[k for k in D if k!='cohort'])
factors={'sample_size':math.sqrt(D['n']/V['n']),'response_balance':math.sqrt(D['p_response']*(1-D['p_response'])/(V['p_response']*(1-V['p_response']))),'response_collinearity':math.sqrt(V['response_VIF']/D['response_VIF']),'residual_sigma':V['residual_sigma']/D['residual_sigma']}
compare('exact_SE_ratio_decomposition.tsv',[{'component':k,'factor':v} for k,v in factors.items()],['component'],['factor'])
ratio=V['SE']/D['SE'];assert abs(math.prod(factors.values())-ratio)<1e-10
curve=[]
for z in [D,V]:
 for effect in np.arange(181)/100:curve.append({'cohort':z['cohort'],'assumed_true_effect':float(effect),'power':power(effect,z['SE'],z['df']),'SE_plugin':z['SE'],'df':z['df']})
compare('conditional_power_scenarios.tsv',curve,['cohort','assumed_true_effect'],['power','SE_plugin','df'])
sel=[{'assumed_effect':e,'validation_power':power(e,V['SE'],V['df'])} for e in [.25,.4,D['estimate'],.8,1.0]]
compare('selected_power_scenarios.tsv',sel,['assumed_effect'],['validation_power'])
obs_power=power(V['estimate'],V['SE'],V['df']);compare('observed_power_arithmetic_only.tsv',[{'effect':V['estimate'],'observed_power':obs_power}],[],['effect','observed_power'])
gd=np.array(raw['gene_deltas']['discovery']);gv=np.array(raw['gene_deltas']['validation']);assert gd.shape==(6,65) and gv.shape==(6,32)
sdD=np.r_[np.std(gd,axis=1,ddof=1),np.std(gd.mean(axis=0),ddof=1)];sdV=np.r_[np.std(gv,axis=1,ddof=1),np.std(gv.mean(axis=0),ddof=1)]
sdrows=[{'feature':g,'discovery_delta_SD':float(a),'validation_delta_SD':float(b),'ratio':float(b/a)} for g,a,b in zip(raw['genes']+['Six_gene_mean'],sdD,sdV)]
compare('observed_delta_dispersion.tsv',sdrows,['feature'],['discovery_delta_SD','validation_delta_SD','ratio'])
varrows=[]
for cohort,G in [('discovery',gd),('validation',gv)]:
 C=np.cov(G,ddof=1);diag=float(np.trace(C)/36);cross=float((C.sum()-np.trace(C))/36);scorevar=float(np.var(G.mean(axis=0),ddof=1));assert abs(scorevar-diag-cross)<1e-12
 varrows.append({'cohort':cohort,'variance_of_score_mean':scorevar,'sum_individual_variances_over_k2':diag,'sum_cross_covariances_over_k2':cross,'mean_gene_pair_correlation':float(np.corrcoef(G)[np.triu_indices(6,1)].mean())})
compare('score_variance_identity.tsv',varrows,['cohort'],list(varrows[0])[1:])
att=[]
for cohort,before,after in [('GSE92415',infos['discovery_before_inflammation']['estimate'],D['estimate']),('GSE23597',infos['validation_before_inflammation']['estimate'],V['estimate'])]:att.append({'cohort':cohort,'before':before,'after':after,'attenuation_percent':100*(1-after/before)})
compare('inflammation_conditioning_attenuation.tsv',att,['cohort'],['before','after','attenuation_percent'])
pairs=list(csv.DictReader((O/'GSE73661_IFX_pairs_metadata_independent_v1.tsv').open(),delimiter='\t'));assert len(pairs)==23 and sum(r['endoscopic_healing']=='Yes' for r in pairs)==8
N=23;prop=8/23;df=19
grid=[]
for sigma_src in ['GSE23597','GSE92415']:
 for vif in [1,1.25,1.5,2]:
  sig=infos[sigma_src]['residual_sigma'];se=sig*math.sqrt(vif/(N*prop*(1-prop)))
  grid.append({'sigma_source':sigma_src,'assumed_response_VIF':vif,'residual_sigma':sig,'SE':se,'MDE80':mde(se,df),'power_at_discovery_0_599':power(D['estimate'],se,df),'n':N,'n_healing':8,'n_nonhealing':15,'expected_df':df})
compare('GSE73661_preexpression_precision_scenarios.tsv',grid,['sigma_source','assumed_response_VIF'],list(grid[0])[2:])
future=[]
for n in [23,40,60,80,100,150,200]:
 for sigma_src in ['GSE23597','GSE92415']:
  for vif in [1.25,1.5]:
   sig=infos[sigma_src]['residual_sigma'];se=sig*math.sqrt(vif/(n*prop*(1-prop)))
   future.append({'n':n,'sigma_source':sigma_src,'assumed_response_VIF':vif,'residual_sigma':sig,'assumed_response_fraction':prop,'SE':se,'MDE80':mde(se,n-4),'power_at_0_599':power(D['estimate'],se,n-4)})
compare('future_sample_size_scenarios.tsv',future,['n','sigma_source','assumed_response_VIF'],list(future[0])[3:])
crit=t.ppf(.975,V['df']);direction={'at_assumed_discovery_effect':{'effect':D['estimate'],'two_sided_rejection_power':power(D['estimate'],V['SE'],V['df']),'positive_support_probability':positive_support(D['estimate'],V['SE'],V['df']),'negative_rejection_probability':float(nct.sf(crit,V['df'],-D['estimate']/V['SE']))},'MDE80_two_sided_rejection':V['MDE80_plugin'],'MDE80_positive_support':mde(V['SE'],V['df'],True),'at_zero_effect':{'two_sided_rejection':power(0,V['SE'],V['df']),'positive_support':positive_support(0,V['SE'],V['df'])}}
res={'status':'INDEPENDENT_RECOMPUTATION_PASSED','software':{'numpy':np.__version__,'scipy':scipy.__version__},'methods':'Source original fit matrices/responses exported with R; independent NumPy SVD estimates/residual SD, projection-based response VIF, SciPy noncentral-t CDF/SF, chi-square quantiles and Brent root solver; no GSE73661 expression read.','components':[D,V],'factor_decomposition':factors,'SE_ratio':ratio,'SD_ratios':sdrows,'variance_identity':varrows,'descriptive_attenuation':att,'selected_effect_power_scenarios':sel,'observed_power_arithmetic_only':obs_power,'future_GSE73661_scenarios':grid,'GSE73661_MDE_range':[min(r['MDE80'] for r in grid),max(r['MDE80'] for r in grid)],'future_grid_rows':len(future),'two_sided_vs_positive_support':direction,'table_checks':checks}
with (R/'002B_precision_numeric_independent_v1.json').open('x') as f:json.dump(res,f,indent=2);f.write('\n')
print(json.dumps({'status':res['status'],'components':res['components'],'SE_ratio':ratio,'factors':factors,'directional_rule_distinction':direction,'GSE73661_MDE_range':res['GSE73661_MDE_range'],'checks':checks},indent=2))
