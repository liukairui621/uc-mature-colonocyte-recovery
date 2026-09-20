"""Execute the fixed Stage005 healthy-reference and non-overlapping CT programs."""
from pathlib import Path
import datetime, hashlib, json, platform
import numpy as np
import pandas as pd
from scipy import stats
import scipy, statsmodels, statsmodels.api as sm
from analyze_003F_baseline_adjustment import patient_visits

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/005_health_function'
PLAN=json.loads((ROOT/'planning/analysis_plan_005_health_function.json').read_text())
PRIMARY=PLAN['ct_function']['primary_programs']
CONTEXT=PLAN['ct_function']['context_programs']
def write(df,name): df.to_csv(OUT/(name+'.tsv'),sep='\t',index=False)
def bh(p,n=None):
    p=np.asarray(p,float); finite=np.flatnonzero(np.isfinite(p)); n=len(p) if n is None else n
    out=np.full(len(p),np.nan)
    if len(finite):
        order=finite[np.argsort(p[finite])]; v=p[order]*n/np.arange(1,len(order)+1)
        out[order]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1])
    return out
def mean_ci(x):
    x=np.asarray(x,float); n=len(x); mu=float(x.mean()); se=float(stats.sem(x)) if n>1 else np.nan
    t=stats.t.ppf(.975,n-1) if n>1 else np.nan
    return {'n':n,'mean':mu,'se':se,'lower':mu-t*se,'upper':mu+t*se}
def welch(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);a=x.var(ddof=1)/len(x);b=y.var(ddof=1)/len(y)
    se=np.sqrt(a+b);df=(a+b)**2/(a*a/(len(x)-1)+b*b/(len(y)-1)); est=x.mean()-y.mean();t=stats.t.ppf(.975,df)
    return {'estimate':est,'se':se,'lower':est-t*se,'upper':est+t*se,'p':2*stats.t.sf(abs(est/se),df),'df':df}
def fit(d,model='baseline_ANCOVA'):
    z=d.copy()
    if model=='v31_ANCOVA': z=z[z.library_type.str.lower().eq('3 prime v3.1')]
    cov=['remission_binary']; outcome='delta' if model=='unadjusted_change' else 'post'
    if model!='unadjusted_change': cov+=['baseline']
    if model=='inflammation_adjusted_ANCOVA': cov+=['delta_inflammation_score']
    X=sm.add_constant(z[cov],has_constant='add')
    if len(z)<=X.shape[1]+2 or np.linalg.matrix_rank(X)<X.shape[1]:
        return {'model':model,'n':len(z),'evaluable':False,'p':np.nan,'reason':'rank or residual df'}
    f=sm.OLS(z[outcome],X).fit(); r=f.get_robustcov_results(cov_type='HC3',use_t=True);j=list(X.columns).index('remission_binary');ci=r.conf_int()[j]
    return {'model':model,'n':len(z),'remission':int(z.remission_binary.sum()),'estimate':r.params[j],'se':r.bse[j],'lower':ci[0],'upper':ci[1],'p':r.pvalues[j],'df':f.df_resid,'max_leverage':float(f.get_influence().hat_matrix_diag.max()),'evaluable':True,'reason':''}
def health():
    sc=pd.read_csv(ROOT/'runs/005_bulk_export/sample_scores.tsv',sep='\t')
    reg=pd.read_csv(ROOT/'runs/005_bulk_export/sample_registry.tsv',sep='\t',dtype=str).fillna('')
    reg=reg[reg.disease.isin(['Healthy','non_IBD_control'])].copy()
    reg['control_id']=np.where(reg.subject.ne(''),reg.subject,reg.gsm)
    ctrl=sc.merge(reg[['series','gsm','control_id','title']],on=['series','gsm'],validate='many_to_one')
    ctrl=ctrl[ctrl.program.isin(['CANDIDATE6','COMMON194_INFLAMMATION'])]
    write(ctrl,'healthy_array_scores')
    controls=ctrl.groupby(['series','control_id','program'],as_index=False).agg(score=('score','mean'),n_arrays=('gsm','size'),gsms=('gsm',lambda x:';'.join(x)))
    write(controls,'healthy_control_units')
    p=pd.read_csv(ROOT/'runs/001A_metadata/discovery_pairs_v1.tsv',sep='\t',dtype=str)
    p=p.rename(columns={'post_gsm':'followup_gsm'});p['context']='GSE92415_W6';p['series']='GSE92415';p['outcome']=p.response;p['treatment_path']=p.arm
    q=pd.read_csv(ROOT/'runs/002B_GSE73661_support/support_patient_data.tsv',sep='\t',dtype=str)
    q['context']='GSE73661_IFX_W4W6';q['series']='GSE73661';q['outcome']=q.endoscopic_healing
    v=pd.read_csv(ROOT/'runs/002D_VDZ_support/patient_level_scores.tsv',sep='\t',dtype=str)
    v['context']='GSE73661_VDZ_'+v.visit;v['series']='GSE73661';v['outcome']=v.endoscopic_healing
    fields=['series','context','subject','baseline_gsm','followup_gsm','outcome','treatment_path']
    pairs=pd.concat([p[fields],q[fields],v[fields]],ignore_index=True)
    assert not pairs.duplicated(['context','subject']).any()
    pts=[];groups=[];contrasts=[];paired=[]
    for (context,progseries),pf in pairs.groupby(['context','series'],sort=False):
        for program in ['CANDIDATE6','COMMON194_INFLAMMATION']:
            ss=sc[sc.series.eq(progseries)&sc.program.eq(program)].set_index('gsm').score
            c=controls[controls.series.eq(progseries)&controls.program.eq(program)]
            mu=c.score.mean();sd=c.score.std(ddof=1)
            dd=pf.copy();dd['program']=program;dd['baseline']=ss.loc[dd.baseline_gsm].values;dd['post']=ss.loc[dd.followup_gsm].values;dd['delta']=dd.post-dd.baseline
            dd['control_mean']=mu;dd['control_sd']=sd;dd['baseline_control_centered']=dd.baseline-mu;dd['post_control_centered']=dd.post-mu
            dd['baseline_control_z']=(dd.baseline-mu)/sd;dd['post_control_z']=(dd.post-mu)/sd;pts.append(dd)
            for outcome,gg in dd.groupby('outcome'):
                for visit in ['baseline','post','delta']:
                    groups.append({'context':context,'program':program,'outcome':outcome,'visit':visit,**mean_ci(gg[visit]),'control_mean':mu,'control_sd':sd})
                paired.append({'context':context,'program':program,'outcome':outcome,**mean_ci(gg.delta)})
                for scope in (['main','omit_control7'] if progseries=='GSE73661' else ['main']):
                    cc=c if scope=='main' else c[~c.control_id.isin(['7','7.0'])]
                    contrasts.append({'context':context,'program':program,'outcome':outcome,'scope':scope,'n_patients':len(gg),'n_controls':len(cc),**welch(gg.post,cc.score)})
    pts=pd.concat(pts,ignore_index=True);ct=pd.DataFrame(contrasts);ct['BH_n16']=np.nan
    for scope in ct.scope.unique():
        ix=ct.scope.eq(scope);ct.loc[ix,'BH_n16']=bh(ct.loc[ix,'p'],16)
    assert len(ct[ct.scope.eq('main')])==16
    write(pts,'bulk_patient_scores');write(pd.DataFrame(groups),'bulk_group_levels');write(ct,'healthy_reference_contrasts');write(pd.DataFrame(paired),'paired_change_descriptives')
    strata=pts.groupby(['context','treatment_path','outcome','program']).agg(n=('subject','size'),baseline_mean=('baseline','mean'),post_mean=('post','mean'),delta_mean=('delta','mean')).reset_index()
    write(strata,'bulk_treatment_stratified_descriptives')
    return {'healthy_arrays':ctrl.groupby('series').gsm.nunique().to_dict(),'healthy_units':controls.groupby('series').control_id.nunique().to_dict(),'paired_context_n':pairs.groupby('context').size().to_dict()}
def ct_function():
    counts=pd.read_csv(ROOT/'runs/005_ct_extraction/sample_gene_counts.tsv.gz',sep='\t').set_index('sample_id')
    coverage=pd.read_csv(ROOT/'runs/005_ct_extraction/coverage.tsv',sep='\t').set_index('program')
    pairs=pd.read_csv(ROOT/'runs/003_preexpression/fixed_UC_colonic_site_pairs.tsv',sep='\t')
    pairs=pairs[(counts.loc[pairs.pre_sample_id,'n_cells'].values>=20)&(counts.loc[pairs.post_sample_id,'n_cells'].values>=20)].copy()
    pairs['delta_inflammation_score']=pairs.post_inflammation_score-pairs.pre_inflammation_score
    scores={};genesets={};logexpr=np.log2((counts.drop(columns=['n_cells','total_UMI'])+.5).div(counts.total_UMI+1,axis=0)*1e6)
    for s,gs in PLAN['programs'].items():
        if coverage.loc[s,'coverage']<.8 or coverage.loc[s,'n_measured']<5: continue
        genesets[s]=[g for g in gs if g in logexpr];scores[s]=logexpr[genesets[s]].mean(axis=1)
    def patient_data(score,s):
        d=pairs.copy();d['baseline']=score.loc[d.pre_sample_id].values;d['post']=score.loc[d.post_sample_id].values;d['delta']=d.post-d.baseline
        pt=d.groupby('patient',as_index=False).agg(remission=('remission','first'),library_type=('library_type','first'),n_site_pairs=('site','size'),baseline=('baseline','mean'),post=('post','mean'),delta=('delta','mean'),delta_inflammation_score=('delta_inflammation_score','mean'))
        pt['remission_binary']=pt.remission.eq('Remission').astype(int);pt['program']=s
        return pt
    pts=pd.concat([patient_data(s,x) for x,s in scores.items()],ignore_index=True)
    c=pts[pts.program.eq('CANDIDATE6')].set_index('patient')
    old=patient_visits(pd.read_csv(ROOT/'runs/003_cell_context/site_pair_metrics.tsv',sep='\t'),'all_fixed').set_index('patient')
    assert set(c.index)==set(old.index) and len(c)==16 and c.remission_binary.sum()==6
    err=max(float(np.max(np.abs(c.loc[old.index,new]-old[prior]))) for new,prior in [('baseline','pre_ct_state_score'),('post','post_ct_state_score'),('delta','delta_ct_state_score')])
    assert err<1e-10,err
    write(pts,'ct_patient_scores');write(pairs,'ct_eligible_site_pairs');write(logexpr.reset_index(),'ct_gene_logCPM')
    models=[];levels=[];correlations=[];genes=[];lgo=[]
    for s in PLAN['programs']:
        if s not in scores:
            models.append({'program':s,'model':'baseline_ANCOVA','evaluable':False,'reason':'predefined feature coverage failed','p':np.nan});continue
        d=pts[pts.program.eq(s)]
        for m in ['baseline_ANCOVA','unadjusted_change','v31_ANCOVA','inflammation_adjusted_ANCOVA']:
            models.append({'program':s,**fit(d,m)})
        for outcome,g in d.groupby('remission'):
            for visit in ['baseline','post','delta']:
                levels.append({'program':s,'remission':outcome,'visit':visit,**mean_ci(g[visit])})
        if s in PRIMARY:
            z=d.set_index('patient').loc[c.index]
            r,p=stats.pearsonr(z.delta,c.delta);zz=np.arctanh(np.clip(r,-.999999,.999999));width=1.96/np.sqrt(len(z)-3)
            correlations.append({'program':s,'n':len(z),'pearson_r':r,'lower':np.tanh(zz-width),'upper':np.tanh(zz+width),'p':p})
            for gene in genesets[s]:
                g=patient_data(logexpr[gene],s)
                yes=g[g.remission_binary.eq(1)];no=g[g.remission_binary.eq(0)]
                genes.append({'program':s,'gene':gene,'mean_delta_remission':yes.delta.mean(),'mean_delta_nonremission':no.delta.mean(),'change_contrast':yes.delta.mean()-no.delta.mean(),'mean_logCPM':logexpr.loc[list(set(pairs.pre_sample_id)|set(pairs.post_sample_id)),gene].mean()})
                other=[k for k in genesets[s] if k!=gene]
                g2=patient_data(logexpr[other].mean(axis=1),s);rr=fit(g2)
                lgo.append({'program':s,'omitted_gene':gene,**rr})
    models=pd.DataFrame(models);models['BH_primary_n3']=np.nan;models['BH_context_n2']=np.nan
    for m in models.model.unique():
        ix=models.program.isin(PRIMARY)&models.model.eq(m);models.loc[ix,'BH_primary_n3']=bh(models.loc[ix,'p'],3)
        ix=models.program.isin(CONTEXT)&models.model.eq(m);models.loc[ix,'BH_context_n2']=bh(models.loc[ix,'p'],2)
    cor=pd.DataFrame(correlations);cor['BH_n3']=bh(cor.p,3)
    write(models,'ct_program_models');write(pd.DataFrame(levels),'ct_group_levels');write(cor,'ct_candidate_coordination');write(pd.DataFrame(genes),'ct_all_member_gene_changes');write(pd.DataFrame(lgo),'ct_leave_one_gene_out')
    det=pd.read_csv(ROOT/'runs/005_ct_extraction/sample_gene_detected_cells.tsv.gz',sep='\t').set_index('sample_id')
    ids=sorted(set(pairs.pre_sample_id)|set(pairs.post_sample_id));qc=[]
    for s,gs in genesets.items():
        for g in gs:
            qc.append({'program':s,'gene':g,'fraction_CT_cells_detected':det.loc[ids,g].sum()/counts.loc[ids,'n_cells'].sum(),'n_eligible_samples_detected':int((det.loc[ids,g]>0).sum()),'n_eligible_samples':len(ids)})
    write(pd.DataFrame(qc),'ct_program_detection')
    return {'candidate_score_max_error':err,'patients':len(c),'remission_patients':int(c.remission_binary.sum()),'site_pairs':len(pairs),'primary_programs':models[models.model.eq('baseline_ANCOVA')&models.program.isin(PRIMARY)].to_dict('records')}
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'analysis_summary.json').exists(): raise RuntimeError('Completed Stage005 exists; refuse overwrite')
    summary={'status':'COMPLETE','role':'post-result exploratory extension','health':health(),'ct_function':ct_function(),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'environment':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__},'plan_sha256':hashlib.sha256((ROOT/'planning/analysis_plan_005_health_function.json').read_bytes()).hexdigest()}
    (OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
