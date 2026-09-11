#!/usr/bin/env python3
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/"runs/003D_decontx_batches/batch_manifest.tsv"
VAR=ROOT/"inputs/002B_asset_audit/TAURUS_actual_var_features_full.tsv"
PAIRS=ROOT/"runs/003_preexpression/fixed_UC_colonic_site_pairs.tsv"
DECONTX_AGG=ROOT/"runs/003D_decontx/sample_state_decontx_all.tsv"
RAW_MODELS=ROOT/"runs/003_cell_context/patient_group_models.tsv"
OUT=ROOT/"runs/003E_ambient_negative_controls"
OUT.mkdir(parents=True,exist_ok=True)

candidate=["AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12"]
negative=["PTPRC","COL1A1"]
genes=candidate+negative

def exact_p(y,g):
    y=np.asarray(y,float); g=np.asarray(g,int)
    n=len(y); k=int(g.sum())
    if n<4 or k==0 or k==n: return np.nan
    obs=y[g==1].mean()-y[g==0].mean()
    extreme=0; total=0
    allidx=np.arange(n)
    for comb in itertools.combinations(range(n),k):
        sel=np.fromiter(comb,dtype=int)
        mask=np.ones(n,dtype=bool); mask[sel]=False
        d=y[sel].mean()-y[allidx[mask]].mean()
        extreme += abs(d)>=abs(obs)-1e-12
        total += 1
    return extreme/total

def fit_group(d,value):
    d=d[np.isfinite(d[value])].copy()
    g=(d["remission"]=="Remission").astype(int).to_numpy()
    y=d[value].to_numpy(float)
    n=len(y); k=int(g.sum())
    base={"metric":value,"n":n,"remission":k,"nonremission":n-k}
    if n<4 or k==0 or k==n:
        return {**base,**{z:np.nan for z in ["mean_remission","mean_nonremission","estimate","lower","upper","p_OLS","p_HC3","p_exact_permutation"]}}
    X=np.column_stack([np.ones(n),g])
    inv=np.linalg.inv(X.T@X)
    beta=inv@X.T@y
    e=y-X@beta
    df=n-2
    s2=(e@e)/df
    se=math.sqrt(s2*inv[1,1])
    crit=stats.t.ppf(.975,df)
    h=np.einsum("ij,jk,ik->i",X,inv,X)
    meat=X.T@((e**2/(1-h)**2)[:,None]*X)
    hc3=math.sqrt((inv@meat@inv)[1,1])
    return {**base,
      "mean_remission":float(y[g==1].mean()),"mean_nonremission":float(y[g==0].mean()),
      "estimate":float(beta[1]),"lower":float(beta[1]-crit*se),"upper":float(beta[1]+crit*se),
      "p_OLS":float(2*stats.t.sf(abs(beta[1]/se),df)),
      "p_HC3":float(2*stats.t.sf(abs(beta[1]/hc3),df)),
      "p_exact_permutation":float(exact_p(y,g))}

def bh_fixed(vals,n):
    vals=np.asarray(vals,float); ok=np.isfinite(vals); out=np.full(len(vals),np.nan)
    idx=np.where(ok)[0]; order=idx[np.argsort(vals[idx])]
    running=1.0
    for pos in range(len(order)-1,-1,-1):
        ix=order[pos]; rank=pos+1
        running=min(running,vals[ix]*n/rank)
        out[ix]=min(1.0,running)
    return out

manifest=pd.read_csv(MANIFEST,sep="\t").sort_values("order")
var=pd.read_csv(VAR,sep="\t")
indices={g:int(np.flatnonzero(var["gene_symbol"].eq(g))[0]) for g in genes}
if len(indices)!=len(genes): raise AssertionError("gene mapping failure")
rows=[]
for rec in manifest.to_dict("records"):
    x=sparse.load_npz(ROOT/rec["matrix_path"]).tocsr()
    md=pd.read_csv(ROOT/rec["metadata_path"],sep="\t")
    if x.shape[0]!=len(md) or x.shape[1]!=len(var): raise AssertionError(rec["sample_id"])
    mask=md["final_analysis"].eq("Non ileal CT colonocyte").to_numpy()
    xs=x[mask,:]
    row={"sample_id":rec["sample_id"],"n_cells":int(mask.sum()),
         "raw_CT_UMI":float(xs.data.astype(np.float64).sum(dtype=np.float64))}
    for g,ix in indices.items():
        row[f"{g}_counts"]=float(xs[:,ix].data.astype(np.float64).sum(dtype=np.float64))
    if mask.any():
        one=md.loc[mask].iloc[0]
        for col in ["Patient","Site","Treatment","Remission_status","LibraryType","Batch"]:
            row[col]=one[col]
    else:
        one=md.iloc[0]
        for col in ["Patient","Site","Treatment","Remission_status","LibraryType","Batch"]:
            row[col]=one[col]
    rows.append(row)
sample=pd.DataFrame(rows)
for g in genes:
    sample[f"{g}_logCPM"]=np.log2((sample[f"{g}_counts"]+.5)/(sample["raw_CT_UMI"]+1)*1e6)
sample["Candidate6_score"]=sample[[f"{g}_logCPM" for g in candidate]].mean(axis=1)
sample["NegativeControl2_score"]=sample[[f"{g}_logCPM" for g in negative]].mean(axis=1)
sample.to_csv(OUT/"raw_CT_sample_negative_control_counts.tsv",sep="\t",index=False)

pairs=pd.read_csv(PAIRS,sep="\t")
metrics=[f"{g}_logCPM" for g in genes]+["Candidate6_score","NegativeControl2_score"]
pre=sample[["sample_id","n_cells"]+metrics].copy()
pre.columns=["pre_sample_id","pre_n"]+[f"pre_{m}" for m in metrics]
post=sample[["sample_id","n_cells"]+metrics].copy()
post.columns=["post_sample_id","post_n"]+[f"post_{m}" for m in metrics]
site=pairs.merge(pre,on="pre_sample_id",how="left").merge(post,on="post_sample_id",how="left")
site=site[(site.pre_n>=20)&(site.post_n>=20)].copy()
for m in metrics: site[f"delta_{m}"]=site[f"post_{m}"]-site[f"pre_{m}"]
delta_cols=[f"delta_{m}" for m in metrics]
site.to_csv(OUT/"negative_control_CT_site_deltas.tsv",sep="\t",index=False)
patient=site.groupby(["patient","remission"],as_index=False)[delta_cols].mean()
patient.to_csv(OUT/"negative_control_CT_patient_deltas.tsv",sep="\t",index=False)
model_rows=[fit_group(patient,m) for m in delta_cols]
models=pd.DataFrame(model_rows)
negmask=models.metric.isin([f"delta_{g}_logCPM" for g in negative])
models["BH_negative_genes_n2"]=np.nan
models.loc[negmask,"BH_negative_genes_n2"]=bh_fixed(models.loc[negmask,"p_exact_permutation"].to_numpy(),2)
models.to_csv(OUT/"negative_control_CT_models.tsv",sep="\t",index=False)

raw=pd.read_csv(RAW_MODELS,sep="\t")
ref=raw[(raw.branch=="exclude_predicted_doublets::all_fixed::threshold20")&(raw.metric=="delta_ct_state_score")].iloc[0]
cand=models[models.metric=="delta_Candidate6_score"].iloc[0]
positive_check={
 "expected_estimate":float(ref.estimate),"recomputed_estimate":float(cand.estimate),
 "absolute_difference":abs(float(ref.estimate)-float(cand.estimate)),
 "expected_lower":float(ref.lower),"recomputed_lower":float(cand.lower),
 "expected_upper":float(ref.upper),"recomputed_upper":float(cand.upper),
 "expected_exact_p":float(ref.p_exact_permutation),"recomputed_exact_p":float(cand.p_exact_permutation),
}
if max(positive_check["absolute_difference"],abs(positive_check["expected_lower"]-positive_check["recomputed_lower"]),abs(positive_check["expected_upper"]-positive_check["recomputed_upper"]),abs(positive_check["expected_exact_p"]-positive_check["recomputed_exact_p"]))>1e-10:
    raise AssertionError(positive_check)

dc=pd.read_csv(DECONTX_AGG,sep="\t")
ct=dc[dc.final_analysis=="Non ileal CT colonocyte"].copy()
contmetrics=["contamination_mean","contamination_median","contamination_q90"]
cross=[]
for m in contmetrics:
    a=ct.loc[ct.Remission_status=="Remission",m].to_numpy(float)
    b=ct.loc[ct.Remission_status=="Non_Remission",m].to_numpy(float)
    w=stats.mannwhitneyu(a,b,alternative="two-sided",method="asymptotic",use_continuity=True)
    cross.append({"metric":m,"unit":"sample_state","n_remission":len(a),"n_nonremission":len(b),
                  "mean_remission":a.mean(),"mean_nonremission":b.mean(),
                  "median_remission":np.median(a),"median_nonremission":np.median(b),
                  "wilcoxon_rank_sum_p":float(w.pvalue)})
cross=pd.DataFrame(cross)
cross.to_csv(OUT/"CT_contamination_cross_sectional_sample_state.tsv",sep="\t",index=False)

cp=ct[["sample_id","n_cells"]+contmetrics].copy()
pre=cp.copy(); pre.columns=["pre_sample_id","pre_n"]+[f"pre_{m}" for m in contmetrics]
post=cp.copy(); post.columns=["post_sample_id","post_n"]+[f"post_{m}" for m in contmetrics]
cs=pairs.merge(pre,on="pre_sample_id",how="left").merge(post,on="post_sample_id",how="left")
cs=cs[(cs.pre_n>=20)&(cs.post_n>=20)].copy()
for m in contmetrics: cs[f"delta_{m}"]=cs[f"post_{m}"]-cs[f"pre_{m}"]
cs.to_csv(OUT/"CT_contamination_site_deltas.tsv",sep="\t",index=False)
cdelta=[f"delta_{m}" for m in contmetrics]
cp_patient=cs.groupby(["patient","remission"],as_index=False)[cdelta].mean()
cp_patient.to_csv(OUT/"CT_contamination_patient_deltas.tsv",sep="\t",index=False)
cmodels=pd.DataFrame([fit_group(cp_patient,m) for m in cdelta])
cmodels.to_csv(OUT/"CT_contamination_patient_models.tsv",sep="\t",index=False)

summary={
 "status":"COMPLETE","role":"post-result ambient diagnostic; no confirmatory tests",
 "created_utc":datetime.now(timezone.utc).isoformat(),
 "negative_controls":negative,"candidate_positive_control":candidate,
 "state":"Non ileal CT colonocyte","threshold":20,
 "site_pairs":int(site.shape[0]),"patients":int(patient.shape[0]),
 "remission_patients":int((patient.remission=="Remission").sum()),
 "nonremission_patients":int((patient.remission=="Non_Remission").sum()),
 "positive_control_reproduction":positive_check,
 "negative_control_models":models[negmask|models.metric.eq("delta_NegativeControl2_score")].to_dict("records"),
 "contamination_cross_sectional":cross.to_dict("records"),
 "contamination_patient_models":cmodels.to_dict("records"),
 "interpretation_limits":[
  "A null negative-control result weighs against a broad immune/stromal ambient shift aligned with remission but does not prove absence of gene-specific ambient RNA.",
  "Contamination estimates remain DecontX-model-dependent.",
  "Sample-state cross-sectional Wilcoxon rows are repeated within patients; patient-level Post-minus-Pre models are the aligned analysis."
 ]
}
def json_clean(obj):
    if isinstance(obj,dict): return {k:json_clean(v) for k,v in obj.items()}
    if isinstance(obj,list): return [json_clean(v) for v in obj]
    if isinstance(obj,(float,np.floating)) and not np.isfinite(obj): return None
    if isinstance(obj,np.generic): return obj.item()
    return obj
summary=json_clean(summary)
(OUT/"analysis_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
print(json.dumps(summary,indent=2,allow_nan=False))
