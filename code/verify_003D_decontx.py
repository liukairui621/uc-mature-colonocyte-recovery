#!/usr/bin/env python3
import itertools
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t

root=Path(__file__).resolve().parents[1]
run=root/"runs/003D_decontx_state_family"
sp=pd.read_csv(run/"decontx_state_patient_deltas.tsv",sep="\t")
res=pd.read_csv(run/"decontx_fixed15_state_models.tsv",sep="\t")

def exact_p(y,g):
    y=np.asarray(y,float); g=np.asarray(g,int)
    obs=y[g==1].mean()-y[g==0].mean()
    n=len(y); k=int(g.sum()); ge=0; total=0
    allidx=set(range(n))
    for comb in itertools.combinations(range(n),k):
        sel=np.fromiter(comb,dtype=int)
        rest=np.fromiter(allidx-set(comb),dtype=int)
        diff=y[sel].mean()-y[rest].mean()
        ge += abs(diff)>=abs(obs)-1e-12
        total += 1
    return ge/total,ge,total

finite=[]
maxdiff={"estimate":0.0,"lower":0.0,"upper":0.0,"p_exact":0.0}
permutations={}
for _,row in res.iterrows():
    state=row.final_analysis
    d=sp[sp.final_analysis==state].copy()
    g=(d.remission=="Remission").astype(int).to_numpy()
    y=d.delta_state_score.to_numpy(float)
    if len(y)<4 or g.sum()==0 or g.sum()==len(g):
        assert pd.isna(row.estimate) and pd.isna(row.p_exact_permutation)
        continue
    b=y[g==1].mean()-y[g==0].mean()
    df=len(y)-2
    s2=(((y[g==1]-y[g==1].mean())**2).sum()+((y[g==0]-y[g==0].mean())**2).sum())/df
    se=math.sqrt(s2*(1/(g==1).sum()+1/(g==0).sum()))
    lo=b-t.ppf(.975,df)*se
    hi=b+t.ppf(.975,df)*se
    p,ge,total=exact_p(y,g)
    maxdiff["estimate"]=max(maxdiff["estimate"],abs(b-row.estimate))
    maxdiff["lower"]=max(maxdiff["lower"],abs(lo-row.lower))
    maxdiff["upper"]=max(maxdiff["upper"],abs(hi-row.upper))
    maxdiff["p_exact"]=max(maxdiff["p_exact"],abs(p-row.p_exact_permutation))
    finite.append((state,p))
    permutations[state]={"extreme":ge,"total":total,"p":p}
pvals=np.array([p for _,p in finite])
order=np.argsort(pvals)
q=np.empty_like(pvals)
running=1.0
for rankpos in range(len(pvals)-1,-1,-1):
    idx=order[rankpos]
    rank=rankpos+1
    val=pvals[idx]*15/rank
    running=min(running,val)
    q[idx]=min(1.0,running)
qdiff=0.0
for (state,_),qq in zip(finite,q):
    actual=float(res.loc[res.final_analysis==state,"BH_fixed_state_family_n15"].iloc[0])
    qdiff=max(qdiff,abs(qq-actual))
targets={}
for state in ["Non ileal CT colonocyte","Non ileal TA","Non ileal LGR5pos stem"]:
    d=sp[sp.final_analysis==state].sort_values("patient")
    r=res[res.final_analysis==state].iloc[0]
    targets[state]={
        "n":int(len(d)),"remission":int((d.remission=="Remission").sum()),
        "nonremission":int((d.remission!="Remission").sum()),
        "estimate":float(r.estimate),"ci":[float(r.lower),float(r.upper)],
        "p_exact":float(r.p_exact_permutation),"BH_n15":float(r.BH_fixed_state_family_n15),
        "permutation_extreme":permutations[state]["extreme"],
        "permutation_total":permutations[state]["total"],
        "patients":d[["patient","remission","delta_state_score"]].to_dict("records"),
    }
summary={
    "status":"PASS","independent_implementation":"Python grouped means, pooled-t CI, exhaustive label combinations, manual BH with n=15",
    "n_fixed_states":int(len(res)),"n_evaluable_states":int(len(finite)),
    "max_absolute_differences":maxdiff,"max_BH_absolute_difference":qdiff,
    "targets":targets,
}
(root/"runs/003D_decontx_state_family/arithmetic_verification.json").write_text(json.dumps(summary,indent=2,default=lambda o: o.item() if isinstance(o,np.generic) else str(o))+"\n")
print(json.dumps({k:v for k,v in summary.items() if k!="targets"},indent=2))
for k,v in targets.items():
    print(k,{kk:vv for kk,vv in v.items() if kk!="patients"})
