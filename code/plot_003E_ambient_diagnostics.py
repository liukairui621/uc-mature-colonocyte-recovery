#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
run=root/"runs/003E_ambient_negative_controls"
figdir=root/"figures"; figdir.mkdir(exist_ok=True)
m=pd.read_csv(run/"negative_control_CT_models.tsv",sep="\t")
c=pd.read_csv(run/"CT_contamination_patient_models.tsv",sep="\t")
sel=m[m.metric.isin(["delta_Candidate6_score","delta_NegativeControl2_score","delta_PTPRC_logCPM","delta_COL1A1_logCPM"])].copy()
order=["delta_COL1A1_logCPM","delta_PTPRC_logCPM","delta_NegativeControl2_score","delta_Candidate6_score"]
sel=sel.set_index("metric").loc[order].reset_index()
sel["label"]=["COL1A1","PTPRC","Negative-control mean","Candidate6"]
colors=["#777777","#777777","#777777","#0072B2"]
fig,axes=plt.subplots(1,2,figsize=(10.5,4.8),gridspec_kw={"width_ratios":[1.2,1]})
ax=axes[0]; y=np.arange(len(sel))
ax.axvline(0,color="#777777",ls="--",lw=1)
for i,r in sel.iterrows():
 ax.errorbar(r.estimate,i,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt="o",color=colors[i],ecolor=colors[i],capsize=3,lw=1.4)
 ax.text(r.upper+.12,i,f"exact P={r.p_exact_permutation:.4f}",va="center",fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(sel.label)
ax.set_xlabel("Remission minus non-remission\nPost-minus-Pre logCPM change")
ax.set_title("A  Raw CT-pseudobulk controls",loc="left",fontweight="bold")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.grid(axis="x",color="#dddddd",lw=.6)
ax=axes[1]; y=np.arange(len(c))
labels=["Mean","Median","90th percentile"]
ax.axvline(0,color="#777777",ls="--",lw=1)
for i,r in c.iterrows():
 ax.errorbar(r.estimate,i,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt="o",color="#D55E00",ecolor="#D55E00",capsize=3,lw=1.4)
 ax.text(r.upper+.015,i,f"exact P={r.p_exact_permutation:.4f}",va="center",fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(labels)
ax.set_xlabel("Remission minus non-remission\nPost-minus-Pre contamination change")
ax.set_title("B  DecontX CT contamination",loc="left",fontweight="bold")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.grid(axis="x",color="#dddddd",lw=.6)
fig.suptitle("Ambient-RNA diagnostics use the patient as the statistical unit",fontweight="bold",y=1.02)
fig.tight_layout()
out=figdir/"Stage003E_ambient_negative_control_diagnostics"
fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches="tight")
fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight")
pd.concat([sel.assign(panel="negative_control"),c.assign(panel="contamination")],ignore_index=True,sort=False).to_csv(out.with_suffix(".source.tsv"),sep="\t",index=False)
print(out)
