#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
src=root/"runs/003D_decontx_state_family/decontx_fixed15_state_models.tsv"
outdir=root/"figures"
outdir.mkdir(exist_ok=True)
x=pd.read_csv(src,sep="\t")
x=x[np.isfinite(x["estimate"])].copy()
x=x.sort_values("estimate")
labels=x["final_analysis"].str.replace("Non ileal ","",regex=False)
y=np.arange(len(x))
fig,ax=plt.subplots(figsize=(8.4,6.4))
ax.axvline(0,color="#777777",lw=1,ls="--",zorder=0)
ax.errorbar(x["raw_estimate"],y-0.13,
            xerr=np.vstack([x["raw_estimate"]-x["raw_lower"],x["raw_upper"]-x["raw_estimate"]]),
            fmt="s",ms=5,color="#0072B2",ecolor="#0072B2",capsize=2.5,lw=1.2,label="Raw counts")
ax.errorbar(x["estimate"],y+0.13,
            xerr=np.vstack([x["estimate"]-x["lower"],x["upper"]-x["estimate"]]),
            fmt="o",ms=5,color="#D55E00",ecolor="#D55E00",capsize=2.5,lw=1.2,label="DecontX corrected")
for i,row in x.reset_index(drop=True).iterrows():
    if row["final_analysis"] in ["Non ileal CT colonocyte","Non ileal TA","Non ileal LGR5pos stem"]:
        ax.text(max(row["upper"],row["raw_upper"])+0.08,i+0.13,f"q={row['BH_fixed_state_family_n15']:.3f}",va="center",fontsize=8,color="#D55E00")
ax.set_yticks(y)
ax.set_yticklabels(labels,fontsize=8.5)
ax.set_xlabel("Remission minus non-remission change in six-gene state score\n(estimate and OLS 95% CI)")
ax.set_title("Fixed 15-state family: raw-count and DecontX sensitivity",loc="left",fontweight="bold")
ax.legend(frameon=False,loc="lower right")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.grid(axis="x",color="#dddddd",lw=.6)
fig.tight_layout()
fig.savefig(outdir/"Stage003D_DecontX_raw_vs_corrected.png",dpi=300,bbox_inches="tight")
fig.savefig(outdir/"Stage003D_DecontX_raw_vs_corrected.pdf",bbox_inches="tight")
x.to_csv(outdir/"Stage003D_DecontX_raw_vs_corrected.source.tsv",sep="\t",index=False)
print(outdir/"Stage003D_DecontX_raw_vs_corrected.png")
