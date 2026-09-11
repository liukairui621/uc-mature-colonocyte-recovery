#!/usr/bin/env python3
import json
from pathlib import Path
import pandas as pd
import numpy as np
import scipy.sparse as sparse

root=Path(__file__).resolve().parents[1]
run=root/"runs/003D_decontx"
out=run/"sample_results"
manifest=pd.read_csv(root/"runs/003D_decontx_batches/batch_manifest.tsv",sep="\t")
genes=["AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12"]
status_rows=[]
tables=[]
checks=[]
for rec in manifest.to_dict("records"):
    sid=rec["sample_id"]
    sp=out/f"{sid}.status.json"
    tp=out/f"{sid}.sample_state_decontx.tsv"
    if not sp.exists() or not tp.exists():
        raise FileNotFoundError(f"missing output for {sid}")
    status=json.loads(sp.read_text())
    if status.get("status")!="COMPLETE":
        raise RuntimeError(f"non-complete status for {sid}: {status}")
    tab=pd.read_csv(tp,sep="\t")
    if set(tab["sample_id"])!={sid}:
        raise AssertionError(f"sample id mismatch {sid}")
    n_cells=int(tab["n_cells"].sum())
    corrected=float(tab["decontx_total_UMI"].sum())
    raw=float(status["raw_total_UMI"])
    precise_raw=float(sparse.load_npz(root/rec["matrix_path"]).data.astype(np.float64).sum(dtype=np.float64))
    if n_cells!=int(rec["n_cells"]) or n_cells!=int(status["n_cells"]):
        raise AssertionError(f"cell mismatch {sid}")
    if abs(raw-precise_raw)>1e-6:
        raise AssertionError(f"R versus float64 raw UMI mismatch {sid}: {raw} vs {precise_raw}")
    manifest_raw_delta=raw-float(rec["raw_total_UMI"])
    if abs(manifest_raw_delta)>max(10,1e-6*raw):
        raise AssertionError(f"manifest float32 raw UMI mismatch too large {sid}: {manifest_raw_delta}")
    if abs(corrected-float(status["corrected_total_UMI"]))>max(1e-5,1e-10*corrected):
        raise AssertionError(f"corrected UMI mismatch {sid}")
    if corrected<0 or corrected>raw+1e-6:
        raise AssertionError(f"corrected UMI outside [0,raw] {sid}")
    gene_cols=[f"{g}_decontx_counts" for g in genes]
    vals=tab[gene_cols].to_numpy(float)
    if not np.isfinite(vals).all() or (vals<0).any():
        raise AssertionError(f"invalid gene counts {sid}")
    if (tab[gene_cols].sum(axis=1)>tab["decontx_total_UMI"]+1e-6).any():
        raise AssertionError(f"candidate sum exceeds total {sid}")
    status_rows.append({
        "sample_id":sid,"order":int(rec["order"]),"n_cells":n_cells,
        "n_genes":int(status["n_genes"]),"n_clusters":int(status["n_clusters"]),
        "seed":int(status["seed"]),"celda_version":status["celda_version"],
        "matrix_version":status["matrix_version"],
        "contamination_mean":float(status["contamination_mean"]),
        "contamination_median":float(status["contamination_median"]),
        "contamination_q90":float(status["contamination_q90"]),
        "raw_total_UMI":raw,"manifest_float32_raw_UMI":float(rec["raw_total_UMI"]),
        "manifest_float32_raw_delta":manifest_raw_delta,"corrected_total_UMI":corrected,
        "retained_fraction":corrected/raw,
        "status_path":str(sp.relative_to(root)),
        "state_path":str(tp.relative_to(root)),
    })
    tables.append(tab)
    checks.append({"sample_id":sid,"cell_count_check":"PASS","raw_UMI_check":"PASS",
                   "corrected_UMI_check":"PASS","candidate_bounds_check":"PASS"})
alltab=pd.concat(tables,ignore_index=True)
status_df=pd.DataFrame(status_rows).sort_values("order")
checks_df=pd.DataFrame(checks)
alltab.to_csv(run/"sample_state_decontx_all.tsv",sep="\t",index=False)
status_df.to_csv(run/"sample_decontx_status.tsv",sep="\t",index=False)
checks_df.to_csv(run/"sample_integrity_checks.tsv",sep="\t",index=False)
state_rows=[]
for state,d in alltab.groupby("final_analysis",sort=True):
    w=d["n_cells"].to_numpy(float)
    state_rows.append({
        "final_analysis":state,"major":d["major"].iloc[0],
        "n_samples":d["sample_id"].nunique(),"n_cells":int(w.sum()),
        "cell_weighted_contamination_mean":float(np.average(d["contamination_mean"],weights=w)),
        "sample_median_contamination_mean":float(d["contamination_mean"].median()),
        "sample_min_contamination_mean":float(d["contamination_mean"].min()),
        "sample_max_contamination_mean":float(d["contamination_mean"].max()),
    })
state_df=pd.DataFrame(state_rows)
state_df.to_csv(run/"contamination_by_state.tsv",sep="\t",index=False)
summary={
    "status":"COMPLETE","n_samples":int(len(status_df)),
    "n_cells":int(status_df["n_cells"].sum()),"n_genes":int(status_df["n_genes"].nunique()==1 and status_df["n_genes"].iloc[0]),
    "n_sample_state_rows":int(len(alltab)),"n_states_all":int(alltab["final_analysis"].nunique()),
    "n_non_ileal_epithelial_states":int(alltab.loc[alltab["major"]=="Non_ileal_epithelium","final_analysis"].nunique()),
    "all_status_complete":bool((status_df.shape[0]==manifest.shape[0])),
    "all_integrity_checks_pass":True,
    "cell_weighted_contamination_mean":float(np.average(status_df["contamination_mean"],weights=status_df["n_cells"])),
    "sample_median_contamination_mean":float(status_df["contamination_mean"].median()),
    "sample_contamination_mean_range":[float(status_df["contamination_mean"].min()),float(status_df["contamination_mean"].max())],
    "raw_total_UMI":float(status_df["raw_total_UMI"].sum()),
    "corrected_total_UMI":float(status_df["corrected_total_UMI"].sum()),
    "overall_retained_fraction":float(status_df["corrected_total_UMI"].sum()/status_df["raw_total_UMI"].sum()),
    "celda_versions":sorted(status_df["celda_version"].unique().tolist()),
    "matrix_versions":sorted(status_df["matrix_version"].unique().tolist()),
    "empty_droplet_background_used":False,
    "cluster_labels":"released final_analysis",
    "ambient_pool":"sample library",
}
(run/"consolidation_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
