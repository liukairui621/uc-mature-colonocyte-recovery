from pathlib import Path
import argparse, hashlib, json, datetime, platform
import h5py, numpy as np, pandas as pd
from scipy.sparse import csr_matrix

CANDIDATES=["AQP8","HMGCS2","GUCA2A","CA2","SLC26A3","MS4A12"]

def decode_array(a):
    a=np.asarray(a)
    if a.dtype.kind in "SO":
        return np.array([v.decode() if isinstance(v,(bytes,bytearray)) else str(v) for v in a],dtype=object)
    return a

def read_frame_column(h5, frame, key):
    obj=h5[f"{frame}/{key}"]
    if isinstance(obj,h5py.Group) and "codes" in obj and "categories" in obj:
        cats=decode_array(obj["categories"][:])
        codes=np.asarray(obj["codes"][:],dtype=int)
        out=np.empty(len(codes),dtype=object)
        out[:]=None
        ok=codes>=0
        out[ok]=cats[codes[ok]]
        return out
    return decode_array(obj[:])

def hash_file(path,block=16*1024*1024):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(block),b""): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True)
    ap.add_argument("--pairs",required=True)
    ap.add_argument("--outdir",required=True)
    ap.add_argument("--chunk-cells",type=int,default=25000)
    ap.add_argument("--smoke-cells",type=int,default=0,
                    help="If positive, process only this many rows and write TESTED outputs.")
    args=ap.parse_args()
    inp=Path(args.input); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    pairs=pd.read_csv(args.pairs,sep="\t")
    selected=set(pairs.pre_sample_id)|set(pairs.post_sample_id)
    var=pd.read_csv(inp.parent.parent/"002B_asset_audit/TAURUS_actual_var_features_full.tsv",sep="\t")
    hits={g:var.index[var.gene_symbol.eq(g)].tolist() for g in CANDIDATES}
    if any(len(v)!=1 for v in hits.values()): raise RuntimeError(f"Candidate mapping not unique: {hits}")
    cand_idx=np.array([hits[g][0] for g in CANDIDATES],dtype=int)
    aggregates=[]
    qc=[]
    with h5py.File(inp,"r") as h:
        if "X" not in h or not isinstance(h["X"],h5py.Group):
            raise RuntimeError("Expected sparse X group")
        xg=h["X"]
        enc=xg.attrs.get("encoding-type","")
        if isinstance(enc,bytes): enc=enc.decode()
        if enc not in ("csr_matrix","csr"):
            raise RuntimeError(f"Expected cell-row CSR X, got {enc}")
        n_cells,n_genes=map(int,xg.attrs["shape"])
        obs={
          k:read_frame_column(h,"obs",k) for k in
          ["sample_id","Patient","Disease","Site","Treatment","Remission_status",
           "final_analysis","major","minor","predicted_doublets","total_counts",
           "n_genes_by_counts","pct_counts_mt","LibraryType","Batch"]
        }
        if any(len(v)!=n_cells for v in obs.values()): raise RuntimeError("obs length mismatch")
        n_run=min(n_cells,args.smoke_cells) if args.smoke_cells else n_cells
        max_rowsum_diff=0.0; mismatch_rows=0; noninteger_values=0; selected_cells=0
        for start in range(0,n_run,args.chunk_cells):
            end=min(start+args.chunk_cells,n_run)
            ip=np.asarray(xg["indptr"][start:end+1],dtype=np.int64)
            lo,hi=int(ip[0]),int(ip[-1])
            dat=np.asarray(xg["data"][lo:hi])
            idx=np.asarray(xg["indices"][lo:hi],dtype=np.int32)
            mat=csr_matrix((dat,idx,ip-lo),shape=(end-start,n_genes))
            noninteger_values += int(np.count_nonzero(np.abs(dat-np.rint(dat))>1e-7))
            rowsum=np.asarray(mat.sum(axis=1)).ravel()
            released=np.asarray(obs["total_counts"][start:end],dtype=float)
            dif=np.abs(rowsum-released)
            max_rowsum_diff=max(max_rowsum_diff,float(dif.max(initial=0)))
            mismatch_rows += int(np.count_nonzero(dif>1e-4))
            keep=np.array([(str(s) in selected) and str(d)=="UC" and
                           str(si)!="Terminal_Ileum" and fa is not None
                           for s,d,si,fa in zip(obs["sample_id"][start:end],
                                                obs["Disease"][start:end],
                                                obs["Site"][start:end],
                                                obs["final_analysis"][start:end])])
            if not keep.any(): continue
            loc=np.flatnonzero(keep)
            sub=mat[loc,:][:,cand_idx].toarray().astype(float)
            total=rowsum[loc]
            norm=sub/np.maximum(total[:,None],1.0)*10000.0
            frame=pd.DataFrame({
              "sample_id":np.asarray(obs["sample_id"][start:end],object)[loc],
              "Patient":np.asarray(obs["Patient"][start:end],object)[loc],
              "Site":np.asarray(obs["Site"][start:end],object)[loc],
              "Treatment":np.asarray(obs["Treatment"][start:end],object)[loc],
              "Remission_status":np.asarray(obs["Remission_status"][start:end],object)[loc],
              "final_analysis":np.asarray(obs["final_analysis"][start:end],object)[loc],
              "major":np.asarray(obs["major"][start:end],object)[loc],
              "minor":np.asarray(obs["minor"][start:end],object)[loc],
              "predicted_doublets":np.asarray(obs["predicted_doublets"][start:end])[loc].astype(bool),
              "total_UMI":total,
              "n_genes_by_counts":np.asarray(obs["n_genes_by_counts"][start:end])[loc],
              "pct_counts_mt":np.asarray(obs["pct_counts_mt"][start:end])[loc],
              "LibraryType":np.asarray(obs["LibraryType"][start:end],object)[loc],
              "Batch":np.asarray(obs["Batch"][start:end],object)[loc]
            })
            for j,g in enumerate(CANDIDATES):
                frame[f"{g}_counts"]=sub[:,j]
                frame[f"{g}_norm10k"]=norm[:,j]
            for cellset,fr in [("author_annotated",frame),
                               ("exclude_predicted_doublets",frame.loc[~frame.predicted_doublets])]:
                if fr.empty: continue
                keys=["sample_id","Patient","Site","Treatment","Remission_status",
                      "LibraryType","Batch","final_analysis","major","minor"]
                named={f"{g}_counts":(f"{g}_counts","sum") for g in CANDIDATES}
                named.update({f"{g}_norm10k_sum":(f"{g}_norm10k","sum") for g in CANDIDATES})
                ag=(fr.groupby(keys,dropna=False,observed=True)
                      .agg(n_cells=("total_UMI","size"),total_UMI=("total_UMI","sum"),
                           mean_n_genes=("n_genes_by_counts","mean"),
                           mean_pct_mt=("pct_counts_mt","mean"),**named).reset_index())
                ag.insert(0,"cell_set",cellset)
                aggregates.append(ag)
            q=(frame.groupby(["sample_id","Patient","Site","Treatment","Remission_status",
                              "LibraryType","Batch"],observed=True)
                 .agg(author_annotated_cells=("total_UMI","size"),
                      predicted_doublets=("predicted_doublets","sum"),
                      total_UMI=("total_UMI","sum"),
                      median_n_genes=("n_genes_by_counts","median"),
                      median_pct_mt=("pct_counts_mt","median")).reset_index())
            qc.append(q)
            selected_cells += len(frame)
            print(json.dumps({"chunk_end":end,"n_run":n_run,"selected_cells_cumulative":selected_cells}),flush=True)
    agg=pd.concat(aggregates,ignore_index=True) if aggregates else pd.DataFrame()
    qcdf=pd.concat(qc,ignore_index=True)
    if not qcdf.empty:
        qcdf=(qcdf.groupby(["sample_id","Patient","Site","Treatment","Remission_status",
                            "LibraryType","Batch"],as_index=False,observed=True)
                  .agg(author_annotated_cells=("author_annotated_cells","sum"),
                       predicted_doublets=("predicted_doublets","sum"),
                       total_UMI=("total_UMI","sum"),
                       median_n_genes=("median_n_genes","median"),
                       median_pct_mt=("median_pct_mt","median")))
    prefix="smoke_" if args.smoke_cells else ""
    agg.to_csv(out/f"{prefix}sample_state_candidate_aggregates.tsv.gz",sep="\t",index=False,
               compression="gzip")
    qcdf.to_csv(out/f"{prefix}sample_qc_aggregates.tsv",sep="\t",index=False)
    summary={
      "status":"TESTED" if args.smoke_cells else "COMPLETE",
      "input":str(inp),"input_bytes":inp.stat().st_size,
      "input_sha256":None if args.smoke_cells else hash_file(inp),
      "h5ad_X_encoding":enc,"h5ad_shape":[n_cells,n_genes],
      "processed_rows":n_run,"selected_cells":selected_cells,
      "candidate_indices":hits,"candidate_genes":CANDIDATES,
      "max_abs_X_rowsum_vs_obs_total_counts":max_rowsum_diff,
      "rows_with_rowsum_difference_gt_1e-4":mismatch_rows,
      "noninteger_sparse_values":noninteger_values,
      "cell_sets":["exclude_predicted_doublets","author_annotated"],
      "created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "environment":{"python":platform.python_version(),"numpy":np.__version__,
                     "pandas":pd.__version__,"h5py":h5py.__version__}
    }
    (out/f"{prefix}extraction_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary),flush=True)
if __name__=="__main__": main()
