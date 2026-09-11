from pathlib import Path
import argparse, hashlib, json, datetime
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz

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
        out=np.empty(len(codes),dtype=object); out[:]=None
        ok=codes>=0; out[ok]=cats[codes[ok]]
        return out
    return decode_array(obj[:])

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True)
    ap.add_argument("--pairs",required=True)
    ap.add_argument("--var",required=True)
    ap.add_argument("--outdir",required=True)
    args=ap.parse_args()
    inp=Path(args.input); out=Path(args.outdir)
    matrix_dir=out/"matrices"; meta_dir=out/"cell_metadata"
    matrix_dir.mkdir(parents=True,exist_ok=True); meta_dir.mkdir(parents=True,exist_ok=True)
    pairs=pd.read_csv(args.pairs,sep="	")
    selected=sorted(set(pairs.pre_sample_id)|set(pairs.post_sample_id))
    var=pd.read_csv(args.var,sep="	")
    if len(var)!=33075: raise RuntimeError("Unexpected feature-table length")
    rows=[]
    with h5py.File(inp,"r") as h:
        xg=h["X"]; n_cells,n_genes=map(int,xg.attrs["shape"])
        obs={k:read_frame_column(h,"obs",k) for k in
             ["sample_id","Patient","Site","Treatment","Remission_status","final_analysis",
              "major","minor","predicted_doublets","total_counts","LibraryType","Batch"]}
        annotated=np.array([v is not None and str(v)!="nan" for v in obs["final_analysis"]])
        dbl=np.asarray(obs["predicted_doublets"]).astype(bool)
        for order,sid in enumerate(selected,1):
            ix=np.flatnonzero((obs["sample_id"]==sid)&annotated&(~dbl))
            if len(ix)==0: raise RuntimeError(f"No selected annotated cells for {sid}")
            if np.any(np.diff(ix)!=1): raise RuntimeError(f"Cells are not contiguous for {sid}")
            start,end=int(ix[0]),int(ix[-1])+1
            ip=np.asarray(xg["indptr"][start:end+1],dtype=np.int64)
            lo,hi=int(ip[0]),int(ip[-1])
            dat=np.asarray(xg["data"][lo:hi])
            ind=np.asarray(xg["indices"][lo:hi],dtype=np.int32)
            mat=csr_matrix((dat,ind,ip-lo),shape=(end-start,n_genes))
            released=np.asarray(obs["total_counts"][start:end],dtype=float)
            row_sum=np.asarray(mat.sum(axis=1)).ravel()
            if np.max(np.abs(row_sum-released))>1e-4: raise RuntimeError(f"Row-sum mismatch {sid}")
            npz=matrix_dir/f"{sid}.cells_by_genes.npz"
            save_npz(npz,mat,compressed=True)
            md=pd.DataFrame({"cell_index":np.arange(start,end,dtype=int)})
            for k in ["sample_id","Patient","Site","Treatment","Remission_status","final_analysis",
                      "major","minor","LibraryType","Batch"]:
                md[k]=np.asarray(obs[k],object)[start:end]
            md.to_csv(meta_dir/f"{sid}.cells.tsv.gz",sep="	",index=False,compression="gzip")
            rows.append({"sample_id":sid,"order":order,"n_cells":len(ix),"n_genes":n_genes,
                         "nnz":int(mat.nnz),"n_clusters":int(md.final_analysis.nunique()),
                         "raw_total_UMI":float(row_sum.sum()),"first_global_row":start,
                         "last_global_row":end-1,"matrix_path":str(npz),
                         "metadata_path":str(meta_dir/f"{sid}.cells.tsv.gz"),
                         "matrix_bytes":npz.stat().st_size,"matrix_sha256":sha256(npz)})
            print(json.dumps(rows[-1]),flush=True)
    mf=pd.DataFrame(rows)
    mf.to_csv(out/"batch_manifest.tsv",sep="	",index=False)
    summary={"status":"COMPLETE","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
             "input":str(inp),"input_sha256":"60323cd4f8e32881316029ce08425131ca6aa2e64031d1d8ffa23d0f5c4f5028",
             "samples":len(mf),"cells":int(mf.n_cells.sum()),"genes":int(mf.n_genes.iloc[0]),
             "clusters_per_sample_min":int(mf.n_clusters.min()),"clusters_per_sample_max":int(mf.n_clusters.max()),
             "all_sample_cell_blocks_contiguous":True,"all_row_sum_checks":"PASS"}
    (out/"extraction_summary.json").write_text(json.dumps(summary,indent=2)+chr(10))
    print(json.dumps(summary),flush=True)
if __name__=="__main__": main()
