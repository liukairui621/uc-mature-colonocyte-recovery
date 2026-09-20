"""Extract fixed reference-gene counts in the original released CT compartment."""
from pathlib import Path
import argparse, collections, datetime, json, platform
import h5py, numpy as np, pandas as pd
from scipy.sparse import csr_matrix
from extract_003_h5ad import read_frame_column, hash_file

ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--smoke',action='store_true'); args=ap.parse_args()
    out=ROOT/('runs/005_ct_smoke' if args.smoke else 'runs/005_ct_extraction')
    out.mkdir(parents=True,exist_ok=True)
    if (out/'COMPLETE.json').exists(): raise RuntimeError('Completed extraction exists')
    plan=json.loads((ROOT/'planning/analysis_plan_005_health_function.json').read_text())
    desired=sorted(set(sum(plan['programs'].values(),[])))
    pairs=pd.read_csv(ROOT/'runs/003_preexpression/fixed_UC_colonic_site_pairs.tsv',sep='\t')
    samples=sorted(set(pairs.pre_sample_id)|set(pairs.post_sample_id))
    var=pd.read_csv(ROOT/'inputs/002B_asset_audit/TAURUS_actual_var_features_full.tsv',sep='\t')
    hgnc=pd.read_csv(ROOT/'inputs/annotation/hgnc_complete_set_20260910.txt',sep='\t',dtype=str).fillna('')
    hgnc=hgnc[hgnc.status.eq('Approved')]; approved=set(hgnc.symbol); aliases=collections.defaultdict(set)
    for r in hgnc.to_dict('records'):
        for col in ['prev_symbol','alias_symbol']:
            for old in r[col].split('|'):
                if old: aliases[old].add(r['symbol'])
    def canon(g):
        return g if g in approved else (next(iter(aliases[g])) if len(aliases[g])==1 else None)
    mapped=[canon(str(g)) for g in var.gene_symbol]
    indices=np.array([i for i,g in enumerate(mapped) if g in desired],int)
    genes=sorted(set(mapped[i] for i in indices))
    to_gene=np.array([genes.index(mapped[i]) for i in indices],int)
    projection=csr_matrix((np.ones(len(indices)),(np.arange(len(indices)),to_gene)),shape=(len(indices),len(genes)))
    pd.DataFrame([{'feature_row':int(i),'original_symbol':var.gene_symbol.iloc[i],'canonical_symbol':mapped[i]} for i in indices]).to_csv(out/'feature_mapping.tsv',sep='\t',index=False)
    coverage=pd.DataFrame([{'program':s,'n_members':len(gs),'n_measured':len(set(gs)&set(genes)),'coverage':len(set(gs)&set(genes))/len(gs),'missing':';'.join(sorted(set(gs)-set(genes)))} for s,gs in plan['programs'].items()])
    coverage.to_csv(out/'coverage.tsv',sep='\t',index=False)
    counts=np.zeros((len(samples),len(genes))); detected=np.zeros_like(counts)
    n_cells=np.zeros(len(samples),dtype=int); libraries=np.zeros(len(samples)); sample_index={s:i for i,s in enumerate(samples)}
    source=ROOT/'inputs/singlecell/GSE282122_Zenodo14007626/TAURUS_raw_counts_annotated_final.h5ad'
    max_rowsum_error=0.; n_noninteger=0
    with h5py.File(source,'r') as h:
        x=h['X']; shape=tuple(x.attrs['shape']); assert shape[1]==len(var)
        obs={key:read_frame_column(h,'obs',key) for key in ['sample_id','Disease','final_analysis','predicted_doublets','total_counts']}
        selected=(np.isin(obs['sample_id'],samples) & (obs['Disease']=='UC') & (obs['final_analysis']=='Non ileal CT colonocyte') & ~obs['predicted_doublets'].astype(bool))
        if args.smoke:
            selected_ids=np.flatnonzero(selected)
            assert len(selected_ids)>0
            stop=min(shape[0],int(selected_ids[0])+25000)
        else: stop=shape[0]
        for start in range(0,stop,25000):
            end=min(start+25000,stop); local=np.flatnonzero(selected[start:end])
            if not len(local): continue
            ip=np.array(x['indptr'][start:end+1],dtype=np.int64); lo,hi=int(ip[0]),int(ip[-1])
            vals=np.array(x['data'][lo:hi]); ind=np.array(x['indices'][lo:hi])
            mat=csr_matrix((vals,ind,ip-lo),shape=(end-start,shape[1]))[local,:]
            total=np.asarray(mat.sum(axis=1)).ravel()
            max_rowsum_error=max(max_rowsum_error,float(np.max(np.abs(total-np.asarray(obs['total_counts'][start:end],float)[local]))))
            n_noninteger+=int(np.count_nonzero(np.abs(mat.data-np.rint(mat.data))>1e-7))
            sub=(mat[:,indices]@projection).toarray()
            ids=np.array([sample_index[str(s)] for s in obs['sample_id'][start:end][local]])
            np.add.at(counts,ids,sub); np.add.at(detected,ids,sub>0)
            np.add.at(n_cells,ids,1); np.add.at(libraries,ids,total)
            print(json.dumps({'through_cell':end,'CT_cells':int(n_cells.sum())}),flush=True)
    assert max_rowsum_error<1e-4 and n_noninteger==0
    frame=pd.DataFrame(counts,columns=genes); frame.insert(0,'total_UMI',libraries); frame.insert(0,'n_cells',n_cells); frame.insert(0,'sample_id',samples)
    frame.to_csv(out/'sample_gene_counts.tsv.gz',sep='\t',index=False)
    det=pd.DataFrame(detected,columns=genes);det.insert(0,'sample_id',samples);det.to_csv(out/'sample_gene_detected_cells.tsv.gz',sep='\t',index=False)
    summary={'status':'TESTED' if args.smoke else 'COMPLETE','CT_cells':int(n_cells.sum()),'samples_with_CT':int((n_cells>0).sum()),'features':len(genes),'max_rowsum_error':max_rowsum_error,'noninteger_counts':n_noninteger,'source_path':str(source),'source_bytes':source.stat().st_size,'source_sha256':None if args.smoke else hash_file(source),'plan_sha256':hash_file(ROOT/'planning/analysis_plan_005_health_function.json'),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'h5py':h5py.__version__,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (out/'COMPLETE.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__': main()
