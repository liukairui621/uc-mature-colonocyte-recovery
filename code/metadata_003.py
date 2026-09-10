from pathlib import Path
import pandas as pd, json, hashlib, datetime, platform
p=Path(__file__).resolve().parents[1]
out=p/'runs/003_preexpression'
out.mkdir(parents=True,exist_ok=True)
m=pd.read_csv(p/'inputs/002B_asset_audit/GSE282122_verified_sample_metadata.tsv',sep='\t')
a=pd.read_csv(p/'inputs/002B_asset_audit/paired_sample_list.csv',sep='\t')
author=set(a['sample_id'].astype(str))
u=m[(m['disease']=='UC') & (m['site']!='Terminal_Ileum') & m['treatment'].isin(['Pre','Post'])].copy()
pairs=[]
for (patient,site),g in u.groupby(['patient','site'],sort=True):
    pre=g[g['treatment']=='Pre']; post=g[g['treatment']=='Post']
    if len(pre)==1 and len(post)==1:
        x,y=pre.iloc[0],post.iloc[0]
        status=x['Remission']
        if status not in ['Remission','Non_Remission']:
            continue
        pairs.append({
          'patient':patient,'site':site,'pre_sample_id':x['sample_id_normalized'],
          'post_sample_id':y['sample_id_normalized'],'remission':status,
          'pre_GSM':x['GSM'],'post_GSM':y['GSM'],
          'library_type':x['librarytype'],'batch':x['batch'],
          'pre_inflammation_label':x['inflammation'],
          'post_inflammation_label':y['inflammation'],
          'pre_inflammation_score':x['inflammation score'],
          'post_inflammation_score':y['inflammation score'],
          'both_in_author_paired_list':x['sample_id_normalized'] in author and y['sample_id_normalized'] in author
        })
d=pd.DataFrame(pairs)
d.to_csv(out/'fixed_UC_colonic_site_pairs.tsv',sep='\t',index=False)
pt=(d.groupby(['patient','remission','library_type','batch'],as_index=False)
      .agg(n_site_pairs=('site','size'),n_author_site_pairs=('both_in_author_paired_list','sum')))
pt.to_csv(out/'fixed_patient_structure.tsv',sep='\t',index=False)
cross=(pt.groupby(['remission','library_type','batch'],as_index=False)
         .agg(patients=('patient','nunique'),site_pairs=('n_site_pairs','sum')))
cross.to_csv(out/'outcome_chemistry_batch_crosswalk.tsv',sep='\t',index=False)
summary={
 'status':'COMPLETE','expression_X_accessed':False,'composition_counts_accessed':False,
 'UC_colonic_site_pairs':int(len(d)),'UC_patients':int(d.patient.nunique()),
 'remission_patients':int(pt.remission.eq('Remission').sum()),
 'nonremission_patients':int(pt.remission.eq('Non_Remission').sum()),
 'remission_site_pairs':int(d.remission.eq('Remission').sum()),
 'nonremission_site_pairs':int(d.remission.eq('Non_Remission').sum()),
 'author_both_site_pairs':int(d.both_in_author_paired_list.sum()),
 'author_both_patients':int(d[d.both_in_author_paired_list].patient.nunique()),
 'all_remitters_v31':bool((pt.loc[pt.remission.eq('Remission'),'library_type']=='3 prime v3.1').all()),
 'metadata_boundaries':['Patient is inferential unit; site pairs are repeated measures collapsed equally within patient',
                        'Chemistry and batch are constant within each patient site pair but imbalanced by outcome',
                        'Outcome is clinical remission after adalimumab, not GSE73661 endoscopic healing',
                        'Cell-state counts and expression have not been accessed in this preparation'],
 'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'python':platform.python_version()
}
(out/'metadata_design_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
