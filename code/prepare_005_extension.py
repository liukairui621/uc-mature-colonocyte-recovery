"""Resolve outcome-independent gene sets and lock the post-result extension."""
from pathlib import Path
import collections, datetime, hashlib, json, platform, urllib.request
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/005_preparation'
OUT.mkdir(parents=True, exist_ok=True)
INPUT = ROOT / 'inputs/005_reference_sets'
INPUT.mkdir(parents=True, exist_ok=True)
CAND = ['AQP8','HMGCS2','GUCA2A','CA2','SLC26A3','MS4A12']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    h = pd.read_csv(ROOT/'inputs/annotation/hgnc_complete_set_20260910.txt',sep='\t',dtype=str).fillna('')
    h = h[h.status.eq('Approved')]
    approved = set(h.symbol); aliases = collections.defaultdict(set)
    for r in h.to_dict('records'):
        for col in ['prev_symbol','alias_symbol']:
            for old in r[col].split('|'):
                if old: aliases[old].add(r['symbol'])
    def canon(g):
        if g in approved: return g
        return next(iter(aliases[g])) if len(aliases[g]) == 1 else None
    source = pd.read_csv(ROOT/'runs/001D_annotation/program_membership_canonical.tsv',sep='\t')
    ct_name = 'UC_Maciag2024_TableS2_CT_Colonocytes_reported_marker_mean_proxy'
    sets = {'CT_MARKERS_EXCLUDING_CANDIDATE6': sorted(set(source.loc[source.program.eq(ct_name),'gene'])-set(CAND))}
    manifest=[]
    for name in ['GOBP_INTESTINAL_ABSORPTION','GOBP_BRUSH_BORDER_ASSEMBLY']:
        url = 'https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName='+name+'&fileType=json'
        dest=INPUT/(name+'.json')
        if not dest.exists():
            req=urllib.request.Request(url,headers={'User-Agent':'UCRecoveryResearch/1.0'})
            dest.write_bytes(urllib.request.urlopen(req,timeout=60).read())
        dat=json.loads(dest.read_text())
        obj=dat[name]
        assert obj.get('exactSource','').startswith('GO:'), obj
        genes=obj['geneSymbols']
        sets[name]=sorted({canon(g) for g in genes if canon(g)}-set(CAND))
        manifest.append({'program':name,'url':url,'source_n':len(genes),'source_metadata':obj,'sha256':sha(dest)})
    hall={}
    for line in (ROOT/'inputs/annotation/h.all.v2025.1.Hs.symbols.gmt').read_text().splitlines():
        x=line.split('\t'); hall[x[0]]=x[2:]
    sets['HALLMARK_E2F_TARGETS']=sorted({canon(g) for g in hall['HALLMARK_E2F_TARGETS'] if canon(g)}-set(CAND))
    sets['COMMON194_INFLAMMATION']=sorted(set(pd.read_csv(ROOT/'runs/001D_annotation/frozen_common_inflammation_members.tsv',sep='\t').gene)-set(CAND))
    sets['CANDIDATE6']=CAND
    membership=pd.DataFrame([{'program':s,'gene':g,'role':'primary_function' if s in list(sets)[:3] else 'context'} for s,gs in sets.items() for g in gs])
    membership.to_csv(OUT/'program_membership.tsv',sep='\t',index=False)
    overlap=[]
    for s,gs in sets.items():
        for t,hs in sets.items():
            overlap.append({'program1':s,'program2':t,'n_shared':len(set(gs)&set(hs)),'shared':';'.join(sorted(set(gs)&set(hs)))})
    pd.DataFrame(overlap).to_csv(OUT/'program_overlaps.tsv',sep='\t',index=False)
    (OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    plan={
      'stage':'005','status':'FROZEN','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'label':'Post-result exploratory extension specified after all original cohort results; not a blinded validation or public preregistration',
      'questions':['How far do post-treatment expression levels remain from within-study non-IBD controls?','Do non-overlapping mature-identity, intestinal-absorption and brush-border programs show coordinated within-CT changes?'],
      'unchanged':'Candidate6 membership/weights, original validation result, cohort roles, CT labels and original frozen results remain unchanged.',
      'programs':sets,'membership_sha256':sha(OUT/'program_membership.tsv'),
      'healthy_reference':{
        'datasets':['GSE92415','GSE73661'],
        'contexts':['GSE92415_W6','GSE73661_IFX_W4W6','GSE73661_VDZ_W6','GSE73661_VDZ_W12'],
        'patients':'Original verified paired patient sets; their outcome definitions unchanged. Discovery includes its original placebo and GLM patients; report treatment-stratified descriptives.',
        'controls':'GSE92415 named non-IBD samples; GSE73661 control arrays grouped by released subject ID, repeated Control_7 averaged at gene level. Repeat GSE73661 contrasts after omitting Control_7 as a sensitivity.',
        'normalization':'Reuse canonical deposited log2 matrices without cohort pooling or new batch correction. Candidate6 and fixed common194 gene means.',
        'main_contrasts':'Post-treatment mean minus control mean for Candidate6 and common194 inflammation in each of two outcome groups in four contexts (16 Welch tests, single BH family n=16). Report original log2-score units and 95% CI. Baseline levels and paired changes are descriptive with patient-level 95% CI.',
        'reference_scale':'Within-study control-mean centering; control SD standardization for displays only. No clinical cutoff, equivalence, recovery percentage or normalization claim from nonsignificance.',
        'singlecell_controls':'No healthy single-cell hypothesis tests: only 3 source control donors and sampling/chemistry differences. Bulk healthy comparison and CT function extension have distinct roles.'},
      'ct_function':{
        'primary_programs':list(sets)[:3],
        'context_programs':['COMMON194_INFLAMMATION','HALLMARK_E2F_TARGETS'],
        'cell_sample_rule':'Same corrected h5ad, released non-ileal CT label, exclude predicted doublets, original fixed same-site UC pairs, >=20 CT cells at both visits; equal eligible-site weights within patient.',
        'feature_rule':'Same HGNC 2026-09-10 approved-exact/unique-alias rule; sum count columns resolving to the same canonical symbol. Exclude all six candidate genes from every reference program.',
        'coverage':'Require >=80% of fixed canonical members and >=5 measured genes; absent members not replaced. No filtering or member selection by expression/outcome.',
        'score':'Unweighted mean log2((gene_pseudobulk_count+0.5)/(all_gene_pseudobulk_UMI+1)*1e6).',
        'primary_model':'Follow-up program score ~ remission + baseline program score; HC3 95% CI and residual-t P, BH across three planned programs (n=3 even if non-evaluable).',
        'sensitivities':['Unadjusted patient change contrast','Primary ANCOVA restricted to 10x v3.1','Primary ANCOVA additionally adjusting concurrent source inflammation-score change','Primary ANCOVA with leave-one-gene-out scores to describe component dependence'],
        'coordination':'Pearson correlations of patient program changes with Candidate6 change; Fisher-z 95% CI, BH family n=3. Context scores described separately with BH n=2 for ANCOVA.',
        'gene_display':'Report all member genes and their mean changes, expression/detection summaries and direction concordance; no selective gene-level significance claims.',
        'replication_check':'Reproduce original Candidate6 baseline/follow-up/delta values and n=16, remission n=6 before adopting new results.'},
      'stopping_rule':'Complete these fixed analyses, technical checks and report/manuscript integration. Do not replace a failed program or add a cohort based on resulting P values.',
      'interpretation':'All estimates are transcript-expression associations, not direct transport measurements, mechanistic causation or new clinical prediction. Cross-program coordination may accompany inflammation resolution.',
      'sources':manifest,'environment':{'python':platform.python_version()}
    }
    target=ROOT/'planning/analysis_plan_005_health_function.json'
    if target.exists(): raise RuntimeError('Plan already exists; do not overwrite frozen specification')
    target.write_text(json.dumps(plan,indent=2)+'\n')
    print(json.dumps({'plan':str(target),'sha256':sha(target),'program_sizes':{s:len(g) for s,g in sets.items()}},indent=2))
if __name__=='__main__': main()
