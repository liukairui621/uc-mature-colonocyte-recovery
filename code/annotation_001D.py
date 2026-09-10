from pathlib import Path
import csv,gzip,re,json,collections,datetime,hashlib
P=Path('/root/projects/UC_Treatment_Recovery'); O=P/'runs/001D_annotation'
H=list(csv.DictReader((P/'inputs/annotation/hgnc_complete_set_20260910.txt').open(),delimiter='\t'))
approved={x['symbol']:x['hgnc_id'] for x in H if x['status']=='Approved'}
aliases=collections.defaultdict(set)
for x in H:
 if x['status']!='Approved':continue
 for key in ['prev_symbol','alias_symbol']:
  for s in x[key].split('|'):
   if s:aliases[s].add(x['symbol'])
def canon(s):
 s=s.strip()
 if s in approved:return s,'approved',s
 tar=sorted(aliases.get(s,set()))
 if len(tar)==1:return tar[0],'unique_alias_or_previous',';'.join(tar)
 return '', 'ambiguous_alias' if tar else 'unmapped',';'.join(tar)
def write(name,rows,fields=None):
 if fields is None:fields=list(rows[0])
 with (O/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
policy={'phase':'pre-validation annotation amendment','validation_expression_accessed':False,'dictionary':'HGNC complete set retrieved 2026-09-10; source Last-Modified 2026-09-09','policy':'Approved exact symbol takes precedence. Otherwise use union of previous and alias symbols only if exactly one approved target. Ambiguous or unknown symbols excluded and audited. Original multi-symbol probe annotations remain excluded even if tokens could converge. Mean log2 probes AFTER canonicalization. No expression/outcome-based mapping. Canonicalize program memberships using same dictionary; retain original membership and rejected entries.','candidate_membership':'unchanged six genes','platform_common_policy':'Primary inflammation covariate will use the canonical Hallmark inflammatory-response genes retained as unambiguous probes on ALL three platforms, fixed before validation expression. All members of this frozen common panel required in each actual expression matrix; no cohort-specific averaging. Original Hallmark200 coverage also reported with original >=80% rule.'}
(O/'annotation_policy.json').write_text(json.dumps(policy,indent=2))
universes={}; allrows=[]; platforms={}
for gpl in ['GPL13158','GPL570','GPL6244']:
 with gzip.open(P/('inputs/annotation/GPL6244.relay.annot.gz' if gpl=='GPL6244' else f'inputs/annotation/{gpl}.annot.gz'),'rt') as f:
  for line in f:
   if line.startswith('ID\t'): header=line.rstrip('\n').split('\t');break
  rows=csv.DictReader(f,fieldnames=header,delimiter='\t',quoting=csv.QUOTE_NONE)
  data=[]
  for r in rows:
   probe=r['ID']
   if not probe or probe.startswith('!'):continue
   old=(r.get('Gene symbol') or '').strip()
   pre=bool(old) and not re.search(r'///|//|;|\|',old) and not probe.startswith('AFFX')
   if pre:c,reason,tar=canon(old)
   else:c,reason,tar='','original_ambiguous_blank_or_control',''
   data.append(dict(probe=probe,original_symbol=old,canonical_symbol=c,mapping_status=reason,possible_targets=tar,retained=int(bool(c))))
  assert len({r['probe'] for r in data})==len(data)
  write(f'{gpl}_canonical_probe_mapping.tsv',data)
  universes[gpl]={r['canonical_symbol'] for r in data if r['retained']}
  platforms[gpl]={'probes':len(data),'retained_probes':sum(r['retained'] for r in data),'canonical_genes':len(universes[gpl]),'mapping_status_counts':dict(collections.Counter(r['mapping_status'] for r in data))}
  allrows+=data
# Discovery matrix probe subset limits actual platform presence.
actual_probes={r['probe'] for r in csv.DictReader((P/'runs/001A_qc/probe_gene_mapping.tsv').open(),delimiter='\t')}
d=list(csv.DictReader((O/'GPL13158_canonical_probe_mapping.tsv').open(),delimiter='\t'))
universes['GSE92415_actual']={r['canonical_symbol'] for r in d if r['probe'] in actual_probes and r['retained']=='1'}
# Preserve source and canonical identities; do not fill unresolved members.
mem=list(csv.DictReader((P/'inputs/annotation/program_membership_001C.tsv').open(),delimiter='\t'))
audit=[]; newmem=[]; seen=set()
for r in mem:
 c,reason,tar=canon(r['gene']);a=dict(r);a.update(canonical_symbol=c,mapping_status=reason,possible_targets=tar,retired_comparator=int(r['program'].startswith('IBrD')))
 audit.append(a)
 if r['program'].startswith('IBrD') or not c:continue
 key=(r['program'],c)
 if key in seen:continue
 seen.add(key);z=dict(r);z['gene']=c;newmem.append(z)
write('membership_mapping_audit.tsv',audit)
write('program_membership_canonical.tsv',newmem)
sets=collections.defaultdict(set)
orig=collections.defaultdict(set)
for r in mem:
 if not r['program'].startswith('IBrD'):orig[r['program']].add(r['gene'])
for r in newmem:sets[r['program']].add(r['gene'])
cov=[]
for plat,U in universes.items():
 for prog,S in sets.items():
  cov.append({'platform_or_matrix':plat,'program':prog,'source_n':len(orig[prog]),'canonical_n':len(S),'measured_n':len(U&S),'fraction_source':len(U&S)/len(orig[prog]),'missing_canonical':';'.join(sorted(S-U)),'unresolved_source_symbols':';'.join(sorted(x for x in orig[prog] if not canon(x)[0]))})
write('platform_program_coverage.tsv',cov)
infl=sets['HALLMARK_INFLAMMATORY_RESPONSE']
common=infl.intersection(*[universes[k] for k in ['GPL13158','GPL570','GPL6244']])
assert len(common)>=160
write('frozen_common_inflammation_members.tsv',[{'gene':x,'weight':1} for x in sorted(common)])
assert common<=universes['GSE92415_actual']
summary={'status':'COMPLETE','validation_expression_accessed':False,'platforms':platforms,'common_inflammation_n':len(common),'common_missing_from_original_canonical':sorted(infl-common),'all_reference_sets':len(sets),'IBrD':'Retired from substantive comparison due feature incompatibility and original model unavailable; historical outputs retained only as provenance. No assumption scRNA will reproduce it.','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'annotation_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
for row in cov:
 if row['program'] in ['MATURE_COLONOCYTE6_EXPLORATORY','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_INTERFERON_GAMMA_RESPONSE']:print(row)
