"""Cross-check the Stage005 submission artifacts against saved result tables."""
from pathlib import Path
import hashlib
import json
import re
import sys

import pandas as pd
from docx import Document

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'submission_bmc_genomics'
sys.path.insert(0, str(ROOT / 'manuscript'))
from manuscript_content import TITLE, ABSTRACT, REFERENCES, FIGURE_LEGENDS

def text(path):
    d = Document(path)
    return '\n'.join([p.text for p in d.paragraphs] +
                     [c.text for t in d.tables for row in t.rows for c in row.cells])

main = text(OUT/'01_Main_Manuscript_BMC_Genomics.docx')
supp = text(OUT/'supplementary_files/Supplementary_Materials.docx')
sheet = text(OUT/'04_Submission_Copy_Paste_Sheet.docx')
title = text(OUT/'02_Title_Page.docx')
cover = text(OUT/'03_Cover_Letter_BMC_Genomics.docx')
checks = {}
checks['titles_synchronized'] = all(TITLE in s for s in (main,supp,sheet,title,cover))
checks['abstract_synchronized'] = all(v in main and v in sheet for v in ABSTRACT.values())
body = main.split('References\n')[0]
cited = set()
for block in re.findall(r'\[([\d,\s-]+)\]',body):
    for item in block.split(','):
        item=item.strip()
        if '-' in item:
            a,b=map(int,item.split('-')); cited.update(range(a,b+1))
        else:
            cited.add(int(item))
checks['all_21_references_cited'] = cited == set(range(1,len(REFERENCES)+1)) and len(REFERENCES)==21
checks['all_six_figures_cited'] = all(f'Figure {i}' in body for i in range(1,7))
checks['six_figure_legends'] = len(FIGURE_LEGENDS)==6
checks['figure_inventory_synchronized'] = ('6 main figures' in title and
    json.loads((OUT/'PACKAGE_MANIFEST.json').read_text())['main_figures']==6 and
    'Figure5_healthy_reference.png' in sheet and 'Figure6_CT_functional_programs.png' in sheet)
checks['post_result_role_explicit'] = 'exploratory extension after the original cohort results' in main
checks['primary_negative_retained'] = 'association was not confirmed in GSE23597' in main
checks['chemistry_wording_corrected'] = 'not identifiable' not in main+supp
checks['no_internal_response_wording'] = not re.search(r'reviewer-requested|in response to the reviewer',main+supp,re.I)
health=pd.read_csv(ROOT/'runs/005_health_function/healthy_reference_contrasts.tsv',sep='\t')
health=health[health.scope.eq('main')]
checks['all_16_health_contrasts_in_supplement'] = len(health)==16 and all(
    f'{r.estimate:.3f} ({r.lower:.3f}, {r.upper:.3f})' in supp for r in health.itertuples())
ct=pd.read_csv(ROOT/'runs/005_health_function/ct_program_models.tsv',sep='\t')
ct=ct[ct.model.eq('baseline_ANCOVA') & ~ct.program.eq('CANDIDATE6')]
checks['all_five_reference_models_in_supplement'] = len(ct)==5 and all(
    f'{r.estimate:.3f} ({r.lower:.3f}, {r.upper:.3f})' in supp for r in ct.itertuples())
checks['frozen_005_plan_unchanged'] = hashlib.sha256((ROOT/'planning/analysis_plan_005_health_function.json').read_bytes()).hexdigest() == '4b4213df3666d6d98723f2106dc0dfe7ee39bab7b79d3d0f2c54fe91203833b2'
# TSV content is portable across Windows CRLF and Linux LF; the frozen plan
# above deliberately keeps a stricter byte-level hash check.
checks['new_tables_copied_text_exact'] = all(
    p.read_text(encoding='utf-8')==(OUT/'supplementary_files'/('Stage005_'+p.name)).read_text(encoding='utf-8')
    for p in (ROOT/'runs/005_health_function').glob('*.tsv'))
result={'passed':all(checks.values()),'n_checks':len(checks),'checks':checks,
        'failed':[k for k,v in checks.items() if not v]}
(OUT/'PACKAGE_STAGE005_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
if not result['passed']:
    raise SystemExit(1)
