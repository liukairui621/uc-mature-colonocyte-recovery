from pathlib import Path
import csv,json,datetime
root=Path('/root/projects/UC_Treatment_Recovery')
for branch in ['primary_paired','barcode_sensitivity']:
    folder=root/'runs/001B_discovery'/branch
    src=folder/'GO_BP_coverage.tsv'
    rows=list(csv.DictReader(src.open(),delimiter='\t'))
    updated=[]
    for row in rows:
        updated.append({'set':row['set'],'source_members_in_tested_universe':row['original_size'],'members_used_in_test':row['mapped_size'],'retained_fraction_within_tested_universe':row['coverage'],'scope':'GO membership was built from assayed mapped genes; this is not full-database pathway coverage'})
    with (folder/'GO_BP_membership_scope_v2.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(updated[0]),delimiter='\t');w.writeheader();w.writerows(updated)
note={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'change':'Clarify GO membership denominator; no change to genes tested, statistics, thresholds or models','deprecated_presentation':'*/GO_BP_coverage.tsv original_size and coverage are within tested universe, not global GO coverage','replacement':'*/GO_BP_membership_scope_v2.tsv','Hallmark':'Full original gene-set coverage remains valid'}
(root/'provenance/reviews/go_membership_clarification.json').write_text(json.dumps(note,indent=2))
print(json.dumps(note))
