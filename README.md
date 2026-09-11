# Current status: Stage003 complete and independently reviewed

Current report: reports/Stage003_cell_context_results_20260910.md.
State ledger: provenance/manifests/stage_003_completion.json.
GSE23597 primary clinical-response validation did not meet its frozen support rule. GSE73661 supplies supportive endoscopic association with construct-overlap limits. GSE282122 supports a within-state epithelial component in a distinct adalimumab clinical-remission context and cannot replace the failed primary validation.

# UC treatment-associated recovery study

Stable project root: /root/projects/UC_Treatment_Recovery

The agreed design is archived in planning/Study_Design_v1.md. Execution results are stored on this server. Discovery uses GSE92415 only. GSE23597 and GSE73661 currently contribute metadata only; no held-out expression has been used for candidate selection.

## Completed first-stage execution

1. python3 code/metadata_001A.py
2. OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 Rscript code/qc_001A.R
3. OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 Rscript code/discovery_001B.R
4. OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 Rscript code/pathways_001B.R
5. python3 code/clarify_go_membership.py
6. Rscript code/plot_program_overview_001B.R
7. python3 code/build_stage1_report.py <review-status-text>

Use the input paths and hashes in provenance/manifests. The complete discovery matrix is GSE92415_series_matrix.relay.txt.gz; the other partial transfer is not an analysis input. Source normalization is RMA and is retained. The bootstrap download script is a retrieval utility, not the analysis pipeline.

Do not blindly overwrite this run for a new scientific decision. Create a new run/branch with the same frozen inputs and declared changes. These first-stage results are exploratory; no publication signature has been frozen.

## Next bounded stage

- Derive patient-level scores and assess relationships with changes in inflammatory burden and clinical response in the discovery set.
- Refine broad epithelial annotations with UC-specific cell-state information; distinguish epithelial cell abundance, differentiated functions, and stress/regeneration.
- Compare concrete candidates with existing inflammation, MARS and IBrD concepts before novelty claims.
- Keep baseline prediction, post-treatment association, and endoscopic outcomes separate.
- Freeze candidates and scoring rules before evaluating independent treatment cohorts.

The three immediate follow-up needs are the public subject-count inconsistency, source-level exposure reconstruction for validation metadata, and cell-state-specific rather than generic pathway interpretation.

## Reproducibility and status

provenance/environments contains runtime/package versions; provenance/manifests contains hashes and commands; provenance/reviews contains independent checks. logs retains original transfer/performance restart records. tests/model_smoke_001B.json records checks before full paired statistics.

reports/Stage1_report_20260910.md is the stage report, with direct references to source exports. Core arrays, gene mapping, all contrast tables, and all pathway results are retained, including negative results and sensitivity branches.
