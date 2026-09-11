# Current status: Stage003E ambient diagnostics complete

Current report: reports/Stage003E_ambient_negative_control_review_20260910.md.
State ledger: provenance/manifests/stage_003E_completion.json.
GSE23597 primary clinical-response validation did not meet its frozen support rule. GSE73661 supplies a supportive endoscopic association with construct-overlap limits. In GSE282122, the common-patient decomposition and raw-count state analysis support a composition-beyond CT within-state component in a distinct adalimumab clinical-remission context. A post-result DecontX sensitivity using released state labels does not overturn the CT association (BH over the prespecified 15-state family, 10 evaluable = 0.0169), while TA and LGR5-positive stem do not retain fixed-family support (both BH = 0.1399). Because this model treats cluster-typical and cluster-atypical signals asymmetrically, the result is method-qualified. Fixed PTPRC and COL1A1 controls are not concordant with Candidate6 but do not exclude gene-specific ambient or background-composition effects. The broad cross-lineage claim is withdrawn, and no GSE282122 result can replace the failed primary validation.

# UC treatment-associated recovery study

Stable project root: /root/projects/UC_Treatment_Recovery

The agreed design is archived in planning/Study_Design_v1.md. Execution results are stored on this server. Discovery uses GSE92415 only. Held-out GSE23597 and GSE73661 expression and GSE282122 single-cell counts have now been analyzed only in their frozen validation, supportive, and cell-context roles; none was used to retune the candidate.

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
