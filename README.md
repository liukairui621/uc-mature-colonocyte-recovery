# UC mature colonocyte recovery

Reproducible analysis for **Recovery and incomplete normalization of a colonocyte expression signature in treated ulcerative colitis**.

## Study question

The project asks whether a six-gene mature absorptive colonocyte score observed in longitudinal UC biopsies reflects a change in captured colonocyte abundance, a change in expression within a defined colonocyte state, or both. The score contains **AQP8, HMGCS2, GUCA2A, CA2, SLC26A3, and MS4A12**.

## Evidence architecture

| Dataset | Role | Main result |
|---|---|---|
| GSE92415 | Exploratory discovery | Clinical-response association in the conditional discovery model |
| GSE23597 | Locked primary validation | Did not meet the prespecified support criterion |
| GSE73661 IFX | Supportive endoscopic analysis | Positive association with endoscopic healing |
| GSE73661 VDZ | Cross-drug extension | Positive, imprecise estimates at weeks 6 and 12 |
| GSE282122 | Patient-level single-cell decomposition | Remission association within CT colonocytes without a matching increase in CT epithelial fraction |

Results are reported separately for each cohort; cross-platform effects and P values are not pooled.

## Healthy-reference and program extension

The post-result Stage005 extension used the original candidate and paired patient sets. Early infliximab and vedolizumab healers retained lower Candidate6 scores than within-study non-IBD controls. Three programs excluding all six candidate genes showed no positive baseline-adjusted remission association or clear change correlation with Candidate6. The complete results, including negative comparisons, are in `runs/005_health_function/` and `reports/Stage005_health_function_20260920.md`.

## Data

All inputs are public and de-identified.

- [GSE92415](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92415)
- [GSE23597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23597)
- [GSE73661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73661)
- [GSE282122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282122)
- [Zenodo record 14007626](https://zenodo.org/records/14007626)

Raw and large processed matrices are not duplicated in this repository. Download manifests and provenance records identify the source files and checksums used.

## Repository map

- `code/`: analysis scripts for metadata, bulk cohorts, single-cell pseudobulk, DecontX and evidence ledgers
- `config/`: fixed analysis configurations
- `planning/`: timestamped internal specifications and amendments
- `provenance/`: source and execution records
- `results/publication_tables/`: compact, machine-readable result tables used in the manuscript
- `figures/publication/`: manuscript and supplementary figures
- `reports/`: complete analysis reports, including negative and sensitivity results
- `manuscript/`: current manuscript text synchronized with the statistical outputs
- `renv.lock`: R dependency lock
- `signature_lock.json`: frozen primary validation specification

## Reproduction

1. Clone the repository and restore the R environment with `renv::restore()`.
2. Download public source files listed in the provenance manifests into a local `inputs/` directory.
3. Run numbered scripts in `code/` in stage order. Each script writes to a stage-specific directory under `runs/`.
4. Run `code/plot_publication_main_figures.py`, `code/figure_003_cell_context.R`, and the stage-specific plotting scripts to regenerate the final figures.
5. Compare generated compact tables with `results/publication_tables/`.

The executed snapshot retains the original absolute Linux analysis root so that hashes and provenance remain interpretable. To work elsewhere without editing the audit snapshot, run `python tools/make_portable_copy.py /path/to/new/working-copy` and use the new copy. Rewritten paths intentionally invalidate frozen file hashes in that portable copy.

Large `inputs/` and `runs/` directories are intentionally ignored. Some single-cell steps require substantial memory and disk space; the exact processed input is recorded in the stage-003 provenance.

## Reproducibility boundaries

The candidate was assembled after inspection of the discovery dataset. The GSE23597 primary model was fixed internally before validation expression access; this was not a public preregistration. The VDZ matrix and sample-level scores had already been accessed before its outcome-analysis rules were fixed. These timing distinctions are preserved in the repository.

## License

Code is released under the MIT License. Dataset rights remain with the original data producers and repositories.
