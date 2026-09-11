# UC mature colonocyte recovery

Reproducible analysis for **Longitudinal bulk and single-cell transcriptomic reanalysis characterizes mucosal healing-associated recovery of a mature absorptive colonocyte score in ulcerative colitis**.

## Study question

The project asks whether a six-gene mature absorptive colonocyte score observed in longitudinal UC biopsies reflects a change in captured colonocyte abundance, a change in expression within a defined colonocyte state, or both. The score contains **AQP8, HMGCS2, GUCA2A, CA2, SLC26A3, and MS4A12**.

## Evidence architecture

| Dataset | Role | Main result |
|---|---|---|
| GSE92415 | Exploratory discovery | Clinical-response association in the conditional discovery model |
| GSE23597 | Locked primary validation | Did not meet the internally fixed support criterion |
| GSE73661 IFX | Supportive endoscopic analysis | Positive association with endoscopic healing |
| GSE73661 VDZ | Cross-drug extension | Positive, imprecise estimates at weeks 6 and 12 |
| GSE282122 | Patient-level single-cell decomposition | Larger CT-state score rebound from a lower baseline without a matching increase in CT epithelial fraction; attenuated after baseline adjustment |

The positive supportive analyses do not replace the negative primary validation. Cross-platform effects and P values are not pooled.

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
4. Compare generated compact tables with `results/publication_tables/`.

The executed snapshot retains the original absolute Linux analysis root so that hashes and provenance remain interpretable. To work elsewhere without editing the audit snapshot, run `python tools/make_portable_copy.py /path/to/new/working-copy` and use the new copy. Rewritten paths intentionally invalidate frozen file hashes in that portable copy.

Large `inputs/` and `runs/` directories are intentionally ignored. Some single-cell steps require substantial memory and disk space; the exact processed input is recorded in the stage-003 provenance.

## Reproducibility boundaries

The candidate was assembled after inspection of the discovery dataset. The GSE23597 primary model was fixed internally before validation expression access; this was not a public preregistration. The VDZ matrix and sample-level scores had already been accessed before its outcome-analysis rules were fixed. These timing distinctions are preserved in the repository.

## License

Code is released under the MIT License. Dataset rights remain with the original data producers and repositories.
