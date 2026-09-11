# Reproducing the UC mature colonocyte analysis

The compact tables and figures in this repository were generated from public GEO and Zenodo inputs. Large matrices are excluded from Git. The manifests under `provenance/` record the source files, checksums and executed analysis root.

## Environments

Restore the R 4.3.3 environment from the project root:

```r
renv::restore()
```

Create a Python 3.12 environment and install the pinned packages:

```bash
python -m pip install -r requirements.txt
```

## Analysis order

| Stage | Main command | Output role |
|---|---|---|
| 001A | `Rscript code/qc_001A.R` | Discovery expression QC and verified pairs |
| 001B–001C | `Rscript code/discovery_001D.R` followed by `Rscript code/programs_001C.R` | Exploratory transcriptome and 33-program comparison |
| 002A | `Rscript code/validation_002A.R` | Locked GSE23597 clinical-response validation |
| 002B | `Rscript code/support_73661_002B.R` | GSE73661 infliximab endoscopic analysis |
| 002D | `Rscript code/support_002D_VDZ.R` | Vedolizumab week-6 and week-12 extensions |
| 003 | `python code/extract_003_h5ad.py` followed by `Rscript code/analyze_003_cell_context.R` | Patient-level single-cell composition/state decomposition |
| 003D | `python code/run_003D_decontx_parallel.py` followed by `Rscript code/analyze_003D_decontx_states.R` | Per-sample DecontX and fixed 15-state family |
| 003E | `python code/analyze_003E_ambient_negative_controls.py` | Ambient-RNA negative controls |
| 003F | `python code/analyze_003F_baseline_adjustment.py` | Baseline, inflammation and chemistry-aware qualification |
| 004A | `python code/build_004A_external_analysis_ledger.py` | Complete external analysis ledger |

Some earlier scripts preserve the absolute Linux path used for the audited execution. To create a writable path-adjusted copy without modifying the frozen snapshot, run:

```bash
python tools/make_portable_copy.py /path/to/new/working-copy
```

The path-adjusted copy intentionally does not retain the original frozen hashes. Scientific definitions, sample rules and result history remain documented in `planning/`, `reports/` and `provenance/`.

## Corrective Stage003F analysis

Stage003F starts from `runs/003_cell_context/site_pair_metrics.tsv`. It retains the Stage003 cell set, 20-cell threshold, same-site pairing and equal-site patient aggregation. It produces:

- patient-level baseline, follow-up and change scores;
- Welch baseline and follow-up group comparisons;
- HC3 baseline-adjusted ANCOVA;
- concurrent inflammation-change sensitivity models;
- all-patient, 10x v3.1-only and released author-paired subsets;
- outcome-by-chemistry and outcome-by-batch tables.

All remitters used 10x v3.1 and the v3 stratum contained only non-remitters. A chemistry-stratified remission permutation is therefore not identifiable. The v3.1-only subset is reported without presenting it as an independent validation.
