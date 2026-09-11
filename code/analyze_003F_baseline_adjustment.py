#!/usr/bin/env python3
"""Post-result qualification of the Stage003 CT-colonocyte association.

This script preserves the original Stage003 estimand and adds patient-level
baseline, follow-up, ANCOVA, inflammation-change, and technical-subset
diagnostics. It does not redefine the frozen exploratory analysis.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "runs" / "003_cell_context" / "site_pair_metrics.tsv"
ORIGINAL_RESULT = ROOT / "runs" / "003_cell_context" / "joint_primary_patient_models.tsv"
PLAN = ROOT / "planning" / "analysis_plan_003F_baseline_adjustment.json"
OUT = ROOT / "runs" / "003F_baseline_adjustment"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bool_series(x: pd.Series) -> pd.Series:
    return x.astype(str).str.lower().isin({"true", "1", "yes"})


def patient_visits(site: pd.DataFrame, scope: str) -> pd.DataFrame:
    z = site.loc[
        (site["cell_set"] == "exclude_predicted_doublets")
        & (site["scope"] == scope)
        & (site["threshold"] == 20)
        & bool_series(site["state_evaluable"])
    ].copy()
    z["delta_inflammation_score"] = z["post_inflammation_score"] - z["pre_inflammation_score"]
    cols = [
        "pre_ct_state_score",
        "post_ct_state_score",
        "delta_ct_state_score",
        "pre_inflammation_score",
        "post_inflammation_score",
        "delta_inflammation_score",
    ]
    out = z.groupby("patient", as_index=False).agg(
        remission=("remission", "first"),
        library_type=("library_type", "first"),
        batch=("batch", "first"),
        n_eligible_site_pairs=("site", "size"),
        **{c: (c, "mean") for c in cols},
    )
    out["remission_binary"] = (out["remission"] == "Remission").astype(int)
    out["scope"] = scope
    return out


def robust_model(data: pd.DataFrame, outcome: str, covariates: list[str], model: str, scope: str) -> dict:
    d = data[[outcome, *covariates]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    X = sm.add_constant(d[covariates], has_constant="add")
    fit = sm.OLS(d[outcome], X).fit()
    robust = fit.get_robustcov_results(cov_type="HC3", use_t=True)
    names = list(X.columns)
    j = names.index("remission_binary")
    ci = robust.conf_int(alpha=0.05)[j]
    return {
        "scope": scope,
        "model": model,
        "outcome": outcome,
        "formula": f"{outcome} ~ " + " + ".join(covariates),
        "n": int(len(d)),
        "remission": int(d["remission_binary"].sum()),
        "nonremission": int((1 - d["remission_binary"]).sum()),
        "parameters": int(len(names)),
        "residual_df": float(fit.df_resid),
        "estimate": float(robust.params[j]),
        "se_HC3": float(robust.bse[j]),
        "lower_HC3": float(ci[0]),
        "upper_HC3": float(ci[1]),
        "p_HC3": float(robust.pvalues[j]),
        "r_squared": float(fit.rsquared),
    }


def group_description(data: pd.DataFrame, variable: str, scope: str) -> dict:
    rem = data.loc[data["remission_binary"] == 1, variable].dropna().to_numpy(float)
    non = data.loc[data["remission_binary"] == 0, variable].dropna().to_numpy(float)
    test = stats.ttest_ind(rem, non, equal_var=False)
    return {
        "scope": scope,
        "variable": variable,
        "n_remission": int(len(rem)),
        "n_nonremission": int(len(non)),
        "mean_remission": float(np.mean(rem)),
        "sd_remission": float(np.std(rem, ddof=1)),
        "mean_nonremission": float(np.mean(non)),
        "sd_nonremission": float(np.std(non, ddof=1)),
        "mean_difference": float(np.mean(rem) - np.mean(non)),
        "welch_t": float(test.statistic),
        "welch_df": float(test.df),
        "welch_p": float(test.pvalue),
    }


def correlation(data: pd.DataFrame, x: str, y: str, label: str) -> dict:
    d = data[[x, y]].dropna()
    r, p = stats.pearsonr(d[x], d[y])
    return {"analysis": label, "x": x, "y": y, "n": int(len(d)), "pearson_r": float(r), "p": float(p)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    site = pd.read_csv(INPUT, sep="\t", low_memory=False)
    scopes = ["all_fixed", "v31_only", "author_both"]
    patients = pd.concat([patient_visits(site, s) for s in scopes], ignore_index=True)
    patients.to_csv(OUT / "patient_visit_metrics.tsv", sep="\t", index=False)

    descriptions = []
    models = []
    for scope in scopes:
        d = patients.loc[patients["scope"] == scope].copy()
        for variable in ["pre_ct_state_score", "post_ct_state_score", "delta_ct_state_score"]:
            descriptions.append(group_description(d, variable, scope))
        models.append(
            robust_model(
                d,
                "delta_ct_state_score",
                ["remission_binary"],
                "unadjusted_change_HC3",
                scope,
            )
        )
        models.append(
            robust_model(
                d,
                "post_ct_state_score",
                ["remission_binary", "pre_ct_state_score"],
                "baseline_adjusted_ANCOVA",
                scope,
            )
        )
        models.append(
            robust_model(
                d,
                "delta_ct_state_score",
                ["remission_binary", "pre_ct_state_score", "delta_inflammation_score", "pre_inflammation_score"],
                "baseline_candidate_inflammation_change_and_baseline_inflammation",
                scope,
            )
        )
        models.append(
            robust_model(
                d,
                "delta_ct_state_score",
                ["remission_binary", "pre_ct_state_score", "delta_inflammation_score"],
                "baseline_and_inflammation_change",
                scope,
            )
        )

    pd.DataFrame(descriptions).to_csv(OUT / "group_visit_descriptions.tsv", sep="\t", index=False)
    model_table = pd.DataFrame(models)
    model_table.to_csv(OUT / "baseline_adjusted_models.tsv", sep="\t", index=False)

    main_data = patients.loc[patients["scope"] == "all_fixed"].copy()
    correlations = pd.DataFrame(
        [
            correlation(main_data, "pre_ct_state_score", "delta_ct_state_score", "baseline_candidate_vs_candidate_change"),
            correlation(main_data, "delta_ct_state_score", "delta_inflammation_score", "candidate_change_vs_inflammation_change"),
        ]
    )
    correlations.to_csv(OUT / "correlations.tsv", sep="\t", index=False)

    fixed = pd.read_csv(ROOT / "runs" / "003_preexpression" / "fixed_patient_structure.tsv", sep="\t")
    chemistry = (
        fixed.groupby(["remission", "library_type"], observed=True)
        .size()
        .rename("patients")
        .reset_index()
    )
    chemistry.to_csv(OUT / "outcome_by_chemistry.tsv", sep="\t", index=False)
    batch = (
        fixed.groupby(["remission", "library_type", "batch"], observed=True)
        .size()
        .rename("patients")
        .reset_index()
    )
    batch.to_csv(OUT / "outcome_by_chemistry_batch.tsv", sep="\t", index=False)

    original = pd.read_csv(ORIGINAL_RESULT, sep="\t")
    original_estimate = float(original.loc[original["metric"] == "delta_ct_state_score", "estimate"].iloc[0])
    reproduced_estimate = float(
        model_table.loc[
            (model_table["scope"] == "all_fixed") & (model_table["model"] == "unadjusted_change_HC3"),
            "estimate",
        ].iloc[0]
    )
    raw_change_reproduction_error = abs(reproduced_estimate - original_estimate)
    if raw_change_reproduction_error > 1e-12:
        raise RuntimeError(
            f"Patient reconstruction did not reproduce Stage003: error={raw_change_reproduction_error:.3g}"
        )
    if len(main_data) != 16 or int(main_data["remission_binary"].sum()) != 6:
        raise RuntimeError("Unexpected all-fixed CT analysis population")

    summary = {
        "status": "COMPLETE",
        "role": "Post-result corrective sensitivity analysis; not preregistered and not a replacement validation",
        "main_scope": "exclude_predicted_doublets::all_fixed::threshold20",
        "patient_unit": True,
        "all_remitters_v31": bool(
            fixed.loc[fixed["remission"] == "Remission", "library_type"].str.lower().eq("3 prime v3.1").all()
        ),
        "stratified_permutation_identifiable": False,
        "stratified_permutation_reason": "The v3 stratum contains no remitters; report v3.1-only estimates rather than fabricating unsupported cross-stratum label exchanges.",
        "stage003_raw_change_estimate": original_estimate,
        "stage003F_reproduced_raw_change_estimate": reproduced_estimate,
        "raw_change_reproduction_error": raw_change_reproduction_error,
        "population_check": {"patients": int(len(main_data)), "remitters": int(main_data["remission_binary"].sum())},
        "inflammation_adjustment_boundary": "Concurrent inflammation change may be part of mucosal recovery; this model describes shared variation and is not a controlled direct-effect estimate.",
        "inputs": {
            str(INPUT.relative_to(ROOT)): sha256(INPUT),
            str(ORIGINAL_RESULT.relative_to(ROOT)): sha256(ORIGINAL_RESULT),
            str(PLAN.relative_to(ROOT)): sha256(PLAN),
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "sessionInfo.txt").write_text(
        "\n".join([f"Python {platform.python_version()}", f"numpy {np.__version__}", f"pandas {pd.__version__}", f"scipy {scipy.__version__}", f"statsmodels {statsmodels.__version__}"]) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"out": str(OUT), "patients": len(main_data), "models": len(models)}, indent=2))


if __name__ == "__main__":
    main()
