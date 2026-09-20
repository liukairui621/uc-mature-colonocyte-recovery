#!/usr/bin/env python3
"""Generate the main publication figures with explicit inter-panel gutters."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import scipy.stats as st


ROOT = Path(__file__).resolve().parents[1]
BLUE = "#2F6690"
GREEN = "#2A9D8F"
ORANGE = "#E9C46A"
GRAY = "#6C757D"
RED = "#C8553D"
INK = "#172A3A"
MUTED = "#34495E"


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def save(fig: plt.Figure, name: str, destinations: list[Path]) -> None:
    for dest in destinations:
        dest.mkdir(parents=True, exist_ok=True)
        fig.savefig(dest / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
        fig.savefig(dest / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_box(ax, x, y, w, h, color, heading, detail, heading_size=11.2, detail_size=9.1):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        facecolor=color,
        edgecolor=color,
        alpha=0.13,
        linewidth=1.8,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h * 0.66, heading, ha="center", va="center",
            weight="bold", fontsize=heading_size, color=INK, linespacing=1.10)
    ax.text(x + w / 2, y + h * 0.25, detail, ha="center", va="center",
            fontsize=detail_size, color=MUTED, linespacing=1.12)


def figure1(destinations: list[Path]) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 10.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.985, "Patient-level evidence for colonocyte expression recovery",
            ha="center", va="top", weight="bold", fontsize=13, color=INK)
    x, w, h = 0.08, 0.84, 0.088
    stages = [
        (0.850, BLUE, "GSE92415 | Exploratory cohort | 65 paired patients",
         "Candidate6 derived from longitudinal expression change"),
        (0.715, GRAY, "Fixed analysis",
         "Genes, equal weights, aliases and primary model specified"),
        (0.580, RED, "GSE23597 | Clinical-response validation | 32 pairs",
         "Response coefficient 0.401; 95% CI -0.481 to 1.282"),
        (0.445, ORANGE, "GSE73661 | Endoscopic healing",
         "IFX: positive association | VDZ: estimates positive, CIs spanning zero"),
        (0.310, GREEN, "GSE282122 | Longitudinal CT-colonocyte analysis",
         "Composition and within-compartment expression assessed separately\nLarger six-gene rebound from a lower baseline in remitters"),
        (0.175, BLUE, "Healthy-reference and reference-program extension",
         "Within-study control comparison in bulk biopsies\nThree CT programs excluding all six candidate genes"),
    ]
    for y, color, heading, detail in stages:
        add_box(ax, x, y, w, h, color, heading, detail, heading_size=10.2, detail_size=8.8)
    centers = [s[0] for s in stages]
    for upper, lower in zip(centers[:-1], centers[1:]):
        ax.annotate("", xy=(0.5, lower + h + 0.008), xytext=(0.5, upper - 0.008),
                    arrowprops=dict(arrowstyle="-|>", lw=1.8, color="#465362"))
    conclusion = FancyBboxPatch((0.08, 0.040), 0.84, 0.085,
                                boxstyle="round,pad=0.010,rounding_size=0.012",
                                facecolor="#EAF6F3", edgecolor=GREEN, linewidth=1.5)
    ax.add_patch(conclusion)
    ax.text(0.5, 0.082,
            "The six-gene signature rebounds during healing, with a residual\n"
            "bulk expression deficit after early healing. Broader reference\n"
            "programs show different longitudinal patterns.",
            ha="center", va="center", fontsize=9.5, color=INK, linespacing=1.20)
    ax.annotate("", xy=(0.5, 0.135), xytext=(0.5, 0.175 - 0.008),
                arrowprops=dict(arrowstyle="-|>", lw=1.8, color="#465362"))
    fig.subplots_adjust(left=0.04, right=0.96, bottom=0.025, top=0.97)
    save(fig, "Figure1_study_design", destinations)


def figure2(destinations: list[Path]) -> None:
    forest = pd.DataFrame([
        ["GSE92415 discovery", 0.5994069, 0.1477006, 1.0511132, BLUE],
        ["GSE23597 primary validation", 0.4008231, -0.4806500, 1.2822963, RED],
        ["GSE73661 IFX", 1.1946910, 0.5027108, 1.8866713, GREEN],
        ["GSE73661 VDZ W6", 0.5212629, -0.4354414, 1.4779672, ORANGE],
        ["GSE73661 VDZ W12", 0.8527526, -0.3851187, 2.0906239, ORANGE],
    ], columns=["label", "est", "lo", "hi", "color"])

    # A dedicated blank GridSpec column guarantees a visible gutter after resizing.
    fig = plt.figure(figsize=(13.8, 5.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 0.22, 1.0], wspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 2])
    y = np.arange(len(forest))[::-1]
    for i, row in forest.iterrows():
        yy = y[i]
        ax1.plot([row.lo, row.hi], [yy, yy], color=row.color, lw=2.3)
        ax1.scatter(row.est, yy, s=54, color=row.color, edgecolor="white", zorder=3)
        ax1.text(2.22, yy, f"{row.est:.3f} ({row.lo:.3f}, {row.hi:.3f})",
                 va="center", fontsize=8.5)
    ax1.axvline(0, color="#444", lw=1, ls="--")
    ax1.set_yticks(y, forest.label)
    ax1.set_xlim(-0.9, 3.45); ax1.set_ylim(-0.7, len(forest)-0.3)
    ax1.set_xlabel("Adjusted outcome coefficient (95% CI)")
    ax1.set_title("A  Bulk-cohort outcome associations", loc="left", weight="bold", pad=8)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.text(-0.88, -0.55, "Displayed together; no pooled estimate", fontsize=8, color=GRAY)

    patient_file = ROOT / "runs/002B_GSE73661_support/paired_score_plot_source.tsv"
    if not patient_file.exists():
        candidates = list((ROOT / "runs/002B_figures").glob("*source*.tsv"))
        patient_file = candidates[0] if candidates else None
    plotted = False
    if patient_file and Path(patient_file).exists():
        d = pd.read_csv(patient_file, sep="\t")
        cols = {c.lower(): c for c in d.columns}
        outcome = next((cols[k] for k in cols if "heal" in k), None)
        delta = next((cols[k] for k in cols if "delta" in k and ("candidate" in k or "score" in k)), None)
        if outcome and delta:
            dd = d[[outcome, delta]].dropna().copy()
            dd[outcome] = dd[outcome].astype(str).str.lower().map(
                lambda z: "Healed" if z in {"1", "true", "yes", "healed"} else "Not healed")
            rng = np.random.default_rng(73661)
            for j, group in enumerate(["Not healed", "Healed"]):
                vals = dd.loc[dd[outcome] == group, delta].astype(float).to_numpy()
                ax2.scatter(j + rng.uniform(-.08, .08, len(vals)), vals, s=36, alpha=.82,
                            color=[GRAY, GREEN][j], edgecolor="white", linewidth=.4)
                mean = vals.mean(); se = vals.std(ddof=1) / math.sqrt(len(vals))
                ci = st.t.ppf(.975, len(vals)-1) * se
                ax2.errorbar(j, mean, yerr=ci, fmt="_", markersize=18, lw=2.4,
                             color=INK, capsize=6)
            plotted = True
    if not plotted:
        ax2.text(.5, .5, "Patient-level source table unavailable", ha="center", va="center",
                 transform=ax2.transAxes)
    ax2.axhline(0, color="#999", lw=.8)
    ax2.set_xticks([0, 1], ["Not healed\n(n=15)", "Healed\n(n=8)"])
    ax2.set_ylabel("Candidate6 change (follow-up minus baseline)", labelpad=10)
    ax2.set_title("B  Infliximab patient-level changes", loc="left", weight="bold", pad=8)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.16, top=0.90)
    save(fig, "Figure2_bulk_outcomes", destinations)


def figure3(destinations: list[Path]) -> None:
    visit_path = ROOT / "runs/003F_baseline_adjustment/patient_visit_metrics.tsv"
    if not visit_path.exists():
        visit_path = ROOT / "results/publication_tables/GSE282122_patient_visit_metrics.tsv"
    visits = pd.read_csv(visit_path, sep="\t")
    visits = visits.loc[visits["scope"] == "all_fixed"].copy()
    visits["group"] = visits["remission"].map({"Remission": "Clinical remission", "Non_Remission": "No clinical remission"})
    composition = pd.read_csv(ROOT / "runs/003_cell_context/patient_metrics.tsv", sep="\t")
    composition = composition.loc[
        (composition["cell_set"] == "exclude_predicted_doublets")
        & (composition["scope"] == "all_fixed")
        & (composition["threshold"] == 20)
    ].copy()
    composition["group"] = composition["remission"].map({"Remission": "Clinical remission", "Non_Remission": "No clinical remission"})
    groups = ["No clinical remission", "Clinical remission"]
    colors = {"No clinical remission": GRAY, "Clinical remission": BLUE}
    rng = np.random.default_rng(282122)

    fig = plt.figure(figsize=(15.5, 5.7))
    gs = fig.add_gridspec(1, 5, width_ratios=[0.92, 0.17, 1.35, 0.17, 0.92], wspace=0.02)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 4])]
    for j, group in enumerate(groups):
        vals = composition.loc[composition["group"] == group, "delta_ct_logit_epi"].dropna().to_numpy(float)
        axes[0].scatter(j + rng.uniform(-0.08, 0.08, len(vals)), vals, s=34,
                        color=colors[group], alpha=0.82, edgecolor="white", linewidth=0.4)
        mean = vals.mean(); se = vals.std(ddof=1) / math.sqrt(len(vals)); ci = st.t.ppf(.975, len(vals)-1) * se
        axes[0].errorbar(j, mean, yerr=ci, fmt="_", markersize=16, lw=2.2, color=INK, capsize=5)
    axes[0].axhline(0, color="#888", lw=0.8, ls="--")
    axes[0].set_xticks([0, 1], ["No remission\n(n=13)", "Remission\n(n=6)"])
    axes[0].set_ylabel("Post - Pre empirical logit", labelpad=8)
    axes[0].set_title("A  CT fraction within epithelium", loc="left", weight="bold", pad=8)

    xpos = {("No clinical remission", "Baseline"): 0, ("No clinical remission", "Follow-up"): 1,
            ("Clinical remission", "Baseline"): 3, ("Clinical remission", "Follow-up"): 4}
    for _, row in visits.iterrows():
        group = row["group"]; x0, x1 = xpos[(group, "Baseline")], xpos[(group, "Follow-up")]
        axes[1].plot([x0, x1], [row["pre_ct_state_score"], row["post_ct_state_score"]],
                     color=colors[group], alpha=0.38, lw=1)
        axes[1].scatter([x0, x1], [row["pre_ct_state_score"], row["post_ct_state_score"]],
                        color=colors[group], s=25, edgecolor="white", linewidth=0.35, zorder=3)
    for group in groups:
        subset = visits.loc[visits["group"] == group]
        for visit, col in [("Baseline", "pre_ct_state_score"), ("Follow-up", "post_ct_state_score")]:
            axes[1].scatter(xpos[(group, visit)], subset[col].mean(), marker="D", s=42, color=INK, zorder=4)
    axes[1].set_xticks([0, 1, 3, 4], ["Baseline", "Follow-up", "Baseline", "Follow-up"])
    axes[1].text(0.5, -0.19, "No remission (n=10)", ha="center", va="top",
                 transform=axes[1].get_xaxis_transform(), fontsize=9)
    axes[1].text(3.5, -0.19, "Remission (n=6)", ha="center", va="top",
                 transform=axes[1].get_xaxis_transform(), fontsize=9)
    axes[1].set_ylabel("Candidate6 CT pseudobulk score", labelpad=9)
    axes[1].set_title("B  Baseline and follow-up CT score", loc="left", weight="bold", pad=8)

    for j, group in enumerate(groups):
        vals = visits.loc[visits["group"] == group, "delta_ct_state_score"].dropna().to_numpy(float)
        axes[2].scatter(j + rng.uniform(-0.08, 0.08, len(vals)), vals, s=34,
                        color=colors[group], alpha=0.82, edgecolor="white", linewidth=0.4)
        mean = vals.mean(); se = vals.std(ddof=1) / math.sqrt(len(vals)); ci = st.t.ppf(.975, len(vals)-1) * se
        axes[2].errorbar(j, mean, yerr=ci, fmt="_", markersize=16, lw=2.2, color=INK, capsize=5)
    axes[2].axhline(0, color="#888", lw=0.8, ls="--")
    axes[2].set_xticks([0, 1], ["No remission\n(n=10)", "Remission\n(n=6)"])
    axes[2].set_ylabel("Post - Pre Candidate6 score", labelpad=8)
    axes[2].set_title("C  CT-score change", loc="left", weight="bold", pad=8)
    axes[2].text(0.03, 0.97, "ANCOVA beta=1.442\nHC3 95% CI 0.121 to 2.764",
                 transform=axes[2].transAxes, ha="left", va="top", fontsize=8.2, color=MUTED)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.23, top=0.89)
    save(fig, "Figure3_composition_vs_state", destinations)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", type=Path, default=ROOT / "figures/publication")
    parser.add_argument("--mirror", type=Path, action="append", default=[])
    args = parser.parse_args()
    destinations = [args.dest, *args.mirror]
    style()
    figure1(destinations)
    figure2(destinations)
    figure3(destinations)
    print("PUBLICATION_MAIN_FIGURES_COMPLETE")


if __name__ == "__main__":
    main()
