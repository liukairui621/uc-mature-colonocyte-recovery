"""Draw Figures 1, 2 and 4 using the common manuscript figure system."""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
from scipy import stats as st

ROOT = Path(__file__).resolve().parents[1]
INK = "#1F2D3D"
MUTED = "#536273"
GRID = "#D9E0E6"
BLUE = "#2878B5"
GREEN = "#2A9D8F"
ORANGE = "#E08A1E"
RED = "#C95A49"
GRAY = "#7C8790"


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9.2,
        "axes.titlesize": 10.3, "axes.labelsize": 9.3,
        "xtick.labelsize": 8.4, "ytick.labelsize": 8.4,
        "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def save(fig, stem: str, destinations: list[Path]) -> None:
    for folder in destinations:
        folder.mkdir(parents=True, exist_ok=True)
        fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.08)
        fig.savefig(folder / f"{stem}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def finish_axis(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.6, alpha=0.65, zorder=0)


def add_panel_labels(fig, axes, labels) -> None:
    fig.canvas.draw()
    for ax, label in zip(axes, labels):
        box = ax.get_position()
        fig.text(box.x0 - 0.025, box.y1 + 0.014, label, ha="left", va="bottom",
                 fontsize=11.5, fontweight="bold", color=INK)


def add_box(ax, y, color, heading, detail) -> None:
    x, w, h = 0.08, 0.84, 0.095
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.012",
                         facecolor=color, edgecolor=color, alpha=0.12, linewidth=1.5)
    ax.add_patch(box)
    ax.text(0.5, y + h * 0.66, heading, ha="center", va="center",
            fontsize=10.1, fontweight="bold", color=INK)
    ax.text(0.5, y + h * 0.27, detail, ha="center", va="center",
            fontsize=8.5, color=MUTED, linespacing=1.12)


def figure1(destinations: list[Path]) -> None:
    fig, ax = plt.subplots(figsize=(6.7, 9.2))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    stages = [
        (0.845, BLUE, "GSE92415 | Exploratory discovery", "65 paired patients | Candidate6 assembled"),
        (0.705, GRAY, "Fixed analysis specification", "Genes, weights, aliases and primary model"),
        (0.565, RED, "GSE23597 | Designated validation", "32 paired patients | Clinical-response criterion not met"),
        (0.425, ORANGE, "GSE73661 | Endoscopic healing", "IFX: positive association | VDZ: imprecise positive estimates"),
        (0.285, GREEN, "GSE282122 | Longitudinal single-cell analysis", "Composition and CT-colonocyte expression analysed separately"),
        (0.145, BLUE, "Recovery depth and transcriptional scope", "Non-IBD controls and three non-overlapping CT programs"),
    ]
    for y, color, heading, detail in stages:
        add_box(ax, y, color, heading, detail)
    for (y1, *_), (y2, *__) in zip(stages[:-1], stages[1:]):
        ax.annotate("", xy=(0.5, y2 + 0.112), xytext=(0.5, y1 - 0.016),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5))
    end = FancyBboxPatch((0.08, 0.025), 0.84, 0.075,
                         boxstyle="round,pad=0.012,rounding_size=0.012",
                         facecolor="#EAF5F2", edgecolor=GREEN, linewidth=1.5)
    ax.add_patch(end)
    ax.text(0.5, 0.062, "Focused CT-colonocyte re-expression with incomplete bulk normalization",
            ha="center", va="center", fontsize=9.3, fontweight="bold", color=INK)
    ax.annotate("", xy=(0.5, 0.105), xytext=(0.5, 0.137),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5))
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.02, top=0.985)
    save(fig, "Figure1_study_design", destinations)


def figure2(destinations: list[Path]) -> None:
    forest = pd.DataFrame([
        ("GSE92415 discovery", .5994069, .1477006, 1.0511132, BLUE),
        ("GSE23597 validation", .4008231, -.4806500, 1.2822963, RED),
        ("GSE73661 IFX", 1.1946910, .5027108, 1.8866713, GREEN),
        ("GSE73661 VDZ W6", .5212629, -.4354414, 1.4779672, ORANGE),
        ("GSE73661 VDZ W12", .8527526, -.3851187, 2.0906239, ORANGE),
    ], columns=["label", "est", "lo", "hi", "color"])
    fig = plt.figure(figsize=(12.0, 5.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.55, .28, 1.0], wspace=0)
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 2])
    y = np.arange(len(forest))[::-1]
    for yy, row in zip(y, forest.itertuples()):
        ax1.plot([row.lo, row.hi], [yy, yy], color=row.color, lw=2.0, zorder=2)
        ax1.scatter(row.est, yy, s=48, color=row.color, edgecolor="white", linewidth=.6, zorder=3)
        ax1.text(2.16, yy, f"{row.est:.3f} ({row.lo:.3f}, {row.hi:.3f})", va="center", fontsize=8.1)
    ax1.axvline(0, color=MUTED, lw=.9, ls="--")
    ax1.set_yticks(y, forest.label); ax1.tick_params(axis="y", length=0, pad=7)
    ax1.set_xlim(-.9, 3.35); ax1.set_ylim(-.6, 4.6)
    ax1.set_xlabel("Adjusted outcome coefficient (95% CI)")
    ax1.set_title("Bulk-cohort outcome associations", loc="left", fontweight="bold", pad=13)
    finish_axis(ax1)

    d = pd.read_csv(ROOT / "runs/002B_GSE73661_support/paired_score_plot_source.tsv", sep="\t")
    out_col = next(c for c in d.columns if "heal" in c.lower())
    score_col = next(c for c in d.columns if "delta" in c.lower() and ("candidate" in c.lower() or "score" in c.lower()))
    mapped = d[out_col].astype(str).str.lower().map(lambda z: "Healed" if z in {"1", "true", "yes", "healed"} else "Not healed")
    rng = np.random.default_rng(73661)
    for j, group in enumerate(["Not healed", "Healed"]):
        vals = d.loc[mapped.eq(group), score_col].astype(float).to_numpy()
        ax2.scatter(j + rng.uniform(-.08, .08, len(vals)), vals, s=34,
                    color=[GRAY, GREEN][j], alpha=.85, edgecolor="white", linewidth=.4, zorder=2)
        mean = vals.mean(); ci = st.t.ppf(.975, len(vals)-1) * vals.std(ddof=1) / math.sqrt(len(vals))
        ax2.errorbar(j, mean, yerr=ci, fmt="_", markersize=18, color=INK, lw=2.1, capsize=5, zorder=3)
    ax2.axhline(0, color=MUTED, lw=.8)
    ax2.set_xticks([0, 1], ["Not healed\n(n=15)", "Healed\n(n=8)"])
    ax2.set_ylabel("Candidate6 change")
    ax2.set_title("Infliximab patient-level changes", loc="left", fontweight="bold", pad=13)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(left=.115, right=.985, bottom=.17, top=.86)
    add_panel_labels(fig, [ax1, ax2], ["A", "B"])
    save(fig, "Figure2_bulk_outcomes", destinations)


def figure4(destinations: list[Path]) -> None:
    visits = pd.read_csv(ROOT / "runs/003F_baseline_adjustment/patient_visit_metrics.tsv", sep="\t")
    visits = visits.loc[visits.scope.eq("all_fixed")].copy()
    visits["group"] = visits.remission.map({"Remission": "Remission", "Non_Remission": "No remission"})
    comp = pd.read_csv(ROOT / "runs/003_cell_context/patient_metrics.tsv", sep="\t")
    comp = comp.loc[(comp.cell_set.eq("exclude_predicted_doublets")) & comp.scope.eq("all_fixed") & comp.threshold.eq(20)].copy()
    comp["group"] = comp.remission.map({"Remission": "Remission", "Non_Remission": "No remission"})
    groups = ["No remission", "Remission"]
    colors = {"No remission": GRAY, "Remission": BLUE}
    rng = np.random.default_rng(282122)
    fig = plt.figure(figsize=(13.2, 4.8))
    gs = fig.add_gridspec(1, 5, width_ratios=[.9, .24, 1.35, .24, .9], wspace=0)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 4])]
    for j, group in enumerate(groups):
        vals = comp.loc[comp.group.eq(group), "delta_ct_logit_epi"].dropna().to_numpy(float)
        axes[0].scatter(j + rng.uniform(-.07, .07, len(vals)), vals, s=31, color=colors[group], alpha=.85,
                        edgecolor="white", linewidth=.4)
        ci = st.t.ppf(.975, len(vals)-1) * vals.std(ddof=1) / math.sqrt(len(vals))
        axes[0].errorbar(j, vals.mean(), yerr=ci, fmt="_", markersize=16, lw=2, color=INK, capsize=5)
    axes[0].axhline(0, color=MUTED, lw=.8, ls="--")
    axes[0].set_xticks([0, 1], ["No remission\n(n=13)", "Remission\n(n=6)"])
    axes[0].set_ylabel("Post - Pre empirical logit")
    axes[0].set_title("CT fraction within epithelium", loc="left", fontweight="bold", pad=13)

    xpos = {("No remission", "Baseline"): 0, ("No remission", "Follow-up"): 1,
            ("Remission", "Baseline"): 3, ("Remission", "Follow-up"): 4}
    for row in visits.itertuples():
        x0, x1 = xpos[(row.group, "Baseline")], xpos[(row.group, "Follow-up")]
        axes[1].plot([x0, x1], [row.pre_ct_state_score, row.post_ct_state_score], color=colors[row.group], alpha=.35, lw=1)
        axes[1].scatter([x0, x1], [row.pre_ct_state_score, row.post_ct_state_score], color=colors[row.group],
                        s=24, edgecolor="white", linewidth=.35, zorder=3)
    axes[1].set_xticks([0, 1, 3, 4], ["Baseline", "Follow-up", "Baseline", "Follow-up"], rotation=15)
    axes[1].text(.5, -.20, "No remission (n=10)", ha="center", transform=axes[1].get_xaxis_transform(), fontsize=8.2)
    axes[1].text(3.5, -.20, "Remission (n=6)", ha="center", transform=axes[1].get_xaxis_transform(), fontsize=8.2)
    axes[1].set_ylabel("Candidate6 CT pseudobulk score")
    axes[1].set_title("Baseline and follow-up CT score", loc="left", fontweight="bold", pad=13)

    for j, group in enumerate(groups):
        vals = visits.loc[visits.group.eq(group), "delta_ct_state_score"].dropna().to_numpy(float)
        axes[2].scatter(j + rng.uniform(-.07, .07, len(vals)), vals, s=31, color=colors[group], alpha=.85,
                        edgecolor="white", linewidth=.4)
        ci = st.t.ppf(.975, len(vals)-1) * vals.std(ddof=1) / math.sqrt(len(vals))
        axes[2].errorbar(j, vals.mean(), yerr=ci, fmt="_", markersize=16, lw=2, color=INK, capsize=5)
    axes[2].axhline(0, color=MUTED, lw=.8, ls="--")
    axes[2].set_xticks([0, 1], ["No remission\n(n=10)", "Remission\n(n=6)"])
    axes[2].set_ylabel("Post - Pre Candidate6 score")
    axes[2].set_title("CT-score change", loc="left", fontweight="bold", pad=13)
    axes[2].text(.03, .97, "Adjusted beta=1.442\nHC3 95% CI 0.121 to 2.764",
                 transform=axes[2].transAxes, ha="left", va="top", fontsize=7.8, color=MUTED)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(left=.065, right=.99, bottom=.25, top=.84)
    add_panel_labels(fig, axes, ["A", "B", "C"])
    save(fig, "Figure4_composition_vs_state", destinations)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", type=Path, default=ROOT / "figures/publication")
    parser.add_argument("--mirror", type=Path, action="append", default=[])
    args = parser.parse_args()
    destinations = [args.dest, *args.mirror]
    style(); figure1(destinations); figure2(destinations); figure4(destinations)


if __name__ == "__main__":
    main()
