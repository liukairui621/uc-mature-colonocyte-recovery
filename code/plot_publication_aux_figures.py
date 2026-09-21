"""Draw Figure 5 and all supplementary figures in the common visual system."""
from __future__ import annotations

import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INK = "#1F2D3D"; GRID = "#D9E0E6"; BLUE = "#2878B5"; GREEN = "#2A9D8F"; ORANGE = "#E08A1E"; GRAY = "#7C8790"


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.2, "axes.titlesize": 10.3,
                         "axes.labelsize": 9.3, "xtick.labelsize": 8.3, "ytick.labelsize": 8.3,
                         "axes.linewidth": .8, "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig, stem, destinations):
    for out in destinations:
        out.mkdir(parents=True, exist_ok=True)
        fig.savefig(out / f"{stem}.pdf", bbox_inches="tight", pad_inches=.08)
        fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight", pad_inches=.08)
    plt.close(fig)


def finish(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=.6, alpha=.65)


def panel_labels(fig, axes, labels):
    fig.canvas.draw()
    for ax, label in zip(axes, labels):
        b = ax.get_position()
        fig.text(b.x0 - .025, b.y1 + .014, label, fontsize=11.5, fontweight="bold", color=INK)


def forest(ax, labels, est, lo, hi, colors, xlabel, xlim=None):
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="#65717D", ls="--", lw=.9)
    for yy, e, l, h, color in zip(y, est, lo, hi, colors):
        ax.plot([l, h], [yy, yy], color=color, lw=1.8)
        ax.scatter(e, yy, color=color, s=42, edgecolor="white", linewidth=.5, zorder=3)
    ax.set_yticks(y, labels); ax.tick_params(axis="y", length=0, pad=7)
    ax.set_ylim(-.6, len(labels)-.4); ax.set_xlabel(xlabel)
    if xlim: ax.set_xlim(*xlim)
    finish(ax)


def figure5(destinations):
    d = pd.read_csv(ROOT / "runs/003_cell_context/joint_primary_patient_models.tsv", sep="\t").set_index("metric")
    labels = ["Composition", "Within-group expression", "Total epithelial difference"]
    est = [d.loc["composition_contribution", "estimate"], d.loc["within_state_contribution", "estimate"], 1.3417095305]
    lo = [d.loc["composition_contribution", "lower"], d.loc["within_state_contribution", "lower"], .047]
    hi = [d.loc["composition_contribution", "upper"], d.loc["within_state_contribution", "upper"], 2.636]
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    forest(ax, labels, est, lo, hi, [GRAY, GREEN, BLUE],
           "Remission minus non-remission contribution (95% CI)", (-.9, 3.0))
    ax.set_title("Kitagawa decomposition", loc="left", fontweight="bold", pad=12)
    fig.subplots_adjust(left=.30, right=.97, bottom=.22, top=.83)
    save(fig, "Figure5_kitagawa_decomposition", destinations)


def supplement_s1(destinations):
    d = pd.read_csv(ROOT / "figures/stage003/source_tables/figure003C_sensitivity_models.tsv", sep="\t")
    d = d.loc[d.axis.eq("Six-gene score within CT state") & d.lower_HC3.notna()].copy()
    d = d.drop_duplicates(["branch_label"]).head(8)
    labels = [x.replace(", ", "\n") for x in d.branch_label]
    fig, ax = plt.subplots(figsize=(8.4, max(4.8, .55*len(d)+1.6)))
    forest(ax, labels, d.estimate, d.lower_HC3, d.upper_HC3, [BLUE]*len(d),
           "Remission minus non-remission Candidate6 change (HC3 95% CI)")
    ax.set_title("CT-state sensitivity analyses", loc="left", fontweight="bold", pad=12)
    fig.subplots_adjust(left=.40, right=.96, bottom=.12, top=.90)
    save(fig, "Supplementary_Figure_S1_threshold_sensitivity", destinations)


def supplement_s2(destinations):
    d = pd.read_csv(ROOT / "runs/003D_decontx_state_family/decontx_fixed15_state_models.tsv", sep="\t")
    d = d.loc[d.estimate.notna()].sort_values("estimate")
    labels = d.final_analysis.str.replace("Non ileal ", "", regex=False).tolist()
    y = np.arange(len(d))[::-1]
    fig, ax = plt.subplots(figsize=(9.2, max(5.5, .48*len(d)+1.7)))
    ax.axvline(0, color="#65717D", ls="--", lw=.9)
    ax.errorbar(d.raw_estimate, y+.12, xerr=[d.raw_estimate-d.raw_lower, d.raw_upper-d.raw_estimate],
                fmt="o", color=GRAY, ecolor=GRAY, capsize=2.5, ms=4.5, lw=1.3, label="Raw")
    ax.errorbar(d.estimate, y-.12, xerr=[d.estimate-d.lower, d.upper-d.estimate],
                fmt="o", color=BLUE, ecolor=BLUE, capsize=2.5, ms=4.5, lw=1.3, label="DecontX-corrected")
    ax.set_yticks(y, labels); ax.tick_params(axis="y", length=0, pad=7)
    ax.set_xlabel("Remission minus non-remission Candidate6 change (95% CI)")
    ax.set_title("Fixed epithelial-state family", loc="left", fontweight="bold", pad=12)
    ax.legend(frameon=False, loc="upper right")
    finish(ax); fig.subplots_adjust(left=.31, right=.97, bottom=.12, top=.90)
    save(fig, "Supplementary_Figure_S2_DecontX", destinations)


def supplement_s3(destinations):
    neg = pd.read_csv(ROOT / "runs/003E_ambient_negative_controls/negative_control_CT_models.tsv", sep="\t")
    neg = neg.loc[neg.metric.isin(["delta_PTPRC_logCPM", "delta_COL1A1_logCPM"])].copy()
    cont = pd.read_csv(ROOT / "runs/003E_ambient_negative_controls/CT_contamination_patient_models.tsv", sep="\t")
    fig = plt.figure(figsize=(11.8, 4.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, .28, 1], wspace=0)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 2])]
    forest(axes[0], ["PTPRC", "COL1A1"], neg.estimate, neg.lower, neg.upper, [ORANGE, ORANGE],
           "Remission minus non-remission change\n(log2 CPM; 95% CI)")
    axes[0].set_title("CT negative-control transcripts", loc="left", fontweight="bold", pad=12)
    names = {"delta_contamination_mean": "Mean", "delta_contamination_median": "Median", "delta_contamination_q90": "90th percentile"}
    forest(axes[1], [names[x] for x in cont.metric], cont.estimate, cont.lower, cont.upper,
           [GREEN]*len(cont), "Remission minus non-remission change\n(DecontX contamination; 95% CI)")
    axes[1].set_title("Patient-level contamination summaries", loc="left", fontweight="bold", pad=12)
    fig.subplots_adjust(left=.17, right=.98, bottom=.21, top=.82)
    panel_labels(fig, axes, ["A", "B"])
    save(fig, "Supplementary_Figure_S3_ambient_controls", destinations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", type=Path, default=ROOT / "figures/publication")
    parser.add_argument("--mirror", type=Path, action="append", default=[])
    args = parser.parse_args(); destinations = [args.dest, *args.mirror]
    style(); figure5(destinations); supplement_s1(destinations); supplement_s2(destinations); supplement_s3(destinations)


if __name__ == "__main__":
    main()
