"""Draw Figures 3 and 6 using the common manuscript figure system."""
from __future__ import annotations

import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "runs/005_health_function"
INK = "#1F2D3D"; GRID = "#D9E0E6"; BLUE = "#2878B5"; GREEN = "#2A9D8F"; ORANGE = "#E08A1E"


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.2, "axes.titlesize": 10.3,
                         "axes.labelsize": 9.3, "xtick.labelsize": 8.4, "ytick.labelsize": 8.4,
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


def panel_labels(fig, axes):
    fig.canvas.draw()
    for ax, label in zip(axes, "AB"):
        b = ax.get_position()
        fig.text(b.x0 - .025, b.y1 + .015, label, fontsize=11.5, fontweight="bold", color=INK)


def figure3(destinations):
    d = pd.read_csv(DATA / "healthy_reference_contrasts.tsv", sep="\t")
    d = d.loc[d.scope.eq("main") & d.outcome.eq("Yes")]
    order = ["GSE92415_W6", "GSE73661_IFX_W4W6", "GSE73661_VDZ_W6", "GSE73661_VDZ_W12"]
    labels = ["GLM/placebo W6 responders (n=32)", "IFX W4-6 healed (n=8)",
              "VDZ W6 healed (n=6)", "VDZ W12 healed (n=3)"]
    fig = plt.figure(figsize=(7.8, 7.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, .30, 1], hspace=0)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[2, 0])]
    for ax, prog, title, color in zip(axes, ["CANDIDATE6", "COMMON194_INFLAMMATION"],
                                      ["Six-gene colonocyte score", "Inflammatory-response score"],
                                      [BLUE, ORANGE]):
        z = d.loc[d.program.eq(prog)].set_index("context").loc[order]
        y = np.arange(4)[::-1]
        ax.axvline(0, color="#65717D", lw=.9, ls="--")
        ax.errorbar(z.estimate, y, xerr=[z.estimate-z.lower, z.upper-z.estimate], fmt="o",
                    color=color, ecolor=color, capsize=3, ms=5.5, lw=1.8)
        ax.set_yticks(y, labels); ax.tick_params(axis="y", length=0, pad=7)
        ax.set_ylim(-.6, 3.6)
        ax.set_xlabel("Post-treatment minus non-IBD control\n(mean log2-expression score; 95% CI)")
        ax.set_title(title, loc="left", fontweight="bold", pad=13)
        finish(ax)
    fig.subplots_adjust(left=.37, right=.97, bottom=.09, top=.92)
    panel_labels(fig, axes)
    save(fig, "Figure3_healthy_reference", destinations)


def figure6(destinations):
    d = pd.read_csv(DATA / "ct_program_models.tsv", sep="\t")
    d = d.loc[d.model.eq("baseline_ANCOVA")].set_index("program")
    c = pd.read_csv(DATA / "ct_candidate_coordination.tsv", sep="\t").set_index("program")
    order = ["CT_MARKERS_EXCLUDING_CANDIDATE6", "GOBP_INTESTINAL_ABSORPTION", "GOBP_BRUSH_BORDER_ASSEMBLY"]
    labels = ["CT markers excluding Candidate6 (16 genes)", "Intestinal absorption (44 genes)",
              "Brush-border assembly (7 genes)"]
    y = np.arange(3)[::-1]
    fig = plt.figure(figsize=(8.2, 7.3))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, .34, 1], hspace=0)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[2, 0])]
    z = d.loc[order]
    axes[0].axvline(0, color="#65717D", ls="--", lw=.9)
    axes[0].errorbar(z.estimate, y, xerr=[z.estimate-z.lower, z.upper-z.estimate], fmt="o",
                     color=BLUE, ecolor=BLUE, capsize=3, ms=5.5, lw=1.8)
    axes[0].set_yticks(y, labels); axes[0].set_ylim(-.6, 2.6)
    axes[0].set_xlabel("Remission coefficient\n(mean log2-CPM score; HC3 95% CI)")
    axes[0].set_title("Baseline-adjusted remission association", loc="left", fontweight="bold", pad=13)
    axes[0].tick_params(axis="y", length=0, pad=7); finish(axes[0])
    z = c.loc[order]
    axes[1].axvline(0, color="#65717D", ls="--", lw=.9)
    axes[1].errorbar(z.pearson_r, y, xerr=[z.pearson_r-z.lower, z.upper-z.pearson_r], fmt="o",
                     color=GREEN, ecolor=GREEN, capsize=3, ms=5.5, lw=1.8)
    axes[1].set_yticks(y, labels); axes[1].set_ylim(-.6, 2.6); axes[1].set_xlim(-.8, .8)
    axes[1].set_xlabel("Correlation with Candidate6 change\n(Pearson r; Fisher-z 95% CI)")
    axes[1].set_title("Patient-level longitudinal coordination", loc="left", fontweight="bold", pad=13)
    axes[1].tick_params(axis="y", length=0, pad=7); finish(axes[1])
    fig.subplots_adjust(left=.43, right=.96, bottom=.10, top=.91)
    panel_labels(fig, axes)
    save(fig, "Figure6_reference_programs", destinations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", type=Path, default=ROOT / "figures/publication")
    parser.add_argument("--mirror", type=Path, action="append", default=[])
    args = parser.parse_args()
    destinations = [args.dest, *args.mirror]
    style(); figure3(destinations); figure6(destinations)


if __name__ == "__main__":
    main()
