#!/usr/bin/env python3
"""
Step 3 - Plot the Pourbaix diagram(s).

Reads each stability grid from step 2 and colours every cell by its most stable
species. Cells sharing a winner form contiguous coloured regions: the stability
domains of the Pourbaix (potential-pH) diagram.

Produces, for the configured datasets:
  - one diagram per dataset:        results/pourbaix_diagram_<dataset>.png
  - one side-by-side comparison:    results/pourbaix_diagram_comparison.png

Input : results/stability_grid_<dataset>.csv
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent

# Dataset -> panel title (used in the comparison figure).
DATASETS = {
    "b3lyp_tzvp": "B3LYP / def2-TZVP\n(Gibbs free energies)",
    "pbe_bhattacharyya": "PBE-D3 / def2-TZVP\n(Bhattacharyya et al., electronic energies)",
}

# Diagram window (match step 2).
PH_MIN, PH_MAX = 0.0, 14.0
U_MIN, U_MAX = 0.8, 2.2

# Styling
X_LABEL, Y_LABEL = "pH", "U (V vs. SHE)"
FONT_FAMILY = "Times New Roman"
LABEL_FONTSIZE = 34
TICKLABEL_FONTSIZE = 30
LEGEND_FONTSIZE = 20
TITLE_FONTSIZE = 30
X_TICK_STEP, Y_TICK_STEP = 2.0, 0.2
TICK_LENGTH, TICK_WIDTH, TICK_DIRECTION = 7, 1.5, "in"
FIGURE_DPI = 600

PALETTE = [
    "#008085", "#00C142", "#0042CA", "#0000FE", "#81007F",
    "#FE2900", "#FF5201", "#FF8C00", "#FFD700", "#00FF00",
    "#00FFFF", "#0000FF", "#8B00FF", "#FF1493", "#FF69B4",
    "#FF6347", "#DC143C", "#8B4513", "#A9A9A9", "#000000",
]


def draw_pourbaix(ax, df):
    """Draw one Pourbaix diagram onto `ax`; return legend handles (largest domain first)."""
    species = sorted(df["most_stable"].unique())
    code_to_int = {s: i for i, s in enumerate(species)}
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(species))]

    grid = (df.assign(c=df["most_stable"].map(code_to_int))
              .pivot(index="U_center", columns="pH_center", values="c").sort_index())
    ax.pcolormesh(grid.columns.to_numpy(), grid.index.to_numpy(), grid.to_numpy(),
                  cmap=ListedColormap(colors), vmin=0, vmax=len(species) - 1, shading="auto")

    ax.set_xlabel(X_LABEL, fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.set_ylabel(Y_LABEL, fontsize=LABEL_FONTSIZE, fontweight="bold")
    yticks = np.round(np.arange(U_MIN, U_MAX + 1e-9, Y_TICK_STEP), 2)
    ax.set_xticks(np.arange(int(PH_MIN), int(PH_MAX) + 1, int(X_TICK_STEP)))
    ax.set_yticks(yticks[yticks <= U_MAX + 1e-9])
    ax.tick_params(axis="both", which="major", labelsize=TICKLABEL_FONTSIZE,
                   length=TICK_LENGTH, width=TICK_WIDTH, direction=TICK_DIRECTION)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    ax.set_xlim(PH_MIN, PH_MAX)
    ax.set_ylim(U_MIN, U_MAX)

    order = df["most_stable"].value_counts().index.tolist()
    return [mpatches.Patch(facecolor=colors[code_to_int[s]], edgecolor="black", label=s)
            for s in order]


def plot_one(name):
    df = pd.read_csv(ROOT / "results" / f"stability_grid_{name}.csv")
    fig, ax = plt.subplots(figsize=(12, 12), dpi=FIGURE_DPI)
    handles = draw_pourbaix(ax, df)
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=LEGEND_FONTSIZE, framealpha=0.0)
    out = ROOT / "results" / f"pourbaix_diagram_{name}.png"
    fig.savefig(out, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out.name}")


def plot_comparison():
    names = list(DATASETS)
    fig, axes = plt.subplots(1, len(names), figsize=(11 * len(names), 11), dpi=FIGURE_DPI)
    for ax, name in zip(np.atleast_1d(axes), names):
        df = pd.read_csv(ROOT / "results" / f"stability_grid_{name}.csv")
        handles = draw_pourbaix(ax, df)
        ax.set_title(DATASETS[name], fontsize=TITLE_FONTSIZE, fontweight="bold", pad=14)
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
                  fontsize=LEGEND_FONTSIZE, framealpha=0.0)
    fig.subplots_adjust(wspace=0.55)
    out = ROOT / "results" / "pourbaix_diagram_comparison.png"
    fig.savefig(out, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out.name}")


def main():
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = [FONT_FAMILY, "DejaVu Serif", "serif"]
    for name in DATASETS:
        plot_one(name)
    plot_comparison()


if __name__ == "__main__":
    main()
