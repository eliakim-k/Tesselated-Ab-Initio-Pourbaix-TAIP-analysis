#!/usr/bin/env python3
"""
Step 3 - Plot the Pourbaix diagram.

Reads the stability grid from step 2 and colors every cell by its most stable
species. Cells sharing a winner form contiguous colored regions: the stability
domains of the Pourbaix (potential-pH) diagram.

The legend uses the short species codes to keep the figure uncluttered; the
code -> cluster-formula key is given in the README.

Input : results/stability_grid_<dataset>.csv
Output: results/pourbaix_diagram_<dataset>.png
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
DATASETS = ["bhattacharyya"]

# Diagram window (match step 2).
PH_MIN, PH_MAX = 0.0, 14.0
U_MIN, U_MAX = 0.8, 2.2

# Styling
X_LABEL, Y_LABEL = "pH", "U (V vs. SHE)"
FONT_FAMILY = "Times New Roman"
LABEL_FONTSIZE = 38
TICKLABEL_FONTSIZE = 38
LEGEND_FONTSIZE = 24
X_TICK_STEP, Y_TICK_STEP = 2.0, 0.2
TICK_LENGTH, TICK_WIDTH, TICK_DIRECTION = 7, 1.5, "in"
FIGURE_WIDTH, FIGURE_HEIGHT, FIGURE_DPI = 12, 12, 600

PALETTE = [
    "#008085", "#00C142", "#0042CA", "#0000FE", "#81007F",
    "#FE2900", "#FF5201", "#FF8C00", "#FFD700", "#00FF00",
    "#00FFFF", "#0000FF", "#8B00FF", "#FF1493", "#FF69B4",
    "#FF6347", "#DC143C", "#8B4513", "#A9A9A9", "#000000",
]


def draw_pourbaix(ax, df):
    """Draw the Pourbaix diagram onto `ax`; return legend handles (largest domain first)."""
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
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), dpi=FIGURE_DPI)
    handles = draw_pourbaix(ax, df)
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=LEGEND_FONTSIZE, framealpha=0.0)
    out = ROOT / "results" / f"pourbaix_diagram_{name}.png"
    fig.savefig(out, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out.name}")


def main():
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = [FONT_FAMILY, "DejaVu Serif", "serif"]
    for name in DATASETS:
        plot_one(name)


if __name__ == "__main__":
    main()
