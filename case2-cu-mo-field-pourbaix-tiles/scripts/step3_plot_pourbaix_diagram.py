#!/usr/bin/env python3
"""
Step 3 - Plot the field-corrected Pourbaix diagram.

Reads the stability grid from step 2 and colors every cell by its most stable
species. Cells sharing a winner form contiguous colored regions: the stability
domains of the (field-corrected) surface Pourbaix diagram. Potential vs. RHE.

Input : results/stability_grid_cu111.csv
Output: results/pourbaix_diagram_cu111.png
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
INPUT_FILE = ROOT / "results" / "stability_grid_cu111.csv"
OUTPUT_FILE = ROOT / "results" / "pourbaix_diagram_cu111.png"

# Diagram window (match step 2).
PH_MIN, PH_MAX = 0.0, 14.0
U_MIN, U_MAX = -1.0, 1.0

# Styling
X_LABEL, Y_LABEL = "pH", "U (V vs. RHE)"
FONT_FAMILY = "Times New Roman"
LABEL_FONTSIZE = 38
TICKLABEL_FONTSIZE = 38
LEGEND_FONTSIZE = 24
X_TICK_STEP, Y_TICK_STEP = 2.0, 0.2
TICK_LENGTH, TICK_WIDTH, TICK_DIRECTION = 7, 1.5, "in"
FIGURE_WIDTH, FIGURE_HEIGHT, FIGURE_DPI = 12, 12, 600

# Fixed color per species (stable legend regardless of which domains appear).
SPECIES_COLORS = {
    "H3Mo":    "#333333",
    "H3MoOH":  "#0096FF",
    "H2MoOH2": "#10ED43",
    "MoOOH3":  "#FF1A1A",
}
# Pretty chemical-formula labels (mathtext subscripts) for the legend.
SPECIES_LABELS = {
    "H3Mo":    r"$H_3Mo$",
    "H3MoOH":  r"$H_3MoOH$",
    "H2MoOH2": r"$H_2MoOH_2$",
    "MoOOH3":  r"$MoO(OH)_3$",
}
FALLBACK_PALETTE = ["#008085", "#00C142", "#0042CA", "#81007F", "#FF8C00"]


def main() -> None:
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = [FONT_FAMILY, "DejaVu Serif", "serif"]
    # Render mathtext (the subscripted formulas) in the regular body font.
    plt.rcParams["mathtext.default"] = "regular"

    df = pd.read_csv(INPUT_FILE)
    species = sorted(df["most_stable"].unique())
    cmap_colors = [SPECIES_COLORS.get(s, FALLBACK_PALETTE[i % len(FALLBACK_PALETTE)])
                   for i, s in enumerate(species)]
    code_to_int = {s: i for i, s in enumerate(species)}

    grid = (df.assign(c=df["most_stable"].map(code_to_int))
              .pivot(index="U_center", columns="pH_center", values="c").sort_index())

    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), dpi=FIGURE_DPI)
    ax.pcolormesh(grid.columns.to_numpy(), grid.index.to_numpy(), grid.to_numpy(),
                  cmap=ListedColormap(cmap_colors), vmin=0, vmax=len(species) - 1,
                  shading="auto")

    order = df["most_stable"].value_counts().index.tolist()
    handles = [mpatches.Patch(facecolor=cmap_colors[code_to_int[s]], edgecolor="black",
                              label=SPECIES_LABELS.get(s, s)) for s in order]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=LEGEND_FONTSIZE, framealpha=0.0)

    ax.set_xlabel(X_LABEL, fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.set_ylabel(Y_LABEL, fontsize=LABEL_FONTSIZE, fontweight="bold")
    yticks = np.round(np.arange(U_MIN, U_MAX + 1e-9, Y_TICK_STEP), 2)
    ax.set_xticks(np.arange(int(PH_MIN), int(PH_MAX) + 1, int(X_TICK_STEP)))
    ax.set_yticks(yticks[(yticks >= U_MIN - 1e-9) & (yticks <= U_MAX + 1e-9)])
    ax.tick_params(axis="both", which="major", labelsize=TICKLABEL_FONTSIZE,
                   length=TICK_LENGTH, width=TICK_WIDTH, direction=TICK_DIRECTION)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    ax.set_xlim(PH_MIN, PH_MAX)
    ax.set_ylim(U_MIN, U_MAX)

    fig.savefig(OUTPUT_FILE, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
