#!/usr/bin/env python3
"""
Step 2 - Rank species stability on the potential-pH grid.

This is the heart of the method, and where the "tiles" idea replaces the usual
analytical Pourbaix construction. Instead of solving, pair by pair, for the
equilibrium lines between species (the pKa of every protonation step, the
equilibrium potential of every redox step, and the intersections of all those
Nernst lines), we simply:

    1. lay a fine grid over the potential-pH plane,
    2. evaluate the relative free energy of EVERY species at each grid cell, and
    3. keep whichever species has the lowest free energy there.

No equilibria, no pKa's, no boundary algebra: the phase boundaries appear on
their own as the lines where the lowest-energy species changes.

Relative free energy of species i versus a fixed reference species, as a
function of electrode potential U (V vs. SHE) and pH (Eq. 4 of Bhattacharyya
et al., J. Phys. Chem. C 2021, 125, 4379):

    dG_i(U, pH) = (mu_i - mu_ref)
                + (n_ref - n_i) * MU_H
                - (n_ref - n_i) * kT_ln10 * pH
                - (q_i - q_ref + n_ref - n_i) * (U + SHE_VACUUM)

  mu : species energy (eV); q : net charge; n : hydrogen count.

Only free-energy DIFFERENCES decide the winner, so the choice of reference shifts
every dG by a constant and does not change the diagram.

Input : data/energies_long_format_<dataset>.csv
Output: results/stability_grid_<dataset>.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent
DATASETS = {"bhattacharyya": {"reference": "4C"}}

# Thermodynamic parameters
MU_H = -11.20         # effective chemical potential of exchanged hydrogen (eV)
SHE_VACUUM = 4.28     # absolute potential of the SHE (V)
TEMPERATURE = 298.15  # K
K_B = 8.617333262145e-5  # eV/K

# Diagram window
PH_MIN, PH_MAX = 0.0, 14.0
U_MIN, U_MAX = 0.8, 2.2

# Cell size (0.04 pH x 0.004 V -> 350 x 350 visually square cells).
CELL_PH = 0.04
CELL_U = 0.004


def cell_centers(lo, hi, step):
    edges = np.arange(lo, hi + step, step)
    return edges, (edges[:-1] + edges[1:]) / 2.0


def rank_dataset(long_csv: Path, grid_csv: Path, reference_id: str) -> None:
    df = pd.read_csv(long_csv)

    ref_row = df.loc[df["species_id"] == reference_id]
    if ref_row.empty:
        reference_id = df["species_id"].iloc[0]
        ref_row = df.iloc[[0]]
        print(f"    reference not found; falling back to '{reference_id}'")
    ref = ref_row.iloc[0]
    mu_ref, q_ref, n_ref = ref["energy_eV"], ref["charge"], ref["n_H"]

    kT_ln10 = K_B * TEMPERATURE * np.log(10)

    ids = df["species_id"].to_numpy()
    mu = df["energy_eV"].to_numpy()[:, None, None]
    q = df["charge"].to_numpy()[:, None, None]
    n = df["n_H"].to_numpy()[:, None, None]

    ph_edges, ph_c = cell_centers(PH_MIN, PH_MAX, CELL_PH)
    u_edges, u_c = cell_centers(U_MIN, U_MAX, CELL_U)
    PH, U = np.meshgrid(ph_c, u_c, indexing="xy")  # (nU, nPH)

    dn = n_ref - n
    dG = ((mu - mu_ref) + dn * MU_H - dn * kT_ln10 * PH[None]
          - (q - q_ref + dn) * (U[None] + SHE_VACUUM))

    winner = np.argmin(dG, axis=0)
    dG_winner = np.take_along_axis(dG, winner[None], axis=0)[0]

    iU, iPH = np.meshgrid(np.arange(len(u_c)), np.arange(len(ph_c)), indexing="ij")
    out = pd.DataFrame({
        "pH_min": ph_edges[iPH].ravel(), "pH_max": ph_edges[iPH + 1].ravel(),
        "pH_center": PH.ravel(),
        "U_min": u_edges[iU].ravel(), "U_max": u_edges[iU + 1].ravel(),
        "U_center": U.ravel(),
        "most_stable": ids[winner].ravel(),
        "dG_most_stable_eV": dG_winner.ravel(),
    })
    grid_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(grid_csv, index=False)
    print(f"[{long_csv.name}] -> {grid_csv.name}: "
          f"{len(ph_c)}x{len(u_c)} cells, {len(ids)} species, ref {reference_id}")
    print("    domains:", out["most_stable"].value_counts().to_dict())


def main() -> None:
    for name, cfg in DATASETS.items():
        rank_dataset(ROOT / "data" / f"energies_long_format_{name}.csv",
                     ROOT / "results" / f"stability_grid_{name}.csv",
                     cfg["reference"])


if __name__ == "__main__":
    main()
