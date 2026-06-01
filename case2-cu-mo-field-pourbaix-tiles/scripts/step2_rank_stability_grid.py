#!/usr/bin/env python3
"""
Step 2 - Rank species stability on the potential-pH grid (field-corrected CHE).

Like the iridium case, this replaces the analytical Pourbaix construction with a
brute-force sweep: lay a fine grid over the potential-pH plane, evaluate the free
energy of every species at every cell, and keep the lowest. No pairwise
equilibria, no pKa's, no boundary algebra - the boundaries emerge on their own.

What is specific to this case is HOW the potential and pH enter. Two effects:

1. Proton-coupled electron transfer (computational hydrogen electrode):
       G(U) = G_system - G_Mo - m*G_H2O + (m - n/2)*G_H2 + (n - 2m)*q*U
   m = water molecules exchanged, n = H atoms, q = electrons per proton (= 1).

2. The local interfacial electric field. A Helmholtz double-layer model maps the
   electrode potential and pH onto the field felt by the adsorbate:
       E(U, pH) = C_H * (U - kT*ln10*pH/e - U_PZC) / (eps * eps0)      [V/m]
   and the field tilts each species' energy through its fitted dipole/polarizability:
       G_field = G(U) + mu*E - (alpha/2)*E**2

The species with the lowest G_field wins the cell. Potential is referenced to RHE.

Inputs : data/field_response_fit.csv     (species, mu, alpha)
         data/formation_terms_cu111.csv  (G_system, G_Mo, G_H2O, G_H2, m, n per species)
Output : results/stability_grid_cu111.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent
FIT_FILE = ROOT / "data" / "field_response_fit.csv"
TERMS_FILE = ROOT / "data" / "formation_terms_cu111.csv"
OUTPUT_FILE = ROOT / "results" / "stability_grid_cu111.csv"

# Double-layer / field model
C_H = 0.25            # Helmholtz capacitance (F/m^2)
U_PZC = -0.70         # potential of zero charge (V)
EPSILON = 2.0         # relative permittivity of the inner layer
EPSILON_0 = 8.854187817e-12  # vacuum permittivity (F/m)
K_B = 1.380649e-23    # Boltzmann constant (J/K)
T = 298.15            # K
E_CHARGE = 1.602176634e-19   # elementary charge (C)
Q = 1.0               # electrons per proton-coupled transfer

# Diagram window (U vs. RHE)
PH_MIN, PH_MAX = 0.0, 14.0
U_MIN, U_MAX = -1.0, 1.0
CELL_PH = 0.01
CELL_U = 0.01


def cell_centres(lo, hi, step):
    edges = np.arange(lo, hi + step, step)
    return edges, (edges[:-1] + edges[1:]) / 2.0


def interfacial_field(U, pH):
    """Local field in V/Angstrom from the Helmholtz double-layer model."""
    E_Vm = C_H * (U - K_B * T * np.log(10) * pH / E_CHARGE - U_PZC) / (EPSILON * EPSILON_0)
    return E_Vm / 1e10


def main() -> None:
    fit = pd.read_csv(FIT_FILE).set_index("species")
    terms = pd.read_csv(TERMS_FILE).set_index("Name")

    species = list(fit.index)
    # Stack per-species parameters (length S), broadcast over the grid later.
    mu = fit["mu"].to_numpy()[:, None, None]
    alpha = fit["alpha"].to_numpy()[:, None, None]
    G_sys = terms.loc["G_system", species].to_numpy(float)[:, None, None]
    G_Mo = terms.loc["G_Mo", species].to_numpy(float)[:, None, None]
    G_H2O = terms.loc["G_H2O", species].to_numpy(float)[:, None, None]
    G_H2 = terms.loc["G_H2", species].to_numpy(float)[:, None, None]
    m = terms.loc["m", species].to_numpy(float)[:, None, None]
    n = terms.loc["n", species].to_numpy(float)[:, None, None]

    ph_edges, ph_c = cell_centres(PH_MIN, PH_MAX, CELL_PH)
    u_edges, u_c = cell_centres(U_MIN, U_MAX, CELL_U)
    PH, U = np.meshgrid(ph_c, u_c, indexing="xy")  # (nU, nPH)

    E_A = interfacial_field(U[None], PH[None])      # (1, nU, nPH)
    G = (G_sys - G_Mo - m * G_H2O + (m - n / 2.0) * G_H2 + (n - 2 * m) * Q * U[None])
    G_field = G + mu * E_A - 0.5 * alpha * E_A**2   # (S, nU, nPH)

    winner = np.argmin(G_field, axis=0)
    G_win = np.take_along_axis(G_field, winner[None], axis=0)[0]
    ids = np.array(species)

    iU, iPH = np.meshgrid(np.arange(len(u_c)), np.arange(len(ph_c)), indexing="ij")
    out = pd.DataFrame({
        "pH_min": ph_edges[iPH].ravel(), "pH_max": ph_edges[iPH + 1].ravel(),
        "pH_center": PH.ravel(),
        "U_min": u_edges[iU].ravel(), "U_max": u_edges[iU + 1].ravel(),
        "U_center": U.ravel(),
        "most_stable": ids[winner].ravel(),
        "G_field_most_stable_eV": G_win.ravel(),
    })
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_FILE, index=False)
    print(f"Grid: {len(ph_c)} pH x {len(u_c)} U = {len(out)} cells, {len(species)} species")
    print("Stability domains (cell count):", out["most_stable"].value_counts().to_dict())
    print(f"Wrote {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
