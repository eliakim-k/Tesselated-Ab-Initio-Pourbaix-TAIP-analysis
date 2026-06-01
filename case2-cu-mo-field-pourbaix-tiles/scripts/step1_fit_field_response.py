#!/usr/bin/env python3
"""
Step 1 - Fit each species' response to an electric field.

Under an applied field E, the free energy of an adsorbate shifts quadratically:

    dG_rel(E) = mu * E - (alpha / 2) * E**2

where mu is the surface dipole moment and alpha the polarizability of the species
*relative to the bare Mo reference*. We obtain mu and alpha by fitting this form
to a field scan of DFT total energies. They are the only field-specific inputs the
Pourbaix step needs.

The relative energy removes the field response of the metal slab itself:

    dG_rel(E) = [G_species(E) - G_species(0)] - [G_Mo(E) - G_Mo(0)]

Input : data/total_energy_vs_field.csv   (one row per system, one column per field)
Output: data/field_response_fit.csv      (species, mu, alpha, R2)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent
INPUT_FILE = ROOT / "data" / "total_energy_vs_field.csv"
OUTPUT_FILE = ROOT / "data" / "field_response_fit.csv"

# Reference whose field response is subtracted, and systems that are not adsorbates.
REFERENCE = "Mo"
NON_ADSORBATES = ("Cu(111)", "Mo")

# Field window (V/Å) used for the fit.
FIELD_MIN, FIELD_MAX = -0.8, 0.8


def response(E, mu, alpha):
    return mu * E - 0.5 * alpha * E**2


def r_squared(y, y_hat):
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot


def main() -> None:
    df = pd.read_csv(INPUT_FILE)
    field_cols = [c for c in df.columns if c.startswith("Field_")]
    fields = np.array([float(c.split("_")[1]) for c in field_cols])

    energies = {row["System"]: row[field_cols].to_numpy(float) for _, row in df.iterrows()}
    if REFERENCE not in energies:
        raise SystemExit(f"Reference system '{REFERENCE}' not found in {INPUT_FILE}")

    i0 = int(np.where(np.isclose(fields, 0.0))[0][0])  # zero-field index
    d_ref = energies[REFERENCE] - energies[REFERENCE][i0]

    mask = (fields >= FIELD_MIN - 1e-9) & (fields <= FIELD_MAX + 1e-9)
    x = fields[mask]

    rows = []
    for _, row in df.iterrows():
        species = row["System"]
        if species in NON_ADSORBATES:
            continue
        y = (energies[species] - energies[species][i0] - d_ref)[mask]
        (mu, alpha), _ = curve_fit(response, x, y)
        rows.append({"species": species, "mu": mu, "alpha": alpha,
                     "R2": r_squared(y, response(x, mu, alpha))})

    out = pd.DataFrame(rows, columns=["species", "mu", "alpha", "R2"])
    out.to_csv(OUTPUT_FILE, index=False)
    pd.set_option("display.float_format", lambda v: f"{v:.6f}")
    print(out.to_string(index=False))
    print(f"\nWrote {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
