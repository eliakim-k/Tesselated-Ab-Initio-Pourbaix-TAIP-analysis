#!/usr/bin/env python3
"""
Step 1 - Convert the pivot table to long format (and codify the species).

DFT energies arrive as a PIVOT table: one row per cluster, one column per charge
state. That layout reads well for a human but is awkward for the analysis, where
the quantity that actually competes for stability is a single (cluster, charge
state) pair - a "species" - not a whole row.

This script reshapes the pivot table into LONG format: one row per species. It
also assigns every species a short CODE so the later steps never have to carry
long formula strings around:

    species_id = <cluster number><charge letter>
      cluster number : cluster ranked by hydrogen content (1 = most H-rich)
      charge letter  : A=-2, B=-1, C=0, D=+1, E=+2

Charge states with a blank/non-finite energy in the pivot file are skipped (the
curated dataset only tabulates the charge states reported in the study).

Input : data/energies_pivot_<dataset>.csv
Output: data/energies_long_format_<dataset>.csv
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent.parent
DATASETS = ["bhattacharyya"]

CHARGE_TO_LETTER = {-2: "A", -1: "B", 0: "C", 1: "D", 2: "E"}
ENERGY_UNIT_TAG = "eV"  # only the eV energy columns are used


def count_hydrogens(formula: str) -> int:
    """Count H atoms in a cluster formula: H_<k>, OH_<k> (1 each), H2O_<k> (2 each)."""
    n_H = 0
    for m in re.finditer(r"(?:^|_)H_(\d+)", formula):
        n_H += int(m.group(1))
    for m in re.finditer(r"OH_(\d+)", formula):
        n_H += int(m.group(1))
    for m in re.finditer(r"H2O_(\d+)", formula):
        n_H += 2 * int(m.group(1))
    return n_H


def charge_from_column(col: str) -> int:
    m = re.search(r"(-?\d+)", col)
    if not m:
        raise ValueError(f"Could not read a charge from column '{col}'")
    return int(m.group(1))


def pivot_to_long(pivot_csv: Path, long_csv: Path) -> None:
    df = pd.read_csv(pivot_csv)
    formula_col = df.columns[0]

    energy_cols = [c for c in df.columns[1:] if ENERGY_UNIT_TAG in c]
    if not energy_cols:
        raise SystemExit(f"No '{ENERGY_UNIT_TAG}' energy columns found in {pivot_csv}")
    charges = [charge_from_column(c) for c in energy_cols]

    # Rank clusters by hydrogen content (most H-rich first) for the code number.
    clusters = [(str(r[formula_col]), count_hydrogens(str(r[formula_col])), r)
                for _, r in df.iterrows()]
    clusters.sort(key=lambda x: x[1], reverse=True)

    rows, skipped = [], []
    for cluster_rank, (formula, n_H, row) in enumerate(clusters, start=1):
        for col, q in zip(energy_cols, charges):
            energy = pd.to_numeric(row[col], errors="coerce")
            species_id = f"{cluster_rank}{CHARGE_TO_LETTER.get(q, '?')}"
            if not np.isfinite(energy):
                skipped.append(f"{formula}_q{q:+d} ({species_id})")
                continue
            rows.append({"species": f"{formula}_q{q:+d}", "species_id": species_id,
                         "energy_eV": float(energy), "charge": q, "n_H": n_H})

    out = pd.DataFrame(rows, columns=["species", "species_id", "energy_eV", "charge", "n_H"])
    out.to_csv(long_csv, index=False)
    print(f"[{pivot_csv.name}] -> {long_csv.name}: {len(out)} species")
    if skipped:
        print(f"    skipped {len(skipped)} charge state(s) not tabulated in this dataset")


def main() -> None:
    for name in DATASETS:
        pivot = ROOT / "data" / f"energies_pivot_{name}.csv"
        if not pivot.exists():
            print(f"[{name}] no pivot file found; skipping.")
            continue
        pivot_to_long(pivot, ROOT / "data" / f"energies_long_format_{name}.csv")


if __name__ == "__main__":
    main()
