# Data

DFT energies computed in this work for the iridium-oxide cluster geometries of
Bhattacharyya, Poidevin & Auer (*J. Phys. Chem. C* 2021, 125, 4379–4390,
[10.1021/acs.jpcc.0c10092](https://doi.org/10.1021/acs.jpcc.0c10092)).

Geometries are from the paper; the energies were recomputed here at
**PBE-D3 / def2-TZVP** (electronic energies; zero-point and entropic corrections
neglected). Energies are in eV.

The data lives in two layouts; `scripts/step1_pivot_to_long_format.py` converts
the pivot table into the long format consumed by the analysis.

## `energies_pivot_bhattacharyya.csv` — pivot layout (input to step 1)

Human-friendly: one row per cluster, one column per net charge state `q`
(`q-2_eV` … `q+2_eV`).

| Column                | Meaning                                                   |
|-----------------------|-----------------------------------------------------------|
| `cluster`             | cluster formula, e.g. `Ir_3_O_4_OH_5_H2O_5`               |
| `q-2_eV` … `q+2_eV`   | total energy (eV) at net charge −2 … +2                   |

Formula encoding: `Ir_3` = 3 Ir; `O_n`, `OH_n`, `H2O_n` = counts of oxygen,
hydroxyl, water groups. Hydrogen count is `1·(OH) + 2·(H2O)`.

Only the charge states reported in the study are tabulated; the remaining cells
are left blank and `step 1` skips them.

## `energies_long_format_bhattacharyya.csv` — long format (analysis input)

Machine-friendly: one row per **species** = one (cluster, charge state) pair, the
unit that competes for stability. Produced from the pivot file by step 1.

| Column        | Meaning                                                                  |
|---------------|--------------------------------------------------------------------------|
| `species`     | cluster formula with charge suffix, e.g. `Ir_3_O_4_OH_10_q+0`            |
| `species_id`  | short code: number = cluster ranked by H content (1 = most H-rich), letter = charge (A=−2, B=−1, C=0, D=+1, E=+2) |
| `energy_eV`   | total energy (eV)                                                        |
| `charge`      | net charge `q`                                                           |
| `n_H`         | number of hydrogen atoms                                                 |

The short `species_id` exists so steps 2–3 never have to carry (or match on) long
formula strings.
