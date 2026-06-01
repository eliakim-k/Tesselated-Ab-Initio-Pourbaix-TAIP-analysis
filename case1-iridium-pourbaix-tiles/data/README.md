# Data

DFT energies for the iridium-oxide cluster models of Bhattacharyya, Poidevin &
Auer (*J. Phys. Chem. C* 2021, 125, 4379–4390,
[10.1021/acs.jpcc.0c10092](https://doi.org/10.1021/acs.jpcc.0c10092)).

Energies use two layouts. `scripts/step1_pivot_to_long_format.py` converts a
**pivot** table into the **long format** consumed by the analysis. A dataset can
also be supplied directly in long format, in which case it has no pivot file.

## Datasets

| `<name>`              | Level of theory                            | Energy stored                              | Supplied as |
|-----------------------|--------------------------------------------|--------------------------------------------|-------------|
| `b3lyp_tzvp`          | B3LYP / def2-TZVP                          | Gibbs free energies (thermally corrected)  | pivot → long (step 1) |
| `pbe_bhattacharyya`   | PBE-D3 / def2-TZVP (Bhattacharyya et al.)  | electronic energies (ZPE/entropy neglected)| long format (direct) |

The two sets use different functionals **and** different energy definitions, so
the resulting diagrams differ (see the root `README.md`). The `pbe_bhattacharyya`
set is the reference species list from the study, provided ready-made in long
format.

## `energies_pivot_<name>.csv` — pivot layout (input to step 1)

Human-friendly: one row per cluster, one column per net charge state `q`. Energies
in eV (`q-2_eV` … `q+2_eV`); Hartree columns (`*_Eh`), where present, are ignored.

| Column                | Meaning                                                   |
|-----------------------|-----------------------------------------------------------|
| `cluster`             | cluster formula, e.g. `Ir_3_O_4_OH_5_H2O_5`               |
| `q-2_eV` … `q+2_eV`   | total energy (eV) at net charge −2 … +2                   |

Formula encoding: `Ir_3` = 3 Ir; `O_n`, `OH_n`, `H2O_n` = counts of oxygen,
hydroxyl, water groups. Hydrogen count is `1·(OH) + 2·(H2O)`.

Only `b3lyp_tzvp` ships a pivot file.

## `energies_long_format_<name>.csv` — long format (analysis input)

Machine-friendly: one row per **species** = one (cluster, charge state) pair, the
unit that competes for stability. Produced from a pivot file by step 1, or
supplied directly.

| Column        | Meaning                                                                  |
|---------------|--------------------------------------------------------------------------|
| `species`     | cluster formula with charge suffix, e.g. `Ir_3_O_4_OH_10_q+0`            |
| `species_id`  | short code: number = cluster ranked by H content (1 = most H-rich), letter = charge (A=−2, B=−1, C=0, D=+1, E=+2) |
| `energy_eV`   | total energy (eV)                                                        |
| `charge`      | net charge `q`                                                           |
| `n_H`         | number of hydrogen atoms                                                 |

The short `species_id` exists so steps 2–3 never have to carry (or match on) long
formula strings.

## The removed B3LYP data point

In the original B3LYP spreadsheet the `q+1` energy of `Ir_3_O_4_OH_5_H2O_5` was a
verbatim copy of that cluster's `q-1` energy (identical to all decimals in both
Hartree and eV) — a copy/paste error, not a physical result. Left in, that single
charge state (species `1D`) was over-stabilised and spuriously won ~96 % of the
potential–pH plane.

The fix is applied at source: the `q+1` cell for that cluster is left blank in
`energies_pivot_b3lyp_tzvp.csv`, exactly like other genuinely missing charge
states (e.g. the `q+2` state of `Ir_3_O_4_O_10`). `step 1` skips blank energies, so
species `1D` never appears in the long-format table and the analysis needs no
outlier-removal logic. The other four charge states of `Ir_3_O_4_OH_5_H2O_5`
(`1A`, `1B`, `1C`, `1E`) are valid and retained. The `pbe_bhattacharyya` set is
used as provided.
