# Iridium Pourbaix Tiles

Build a **potential–pH (Pourbaix) diagram** for iridium-oxide cluster catalysts
directly from first-principles energies, by **brute-force grid sampling**
("tiles") rather than by analytically solving for the equilibria between species.

<p align="center">
  <img src="results/pourbaix_diagram_comparison.png" width="860" alt="Pourbaix diagrams: B3LYP/TZVP vs PBE (Bhattacharyya)">
</p>

## Why this approach (the motivation)

The textbook way to draw a Pourbaix diagram is to work out, *by hand*, every
equilibrium between every pair of species: the `pKa` of each
protonation/deprotonation step, the equilibrium potential of each redox step, the
Nernst line for each proton-coupled electron transfer, and then the intersections
of all those lines that bound each stability field. With a handful of species
this is tractable. With dozens of charge and protonation states of a nanoparticle
cluster it becomes a bookkeeping nightmare and a frequent source of error.

**This method needs none of that.** The only physics it uses is a single
closed-form expression for the relative free energy of *one* species as a
function of `U` and `pH`. The procedure is then purely mechanical:

1. lay a fine grid over the potential–pH plane;
2. evaluate that free energy for **every** species at **every** grid cell;
3. colour each cell by whichever species is lowest in energy there.

No pairwise equilibria, no `pKa` values, no equilibrium-potential algebra, no
line-intersection geometry. The phase boundaries are never solved for — they
simply *emerge* as the contours where the lowest-energy species changes. Adding a
new candidate species, charge state, or protonation level means adding a row to a
table, not re-deriving the diagram.

## Scientific background and case study

This workflow reproduces and extends the constant-potential, explicit-pH
framework of:

> Bhattacharyya, K.; Poidevin, C.; Auer, A. A. *Structure and Reactivity of
> IrOₓ Nanoparticles for the Oxygen Evolution Reaction in Electrocatalysis: An
> Electronic Structure Theory Study.* **J. Phys. Chem. C** 2021, *125*,
> 4379–4390. DOI: [10.1021/acs.jpcc.0c10092](https://doi.org/10.1021/acs.jpcc.0c10092)

Iridium oxide (IrOₓ) is the benchmark anode catalyst for the oxygen evolution
reaction (OER) in acidic water electrolysers. Bhattacharyya et al. model the
nanoparticle as a finite molecular cluster (`Ir₃(Oμ₁)₄(OHμ₁)₄(H₂O)₆`, here written
as `Ir_3_O_4_OH_5_H2O_5` and relatives). Under operating conditions the cluster
exchanges **electrons** with the electrode (changing its net charge `q`) and
**protons** with the solvent (changing its hydrogen count `n`). The most stable
state therefore moves across the `U`–`pH` plane: protonated, hydroxyl-rich species
dominate at low potential, and deprotonation/oxidation takes over as the potential
rises toward the ~1.53 V vs. SHE OER onset.

## The free-energy expression

At the centre `(U, pH)` of each cell the relative Gibbs free energy of species *i*
is evaluated against a fixed reference species (Eq. 4 of Bhattacharyya et al.):

```
dG_i(U, pH) = (mu_i - mu_ref)
            + (n_ref - n_i) * MU_H
            - (n_ref - n_i) * kT_ln10 * pH
            - (q_i - q_ref + n_ref - n_i) * (U + 4.28)
```

| Symbol    | Meaning                                                                  |
|-----------|--------------------------------------------------------------------------|
| `mu`      | species energy (eV) — see the note on **which** energy, below            |
| `q`       | net charge (electrons exchanged with the electrode)                      |
| `n`       | number of H atoms (protons exchanged with the solvent)                   |
| `MU_H`    | effective chemical potential of the exchanged hydrogen, −11.20 eV        |
| `kT_ln10` | `k_B·T·ln 10` ≈ 0.0592 eV at 298.15 K (Nernstian pH slope)               |
| `4.28`    | absolute potential of the SHE (V), converting "vs. SHE" to an electron energy reference (Isse & Gennaro) |

Because only *differences* between species enter, the choice of reference shifts
every `dG` by a constant and **does not change the diagram**.

## The result is only as good as the energies you feed it

The grid procedure is exact, but the diagram it draws inherits **every**
approximation in the input energies `mu`. Three choices matter, in increasing
order of impact:

- **Basis set** — e.g. def2-SVP vs def2-TZVP shifts relative energies.
- **Exchange–correlation functional** — a GGA (PBE) and a hybrid (B3LYP) order
  the charge/oxidation states differently.
- **Which energy is used as `mu`** — this is the big one. Feeding *bare electronic
  (SCF) energies* gives a different diagram than feeding *thermally-corrected
  Gibbs free energies* (which add zero-point, thermal and entropic terms). The two
  are not interchangeable: relative stabilities, and therefore domain boundaries,
  move. (Bhattacharyya et al. explicitly *neglect* ZPE and entropy, i.e. they use
  electronic energies as a free-energy proxy.)

The two datasets shipped here make this concrete. The diagrams come from
different functionals **and** different energy definitions, and they look
substantially different (see the comparison figure above) even though the grid
machinery is identical:

| Dataset (`<name>`)      | Level of theory                              | Energy used                         |
|-------------------------|----------------------------------------------|-------------------------------------|
| `b3lyp_tzvp`            | B3LYP / def2-TZVP                             | Gibbs free energies (thermally corrected) |
| `pbe_bhattacharyya`     | PBE-D3 / def2-TZVP (Bhattacharyya et al.)    | electronic energies (ZPE/entropy neglected) |

The takeaway: this tool *constructs* the diagram reliably, but interpreting it
requires stating clearly which functional, basis set, and energy definition the
input came from.

## Data: pivot vs. long format

The energies live in two layouts, and `step 1` converts the first into the second.

- **Pivot** (`data/energies_pivot_<name>.csv`) — the human-friendly layout: one row
  per cluster, one column per charge state. This is how DFT results typically
  arrive (used here by `b3lyp_tzvp`).
- **Long format** (`data/energies_long_format_<name>.csv`) — the machine-friendly
  layout: one row per **species** = one (cluster, charge state) pair, which is the
  unit that competes for stability. `step 1` also assigns every species a short
  **code** (e.g. `2C`) so the later steps never carry long formula strings:
  number = cluster ranked by H content (1 = most H-rich), letter = charge
  (A = −2, B = −1, C = 0, D = +1, E = +2).

A dataset may also be supplied **directly in long format**, in which case it has no
pivot file and `step 1` simply skips it. The `pbe_bhattacharyya` set is provided
this way.

## Repository layout

```
iridium-pourbaix-tiles/
├── data/
│   ├── energies_pivot_b3lyp_tzvp.csv             # raw DFT energies, pivot layout
│   ├── energies_long_format_b3lyp_tzvp.csv       # long format + codes (from step 1)
│   ├── energies_long_format_pbe_bhattacharyya.csv # supplied directly in long format
│   └── README.md                                 # data dictionary + provenance
├── scripts/
│   ├── step1_pivot_to_long_format.py             # pivot -> long format + species codes
│   ├── step2_rank_stability_grid.py              # long format -> winning species per cell
│   └── step3_plot_pourbaix_diagram.py            # stability grid -> diagrams + comparison
├── results/
│   ├── pourbaix_diagram_b3lyp_tzvp.png
│   ├── pourbaix_diagram_pbe_bhattacharyya.png
│   └── pourbaix_diagram_comparison.png           # (stability_grid_*.csv are regenerated, not committed)
├── requirements.txt
├── LICENSE
└── README.md
```

## Usage

```bash
pip install -r requirements.txt

python scripts/step1_pivot_to_long_format.py    # converts pivot datasets -> long format
python scripts/step2_rank_stability_grid.py     # results/stability_grid_*.csv
python scripts/step3_plot_pourbaix_diagram.py   # results/pourbaix_diagram_*.png + comparison
```

To add a dataset, either drop an `energies_pivot_<name>.csv` in `data/` and list
`<name>` in step 1's `DATASETS`, or drop a ready-made
`energies_long_format_<name>.csv` directly. Then list `<name>` in steps 2 and 3.
Potential/pH window, grid resolution, reference species, colours and styling are
all exposed as constants at the top of each script.

## A note on the trimmed B3LYP data

This release contains **no outlier-removal step**: the diagram is plotted directly
from the data. One charge state in the raw B3LYP set — the `q+1` state of
`Ir_3_O_4_OH_5_H2O_5` (species `1D`) — was a verbatim duplicate of that cluster's
`q-1` energy, an obvious data-entry error that spuriously dominated the map. It
has been removed at source (the cell is blank in `energies_pivot_b3lyp_tzvp.csv`),
so it never enters the analysis and the scripts stay free of correction logic. See
`data/README.md`.

## Requirements

Python ≥ 3.9 with `numpy`, `pandas`, `matplotlib` (see `requirements.txt`). Figures
use Times New Roman where available, falling back to a generic serif otherwise.

## License

Released under the MIT License (see `LICENSE`).
