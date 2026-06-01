# TAIP map for Iridium oxide clusters

Construct a **potential–pH (Pourbaix) diagram** from first-principles energies by
sampling the potential–pH plane on a grid of tiles and, in each tile, keeping the
species with the lowest free energy.

## The problem

A Pourbaix diagram shows which species is most stable as a function of electrode
potential `U` and `pH`. Built analytically, it requires — for every pair of
candidate species — writing the acid–base (`pKa`) and redox (Nernst) equilibrium
conditions and intersecting the resulting lines to bound each stability field. The
number of pairwise relations grows with the number of charge and protonation
states considered.

## The approach

Each species *i* has a closed-form relative free energy as a function of `U`
(V vs. SHE) and `pH` (Eq. 4 of Bhattacharyya et al.):

```
dG_i(U, pH) = (mu_i - mu_ref)
            + (n_ref - n_i) * MU_H
            - (n_ref - n_i) * kT_ln10 * pH
            - (q_i - q_ref + n_ref - n_i) * (U + 4.28)
```

| Symbol    | Meaning                                                                  |
|-----------|--------------------------------------------------------------------------|
| `mu`      | species energy (eV)                                                      |
| `q`       | net charge (electrons exchanged with the electrode)                      |
| `n`       | number of H atoms (protons exchanged with the solvent)                   |
| `MU_H`    | effective chemical potential of the exchanged hydrogen, −11.20 eV        |
| `kT_ln10` | `k_B·T·ln 10` ≈ 0.0592 eV at 298.15 K (Nernstian pH slope)               |
| `4.28`    | absolute potential of the SHE (V), converting "vs. SHE" to an electron energy reference (Isse & Gennaro) |

Because only *differences* between species enter, the choice of reference species
shifts every `dG` by a constant and leaves the lowest-energy species in each tile
unchanged.

The diagram then follows in three steps:

1. lay a fine grid over the potential–pH plane;
2. evaluate `dG_i` for every species at the center of each tile;
3. assign each tile to the species with the lowest `dG` and color it accordingly.

A phase boundary is simply the locus where the lowest-energy species changes; it is
read off the grid rather than solved for. Since the only physical inputs are the
energies `mu`, the diagram inherits their level of theory — the functional, the
basis set, and whether `mu` are bare electronic or thermally-corrected Gibbs free
energies.

## Application: IrOₓ clusters for the oxygen evolution reaction

We apply the workflow to the iridium-oxide catalyst of:

> Bhattacharyya, K.; Poidevin, C.; Auer, A. A. *Structure and Reactivity of
> IrOₓ Nanoparticles for the Oxygen Evolution Reaction in Electrocatalysis: An
> Electronic Structure Theory Study.* **J. Phys. Chem. C** 2021, *125*,
> 4379–4390. DOI: [10.1021/acs.jpcc.0c10092](https://doi.org/10.1021/acs.jpcc.0c10092)

The nanoparticle is modeled as a finite molecular cluster (`Ir₃(Oμ₁)₄(OHμ₁)₄(H₂O)₆`,
written here as `Ir_3_O_4_OH_5_H2O_5` and relatives). The cluster geometries are
taken from Bhattacharyya et al.; the energies used here were computed in this work
with molecular (Gaussian-basis) DFT (**ORCA**, PBE-D3/def2-TZVP). Each species carries a net
charge `q` and a hydrogen count `n`; the input table lists the charge/protonation
states considered.

<p align="center">
  <img src="results/pourbaix_diagram_bhattacharyya.png" width="560" alt="Pourbaix diagram of IrOx clusters from the tiles method">
</p>

The legend uses short species **codes** to keep the figure uncluttered. Each code
decodes to a cluster as follows (the implicit Ir₃O₄ core is dropped, matching the
convention of the published labels), with the corresponding region number in
Bhattacharyya et al.:

| Code | Cluster (this work)            | Bhattacharyya region |
|------|--------------------------------|----------------------|
| 1D   | Ir₃(OH)₅(H₂O)₅⁺                 | (1)                  |
| 1E   | Ir₃(OH)₅(H₂O)₅²⁺                | (2)                  |
| 2D   | Ir₃(OH)₆(H₂O)₄⁺                 | (3)                  |
| 3C   | Ir₃(OH)₇(H₂O)₃                  | (4)                  |
| 4B   | Ir₃(OH)₁₀⁻                      | (5)                  |
| 5B   | Ir₃(OH)₃O₇⁻                     | (7)                  |
| 6B   | Ir₃(OH)₂O₈⁻                     | (8)                  |
| 4C   | Ir₃(OH)₁₀ (neutral)            | — (not in the published set) |

(Code = cluster ranked by H content + charge letter A=−2, B=−1, C=0, D=+1, E=+2;
see [`data/README.md`](data/README.md). The published labels use the μ₁ on-top
notation, e.g. our Ir₃(OH)₅ ≡ Ir₃(OHμ₁)₅.)

## Comparison with the published diagram

The published potential–pH diagram of Bhattacharyya et al. (numbered regions 1–8)
is available from the publisher
([figure](https://pubs.acs.org/cms/10.1021/acs.jpcc.0c10092/asset/images/large/jp0c10092_0004.jpeg);
[paper](https://doi.org/10.1021/acs.jpcc.0c10092)). Mapping the most-stable species
of the grid above onto those regions:

- Codes `1D`, `1E`, `2D`, `3C`, `4B`, and `6B` occupy the same regions as published
  (1), (2), (3), (4), (5), and (8), respectively.
- Code `5B` (region 7) is present in the input table but is not the lowest-energy
  species in any tile of this map.
- Published region (6), `Ir₃(OH)₆O₄⁻`, is not represented because that species is
  not in the input table.
- A thin neutral `Ir₃(OH)₁₀` (code `4C`) appears along the (4)/(5) boundary.

## Data: pivot vs. long format

The energies live in two layouts, and `step 1` converts between them.

- **Pivot** (`data/energies_pivot_bhattacharyya.csv`) — one row per cluster, one
  column per charge state. This is how DFT results typically arrive.
- **Long format** (`data/energies_long_format_bhattacharyya.csv`) — one row per
  **species** = one (cluster, charge state) pair, the unit ranked in step 2. `step
  1` also assigns every species a short **code** (e.g. `2D`) so the later steps
  carry no long formula strings: number = cluster ranked by H content (1 = most
  H-rich), letter = charge (A = −2, B = −1, C = 0, D = +1, E = +2).

## Repository layout

```
iridium-pourbaix-tiles/
├── data/
│   ├── energies_pivot_bhattacharyya.csv          # DFT energies, pivot layout
│   ├── energies_long_format_bhattacharyya.csv    # long format + codes (from step 1)
│   └── README.md                                 # data dictionary + provenance
├── scripts/
│   ├── step1_pivot_to_long_format.py             # pivot -> long format + species codes
│   ├── step2_rank_stability_grid.py              # long format -> winning species per tile
│   └── step3_plot_pourbaix_diagram.py            # stability grid -> Pourbaix diagram
├── results/
│   └── pourbaix_diagram_bhattacharyya.png        # (stability_grid_*.csv is regenerated, not committed)
├── requirements.txt
├── LICENSE
└── README.md
```

## Usage

```bash
pip install -r requirements.txt

python scripts/step1_pivot_to_long_format.py    # data/energies_long_format_bhattacharyya.csv
python scripts/step2_rank_stability_grid.py     # results/stability_grid_bhattacharyya.csv
python scripts/step3_plot_pourbaix_diagram.py   # results/pourbaix_diagram_bhattacharyya.png
```

Each script has a configuration block at the top — potential/pH window, grid
resolution, reference species, colors, and styling — so the diagram can be re-tuned
without touching the logic. To analyze a different system, drop an
`energies_pivot_<name>.csv` in `data/` and set `<name>` in the scripts.

## Citation

If you use this method or workflow, please cite:

> Kambale, E. M.; Rivera Rocabado, D. S.; Kanematsu, Y.; Ishimoto, T.
> *Field-Dependent Redox Thermodynamics of MoOₘHₙ Species on Cu(111) and Ni(111)
> Surfaces under Alkaline Hydrogen Evolution Conditions.* Preprints.org, 2026.
> DOI: [10.20944/preprints202604.0944.v1](https://doi.org/10.20944/preprints202604.0944.v1)

The cluster geometries are from Bhattacharyya, Poidevin & Auer, *J. Phys. Chem. C*
2021, 125, 4379 ([10.1021/acs.jpcc.0c10092](https://doi.org/10.1021/acs.jpcc.0c10092));
the energies used here were computed in this work.

## License

Released under the MIT License (see `LICENSE`).
