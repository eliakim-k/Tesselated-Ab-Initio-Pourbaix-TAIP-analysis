# Cu–Mo Field Pourbaix Tiles

Build a **field-corrected surface Pourbaix (potential–pH) diagram** for a
supported catalyst directly from first-principles energies, by **brute-force grid
sampling** ("tiles") rather than by analytically solving for the equilibria
between species.

This is the field-aware sibling of the
[iridium-pourbaix-tiles](../iridium-pourbaix-tiles) workflow: the same
tile-and-rank idea, but each species now also responds to the local interfacial
electric field through a fitted dipole and polarizability.

<p align="center">
  <img src="results/pourbaix_diagram_cu111.png" width="560" alt="Field-corrected Pourbaix diagram of MoOmHn on Cu(111)">
</p>

## Why this approach (the motivation)

The textbook way to draw a Pourbaix diagram is to work out, *by hand*, every
equilibrium between every pair of species: the `pKa` of each
protonation/deprotonation step, the equilibrium potential of each redox step, the
Nernst line for each proton-coupled electron transfer, and the intersections of
all those lines. This method needs none of that. The only physics it uses is a
single closed-form expression for the free energy of *one* species as a function
of `U` and `pH`. The procedure is then mechanical:

1. lay a fine grid over the potential–pH plane;
2. evaluate that free energy for **every** species at **every** grid cell;
3. colour each cell by whichever species is lowest in energy there.

No pairwise equilibria, no `pKa` values, no boundary algebra — the phase
boundaries emerge on their own as the contours where the lowest-energy species
changes.

## Case study: MoOₘHₙ on Cu(111) under alkaline HER

We apply the method to the molybdenum hydride/oxyhydride species that form on a
copper support under hydrogen-evolution conditions:

> Kambale, E. M.; Rivera Rocabado, D. S.; Kanematsu, Y.; Ishimoto, T.
> *Field-Dependent Redox Thermodynamics of MoOₘHₙ Species on Cu(111) and Ni(111)
> Surfaces under Alkaline Hydrogen Evolution Conditions.* Preprints.org, 2026.
> DOI: [10.20944/preprints202604.0944.v1](https://doi.org/10.20944/preprints202604.0944.v1)

The geometries and energies are from periodic plane-wave DFT calculations
(**VASP**). Four adsorbates compete across the `U`–`pH` plane — `H3Mo`, `H3MoOH`,
`H2MoOH2` and `MoOOH3` — differing in how many water molecules (`m`) and hydrogens
(`n`) they carry. As the potential and pH change, the surface switches between the
reduced hydride (`H3Mo`) at low potential and progressively oxidised oxyhydrides
toward higher potential.

### Surface species

The four adsorbed structures on Cu(111) (side and top views):

<table>
<tr>
  <td align="center"><img src="structures/H3Mo.png" width="200"><br><b>H₃Mo</b></td>
  <td align="center"><img src="structures/H3MoOH.png" width="200"><br><b>H₃MoOH</b></td>
</tr>
<tr>
  <td align="center"><img src="structures/H2MoOH2.png" width="200"><br><b>H₂MoOH₂</b></td>
  <td align="center"><img src="structures/MoOOH3.png" width="200"><br><b>MoO(OH)₃</b></td>
</tr>
</table>

<sub>Teal = Mo, red = O, white = H, copper = Cu(111) substrate.</sub>

## The free-energy expression

Two ingredients set each species' free energy at a grid point `(U, pH)`.

**1. Proton-coupled electron transfer (computational hydrogen electrode):**

```
G(U) = G_system - G_Mo - m*G_H2O + (m - n/2)*G_H2 + (n - 2m)*q*U
```

with `m` = water molecules exchanged, `n` = H atoms, `q` = 1 electron per proton.
Potential `U` is referenced to RHE.

**2. The local interfacial electric field (the field correction).** A Helmholtz
double-layer model turns the electrode potential and pH into the field felt by the
adsorbate, which then tilts its energy via the dipole `mu` and polarizability
`alpha` fitted in step 1:

```
E(U, pH) = C_H * (U - kT*ln10*pH/e - U_PZC) / (eps * eps0)        [V/m]
G_field  = G(U) + mu*E - (alpha/2)*E**2
```

| Symbol      | Meaning                                            | Value          |
|-------------|----------------------------------------------------|----------------|
| `C_H`       | Helmholtz capacitance                              | 0.25 F/m²      |
| `U_PZC`     | potential of zero charge                           | −0.70 V        |
| `eps`,`eps0`| inner-layer / vacuum permittivity                  | 2, 8.854e−12 F/m |
| `mu`        | surface dipole of the species (rel. to Mo)         | fitted (step 1)|
| `alpha`     | polarizability of the species (rel. to Mo)         | fitted (step 1)|

The species with the lowest `G_field` wins each cell.

## How the field response is obtained

`mu` and `alpha` are **not** assumed — they are fitted, in step 1, to a DFT (VASP)
field scan: under an applied field `E`, each species' energy (relative to the bare
Mo reference) follows `dG_rel(E) = mu*E − (alpha/2)*E²`. This keeps the workflow
self-consistent: the polarization parameters come from the same level of theory as
the energies.

## Interpreting the diagram

The grid procedure is exact, but the diagram inherits every approximation in the
inputs — the functional, the choice of energy (electronic vs. Gibbs), and the
double-layer parameters (`C_H`, `U_PZC`, `eps`). A diagram should always be read
together with the level of theory and the double-layer model behind it.

## Repository layout

```
cu-mo-field-pourbaix-tiles/
├── data/
│   ├── total_energy_vs_field.csv       # DFT (VASP) total energies vs applied field
│   ├── formation_terms_cu111.csv       # G_system, G_Mo, G_H2O, G_H2, m, n per species
│   ├── field_response_fit.csv          # mu, alpha per species (from step 1)
│   └── README.md                       # data dictionary + provenance
├── scripts/
│   ├── step1_fit_field_response.py     # field scan -> mu, alpha per species
│   ├── step2_rank_stability_grid.py    # mu/alpha + formation terms -> winner per cell
│   └── step3_plot_pourbaix_diagram.py  # stability grid -> field-corrected Pourbaix diagram
├── structures/                         # DFT structure images of the four adsorbates
├── results/
│   └── pourbaix_diagram_cu111.png      # (stability_grid_*.csv is regenerated, not committed)
├── requirements.txt
├── LICENSE
└── README.md
```

## Usage

```bash
pip install -r requirements.txt

python scripts/step1_fit_field_response.py     # data/field_response_fit.csv
python scripts/step2_rank_stability_grid.py     # results/stability_grid_cu111.csv
python scripts/step3_plot_pourbaix_diagram.py   # results/pourbaix_diagram_cu111.png
```

Each script has a configuration block at the top — potential/pH window, grid
resolution, double-layer constants, colours — so the diagram can be re-tuned
without touching the logic.

## Requirements

Python ≥ 3.9 with `numpy`, `pandas`, `scipy`, `matplotlib` (see `requirements.txt`).
The figure uses Times New Roman where available and falls back to a generic serif.

## License

Released under the MIT License (see `LICENSE`).
