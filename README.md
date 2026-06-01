# Tessellated Ab-Initio Pourbaix (TAIP) Analysis

Construct **surface and cluster Pourbaix (potential–pH) diagrams from
first-principles energies** by tiling the potential–pH plane into a fine grid and,
in each tile, keeping the species with the lowest free energy.

## The problem

A Pourbaix diagram shows which species is most stable as a function of electrode
potential `U` and `pH`. Built analytically, it requires — for every pair of
candidate species — writing the acid–base (`pKa`) and redox (Nernst) equilibrium
conditions and intersecting the resulting lines to bound each stability field. The
number of pairwise relations grows with the number of charge and protonation
states considered.

## The approach

Each candidate species has a single closed-form free energy as a function of `U`
and `pH` (the exact expression for each case study is given in that case's README).
The diagram then follows in three steps:

1. lay a fine grid over the potential–pH plane;
2. evaluate that free energy for **every** species at the center of each tile;
3. assign each tile to the species with the lowest free energy and color it accordingly.

A phase boundary is the locus where the lowest-energy species changes; it is read
off the grid rather than solved for. Adding a candidate species, charge state, or
protonation level means adding a row to the input table, not re-deriving the
diagram.

## The two case studies

Each case is a self-contained, runnable repository (its own data, three staged
scripts, README).

| | [`case1-iridium-pourbaix-tiles`](case1-iridium-pourbaix-tiles) | [`case2-cu-mo-field-pourbaix-tiles`](case2-cu-mo-field-pourbaix-tiles) |
|---|---|---|
| System | IrOₓ nanoparticle clusters | MoOₘHₙ on Cu(111) |
| Competing states | discrete charge + protonation states | hydride/oxyhydride adsorbates |
| Potential reference | SHE | RHE |
| Extra physics | — | interfacial-field correction (dipole μ, polarizability α, Helmholtz double layer) |
| Case study | Bhattacharyya, Poidevin & Auer, *J. Phys. Chem. C* 2021 | Kambale, Rivera Rocabado, Kanematsu & Ishimoto, *Preprints.org* 2026 |

<p align="center">
  <img src="case1-iridium-pourbaix-tiles/results/pourbaix_diagram_bhattacharyya.png" width="380" alt="Case 1 diagram">
  &nbsp;&nbsp;
  <img src="case2-cu-mo-field-pourbaix-tiles/results/pourbaix_diagram_cu111.png" width="380" alt="Case 2 diagram">
</p>

Case 2 is the field-aware counterpart of Case 1: the same tile-and-rank engine,
with each species additionally responding to the local electric field.

## Grid resolution and cost

The diagram is a discretization, so a phase boundary is resolved only to the width
of one tile; at coarse resolution boundaries appear as a staircase, and refining
the grid sharpens them. The ranking cost is

```
work  ∝  (number of tiles) × (number of species)
number of tiles  =  (ΔU / δU) × (ΔpH / δpH)
```

so halving the tile size in both axes multiplies the number of tiles (and the
runtime) by four. The thermodynamics in each tile is unchanged by the grid
spacing; only the boundary resolution changes. Representative single-core timings
for the ranking step:

| grid (pH × U) | tiles | species | time (1 core) |
|---------------|-------|---------|---------------|
| 350 × 350     | 122,500   | 77 | ~0.1 s  |
| 1400 × 200    | 280,000   | 4  | ~0.01 s |
| 1400 × 1400   | 1,960,000 | 8  | ~0.17 s |
| 2800 × 2800   | 7,840,000 | 8  | ~0.85 s |

A finer diagram is obtained by setting a smaller `δU`/`δpH` (more, smaller tiles).

## Performance: vectorization

The ranking is independent across tiles. Rather than loop over tiles in Python,
the scripts evaluate the free energy as one NumPy broadcast over a
`(species, U, pH)` array and take an `argmin` along the species axis, which is why
the timings above are sub-second. If the `(species, U, pH)` array exceeds memory,
the grid can be evaluated in chunks (split into blocks and rank each block), still
vectorized per block.

## Repository layout

```
TAIP/
├── README.md                          # this page
├── case1-iridium-pourbaix-tiles/      # charge/protonation-state Pourbaix (SHE)
└── case2-cu-mo-field-pourbaix-tiles/  # field-corrected Pourbaix (RHE)
```

(If your case folders are named without the `case1-`/`case2-` prefix, adjust the
links above accordingly.)

## Citation

If you use the TAIP method or this code, **please cite the method paper**:

> Kambale, E. M.; Rivera Rocabado, D. S.; Kanematsu, Y.; Ishimoto, T.
> *Field-Dependent Redox Thermodynamics of MoOₘHₙ Species on Cu(111) and Ni(111)
> Surfaces under Alkaline Hydrogen Evolution Conditions.* **Preprints.org** 2026.
> [doi:10.20944/preprints202604.0944.v1](https://doi.org/10.20944/preprints202604.0944.v1)

The Case 1 case-study energies are from:

> Bhattacharyya, K.; Poidevin, C.; Auer, A. A. *Structure and Reactivity of IrOₓ
> Nanoparticles for the Oxygen Evolution Reaction in Electrocatalysis: An Electronic
> Structure Theory Study.* **J. Phys. Chem. C** 2021, *125*, 4379–4390.
> [doi:10.1021/acs.jpcc.0c10092](https://doi.org/10.1021/acs.jpcc.0c10092)

## License

MIT (see each case directory's `LICENSE`).
