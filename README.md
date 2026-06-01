# Tessellated Ab-Initio Pourbaix (TAIP) Analysis

Build **surface and cluster Pourbaix (potential–pH) diagrams from first-principles
energies** by tessellating the potential–pH plane into a fine grid of tiles and,
in each tile, simply keeping the species with the lowest free energy.

The point of the approach is what it *avoids*: there is no analytical construction
of the equilibria between species — no `pKa` of each protonation step, no
equilibrium potential of each redox step, no Nernst-line intersections to solve.
You only need a single closed-form free energy per species as a function of `U`
and `pH`; the phase boundaries then **emerge on their own** as the contours where
the winning species changes.

## The two case studies

| | Case 1 — [`case1-iridium-pourbaix-tiles`](case1-iridium-pourbaix-tiles) | Case 2 — [`case2-cu-mo-field-pourbaix-tiles`](case2-cu-mo-field-pourbaix-tiles) |
|---|---|---|
| System | IrOₓ nanoparticle clusters | MoOₘHₙ on Cu(111) |
| Competing states | discrete charge + protonation states | hydride/oxyhydride adsorbates |
| Potential reference | SHE | RHE |
| Extra physics | — | explicit interfacial-field correction (dipole μ, polarizability α, Helmholtz double layer) |
| Case study | Bhattacharyya, Poidevin & Auer, *J. Phys. Chem. C* 2021 | Kambale, Rivera Rocabado, Kanematsu & Ishimoto, *Preprints.org* 2026 |

Each case is a self-contained, runnable repository (its own data, three staged
scripts, README). Case 2 is the field-aware generalisation of Case 1: same
tile-and-rank engine, with each species additionally responding to the local
electric field.

<p align="center">
  <img src="case1-iridium-pourbaix-tiles/results/pourbaix_diagram_bhattacharyya.png" width="380" alt="Case 1 diagram">
  &nbsp;&nbsp;
  <img src="case2-cu-mo-field-pourbaix-tiles/results/pourbaix_diagram_cu111.png" width="380" alt="Case 2 diagram">
</p>

## How it works, in three steps

Both cases share the same pipeline:

1. **prepare** the species table (energies, charge/H counts, field response);
2. **rank** every species at the centre of every grid tile and record the winner;
3. **plot** the tiles coloured by their winning species.

## Grid resolution, granularity, and cost

The diagram is a **discretization**, and that has a visible consequence: a phase
boundary is only ever resolved to the width of one tile, so at coarse resolution
boundaries look like a **staircase**. Refining the grid smooths them — at a price.

The cost is essentially

```
work  ∝  (number of tiles) × (number of species)
number of tiles  =  (ΔU / δU) × (ΔpH / δpH)
```

so **halving the tile size in both axes quadruples the number of tiles** (and the
runtime). The boundaries get ~2× sharper; the compute grows ~4×. There is no
accuracy gain in the thermodynamics from a finer grid — only a smoother-looking
boundary — so the resolution is a visual/cosmetic choice, not a physical one.

In practice this stays cheap because the ranking is evaluated as a single
vectorized operation (next section). Representative single-core timings for the
ranking step:

| grid (pH × U) | tiles | species | time (1 core) |
|---------------|-------|---------|---------------|
| 350 × 350     | 122,500   | 77 | ~0.1 s  |
| 1400 × 200    | 280,000   | 4  | ~0.01 s |
| 1400 × 1400   | 1,960,000 | 8  | ~0.17 s |
| 2800 × 2800   | 7,840,000 | 8  | ~0.85 s |

So a user who wants a finer diagram simply sets a smaller `δU`/`δpH` (more, smaller
tiles) and pays a little more time — but even multi-million-tile grids finish in
well under a second on a single core.

## Performance: vectorization (and why caching/parallelism aren't needed)

Ranking is *embarrassingly parallel* — every tile is independent. Rather than loop
over tiles in Python, the scripts evaluate the free energy as one **NumPy broadcast
over a `(species, U, pH)` array** and take an `argmin` along the species axis. The
entire grid is computed in a handful of array operations, which is why the timings
above are sub-second.

This is a deliberate simplification of an earlier implementation that used
per-point memoization (`functools.lru_cache`) and a thread pool
(`concurrent.futures`). With full vectorization those add complexity without
benefit at the grid sizes used here: NumPy already runs the arithmetic in
optimized C over the whole array. The only situation that needs more is a grid so
fine (or a species set so large) that the `(species, U, pH)` array exceeds memory
— then evaluate it in **chunks** (split the grid into blocks and rank each block),
still vectorized per block, optionally across processes. For the diagrams here
that is never required.

## Repository layout

```
TAIP/
├── README.md                          # this page
├── case1-iridium-pourbaix-tiles/      # charge/protonation-state Pourbaix (SHE)
└── case2-cu-mo-field-pourbaix-tiles/  # field-corrected Pourbaix (RHE)
```

(If your case folders are named without the `case1-`/`case2-` prefix, adjust the
links above accordingly.)

## Citations

- Bhattacharyya, K.; Poidevin, C.; Auer, A. A. *Structure and Reactivity of IrOₓ
  Nanoparticles for the Oxygen Evolution Reaction in Electrocatalysis: An Electronic
  Structure Theory Study.* **J. Phys. Chem. C** 2021, *125*, 4379–4390.
  [doi:10.1021/acs.jpcc.0c10092](https://doi.org/10.1021/acs.jpcc.0c10092)
- Kambale, E. M.; Rivera Rocabado, D. S.; Kanematsu, Y.; Ishimoto, T.
  *Field-Dependent Redox Thermodynamics of MoOₘHₙ Species on Cu(111) and Ni(111)
  Surfaces under Alkaline Hydrogen Evolution Conditions.* **Preprints.org** 2026.
  [doi:10.20944/preprints202604.0944.v1](https://doi.org/10.20944/preprints202604.0944.v1)

## License

MIT (see each case directory's `LICENSE`).
