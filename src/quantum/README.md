# src/quantum — Semi-empirical + DFT layer

**Campaign 1, Layer 1.4.**

## Scope
Reactive-fragment energies, HOMO/LUMO, electrostatic potential,
CYP-metabolite site prediction.

## Upstream tools (build-vs-buy: BUY)
- **xtb (GFN2-xTB)** (LGPL) — semi-empirical quantum, primary.
- **PySCF** (Apache-2) — DFT for hard reactive fragments.
- **ORCA** (free academic) — optional fallback.
- **FAME 3 / XenoSite / SyGMa** — metabolism-site predictors;
  we re-rank with xTB.

## What we write
- xtb subprocess wrapper + energy/force parser.
- PySCF Python driver for DFT single-point calculations.
- CYP-substrate rank fusion: FAME 3 prior × xTB energy margin.

## Exit test
Reactive-metabolite prediction matches literature on ≥ 15/20
marketed CYP3A4 substrates.

## Viewer
Electrostatic-potential isosurface on the ligand;
HOMO / LUMO density clouds.
