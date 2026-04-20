# src/cg — Coarse-grained layer (Martini 3)

**Campaign 1, Layer 1.6.**

## Scope
Lipid bilayers and large protein complexes at Martini 3 resolution
— ~10× faster than atomistic, required for membrane-embedded drug
targets and whole-organelle physics.

## Upstream tools (build-vs-buy: BUY)
- **Martini 3** force field (Apache-2).
- **martinize2 / vermouth** (Apache-2) — protein → CG topology.
- **OpenMM** with Martini params — same engine as atomistic.

## What we write
- Topology builder glue (POPC / POPE / cholesterol mixes).
- Protein elastic-network generator wrapper.
- Area-per-lipid / membrane-thickness analysis scripts.

## Exit test
10 µs of 200 nm² POPC bilayer in < 48 h on 1 GPU;
area-per-lipid within 2 % of published value.

## Viewer
Martini bilayer with head-groups colour-coded; live area-per-lipid plot.
