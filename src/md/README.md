# src/md — Classical MD engine

**Campaign 1, Layer 1.2.**

## Scope
Parameterised system → equilibrated trajectory. CUDA + Metal + CPU
platforms via OpenMM.

## Upstream tools (build-vs-buy: BUY)
- **OpenMM 8.x** (MIT) — primary MD backend.
- **AMBER ff14SB / ff19SB** — protein force field (bundled with OpenMM).
- **CHARMM36m** — alternative protein FF.
- **TIP3P** water, PME electrostatics.
- **GROMACS** (LGPL) — wrapped for independent-engine blind cross-check
  only, not primary.
- **MACE / MACE-OFF23** (MIT) — ML-potential via `openmm-ml` plugin.
  Deep-ensemble of 3–5 replicas for force-level UQ (Layer 1.5).

## What we write
- OpenMM driver that consumes `src/chem/` output.
- Trajectory reporter → shared-memory ring buffer for viewer.
- Ensemble wrapper for ML-potential UQ.

## Exit test
100 ns ubiquitin MD, RMSD < 3 Å; reproduces published benchmark
numbers on a single H100 in < 24 h.

## Viewer
Real-time protein ribbon + solvent skin; RMSD gauge;
temperature / pressure readout.
