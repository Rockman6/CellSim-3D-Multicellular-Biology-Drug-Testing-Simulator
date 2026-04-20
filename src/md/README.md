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

## Exit tests
- **Layer 1.2 MVP (shipped):** ≥ 8/10 canonical drugs from
  `benchmarks/chembl/smoke_10.smi` run 10 ps vacuum Langevin MD
  with no NaNs, final temperature within 50 K of 300 K setpoint,
  and RMSD < 10 Å vs frame 0. Run:
  ```bash
  conda activate cellsim
  python tests/md/test_ligand_vacuum.py
  ```
- **Full Layer 1.2 gate (deferred):** 100 ns ubiquitin MD in
  explicit TIP3P water with PBC, backbone RMSD stabilising
  < 3 Å vs minimised start over the final 50 ns. Needs PDBFixer
  + ff14SB protein loader (next PR) and GPU (reference: 1 day
  on a single H100).

## Modules
- `simulate.py` — `simulate_ligand(smiles_or_param, *, n_steps,
  dt_fs, temperature_K, platform, …) → TrajectoryResult`. Rebuilds
  the OpenMM System from OpenFF, picks the best platform (CUDA >
  OpenCL > CPU), minimises, sets Maxwell-Boltzmann velocities at
  the setpoint temperature with a seeded RNG, Langevin NVT for the
  requested steps, reports per-frame telemetry.
- `viewer.py` — matplotlib animation (GIF / MP4) with 3D
  ball-and-stick trajectory + T / E / RMSD telemetry panel.

## Viewer quickstart
```bash
conda activate cellsim
python -m src.md.viewer aspirin \
    --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" \
    --save aspirin_md.gif
```
Produces a side-by-side view: animated ligand on the left,
temperature / potential-energy / RMSD time-series on the right with
a marker that tracks the current frame.
