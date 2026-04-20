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

CellSim is non-AI (see `MISSION.md`). ML potentials (MACE / NequIP /
Allegro / OrbNet) are explicitly **not** used as the force path.
UQ comes from Layer 1.6 (parameter sweeps + Sobol + Monte-Carlo),
not from neural ensembles.

## What we write
- OpenMM driver that consumes `src/chem/` output.
- Trajectory reporter → shared-memory ring buffer for viewer.

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
- `viewer.py` — matplotlib animation (GIF / MP4) for a ligand
  trajectory: 3D ball-and-stick + T / E / RMSD panel.
- `protein.py` — `load_protein_pdb(pdb_path, *, ff_protein,
  ff_water, ph, padding_nm, …) → ProteinSystemResult` (PDBFixer +
  AMBER14-all + TIP3P-FB solvation + energy minimisation) and
  `short_protein_md(result, n_steps, …) → ProteinTrajectoryResult`
  (Langevin NVT spot-check that also records per-frame Cα
  positions + Cα RMSD + residue metadata for the viewer).
- `protein_viewer.py` — matplotlib animation of a protein
  trajectory: 3D Cα trace coloured along sequence + T / E /
  Cα-RMSD panel. Run:
  ```bash
  conda activate cellsim
  python -m src.md.protein_viewer benchmarks/md/1ubq.pdb \
      --md-steps 2000 --save-every 50 --save 1ubq_md.gif
  ```

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
