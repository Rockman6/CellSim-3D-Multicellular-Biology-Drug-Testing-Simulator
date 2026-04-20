"""Martini 3 lipid-bilayer builder + OpenMM CG-MD driver.

**Scaffold only.** Planned API documented below.

Pipeline (to implement):

    build_martini3_bilayer(composition, box_nm, temperature_K)
        -> CGSystem(topology, positions, system)

    composition example:
        {"POPC": 0.7, "POPE": 0.2, "CHOL": 0.1}  # mol fractions

    Martini 3 force-field parameters are distributed with `vermouth`.
    Lipid initial placement via a geometric insane-style algorithm
    (regular lattice + random rotation + overlap rejection) — no
    external tool, pure Python. Solvate with Martini W beads.

    run_cg_md(cg_system, *, n_steps, dt_fs=20, temperature_K=310,
              pressure_bar=1.0, platform=None) -> CGTrajectoryResult

    Langevin NVT for equilibration → Martini-style pressure coupling
    for production. Stride-saved area-per-lipid, bilayer thickness,
    lipid lateral diffusion from the trajectory.

Exit gate (Layer 1.5):
  - 10 µs POPC bilayer, 200 nm² patch
  - area-per-lipid within 2 % of the published value (64.3 Å²/lipid
    at 310 K; Marrink 2007; ~65-67 Å² in Martini 3 per Souza 2021)
  - bilayer thickness matches published X-ray (~ 38 Å for POPC)

Non-AI: pure classical FF (Martini 3 parameters from literature),
OpenMM integration, published analysis pipelines. No learned model.

This file will be populated in a future iteration; for now it is a
placeholder so the import path exists and the scaffold test suite
can assert the layer is wired.
"""

from __future__ import annotations


def build_martini3_bilayer(*args, **kwargs):  # pragma: no cover
    raise NotImplementedError(
        "src.cg.bilayer.build_martini3_bilayer is planned but not "
        "yet implemented. See src/cg/README.md for scope.")


def run_cg_md(*args, **kwargs):  # pragma: no cover
    raise NotImplementedError(
        "src.cg.bilayer.run_cg_md is planned but not yet "
        "implemented. See src/cg/README.md for scope.")
