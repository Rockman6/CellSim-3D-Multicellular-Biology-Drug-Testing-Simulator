"""src/cell — Campaign-2 seed: cell-level readouts from Campaign-1 priors.

This is the FIRST consumer of the `src/bridge` rate-law priors — the
point where a molecular free energy becomes a cell-level quantity. It is
deliberately minimal (a single-compartment, single-site occupancy model)
so the end-to-end path exists and is honest before it is elaborated:

    binding ΔG  →  K_d prior (+ accuracy, + trust)  →  target occupancy
    (src/fep, src/dock)   (src/bridge)                 (here)

Everything carried through: the uncertainty on ΔG becomes a CI on
occupancy, and the target-class `trust` verdict rides along so a readout
computed from an untrustworthy K_d is never mistaken for a decision-grade
number.

Non-AI: closed-form Hill occupancy. No learned surrogate.

NOT yet here (later Campaign-2 work; see docs/campaign1_closeout.md):
  - multiple compartments (extracellular / cytoplasm) + membrane
    permeability, so occupancy is computed at the LOCAL concentration;
  - competitive / multi-drug binding at a shared site;
  - dynamics (ODE integration over time).
"""
from __future__ import annotations

from src.cell.occupancy import (
    OccupancyResult,
    occupancy_from_prior,
    hill_occupancy,
)

__all__ = ["OccupancyResult", "occupancy_from_prior", "hill_occupancy"]
