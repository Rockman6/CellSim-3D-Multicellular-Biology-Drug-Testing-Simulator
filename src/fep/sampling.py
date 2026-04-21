"""Alchemical MD sampling + MBAR for Layer 1.3 FEP.

Uses the canonical openmmtools MCMC pattern:

  CompoundThermodynamicState per λ (ThermodynamicState ∘
  AlchemicalState) → LangevinDynamicsMove.apply(state,
  sampler_state) runs MD → SamplerState.reduced_potential across
  compound states populates the u_kn matrix → pymbar.MBAR.

This replaces the hand-rolled (one Context, swap λ in-place,
save/restore velocities) pattern that was intermittent-NaN
on polar / aromatic solvated systems at smoke parameters —
the MCMC-move abstraction reassigns velocities from the
Maxwell-Boltzmann distribution at the start of each move,
decoupling per-λ dynamics from the cross-λ energy-evaluation
loop.

No replica exchange, no REST2 — this is the non-AI baseline
every openmmtools-based FEP pipeline shares.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class AlchemicalSamplingResult:
    ok: bool
    reason: str = ""
    n_windows: int = 0
    n_samples_per_window: int = 0
    dG_kcalmol: Optional[float] = None
    dG_uncertainty_kcalmol: Optional[float] = None
    wall_seconds: Optional[float] = None
    lambda_schedule: list = field(default_factory=list)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] FEP sampling  {self.reason}"
        return (
            f"[OK]   FEP sampling  "
            f"ΔG = {self.dG_kcalmol:+.2f} ± "
            f"{self.dG_uncertainty_kcalmol:.2f} kcal/mol  "
            f"({self.n_windows} windows × "
            f"{self.n_samples_per_window} samples,  "
            f"wall {self.wall_seconds:.1f} s)")


def _default_lambda_schedule(n_windows: int) -> list[float]:
    if n_windows < 2:
        n_windows = 2
    return [i / (n_windows - 1) for i in range(n_windows)]


def sample_alchemical_windows(
    alch_system,
    topology,            # kept for API compat; unused
    positions,
    *,
    n_windows: int = 11,
    n_equilibration_steps: int = 500,
    n_production_steps: int = 2000,
    timestep_fs: float = 2.0,
    temperature_K: float = 298.15,
    friction_ps: float = 1.0,
    seed: int = 1,
    sample_stride: int = 100,
) -> AlchemicalSamplingResult:
    """Run Langevin MD at each λ via openmmtools MCMC moves;
    return ΔG between λ=0 and λ=1 estimated via MBAR."""
    import numpy as np
    from openmm import unit as ommunit
    from openmmtools import mcmc, states
    from openmmtools.alchemy import AlchemicalState

    t0 = time.time()
    schedule = _default_lambda_schedule(n_windows)
    result = AlchemicalSamplingResult(
        ok=False, n_windows=n_windows,
        lambda_schedule=schedule)

    # Build the base ThermodynamicState once; then wrap it with
    # an AlchemicalState per λ via CompoundThermodynamicState.
    try:
        thermo = states.ThermodynamicState(
            system=alch_system,
            temperature=temperature_K * ommunit.kelvin)
        alch_state_template = AlchemicalState.from_system(
            alch_system)
    except Exception as e:
        result.reason = (
            f"ThermodynamicState / AlchemicalState init failed: "
            f"{str(e)[:200]}")
        return result

    compound_states = []
    for lam in schedule:
        alch = AlchemicalState(
            lambda_electrostatics=lam, lambda_sterics=lam)
        compound_states.append(
            states.CompoundThermodynamicState(
                thermodynamic_state=states.ThermodynamicState(
                    system=alch_system,
                    temperature=temperature_K * ommunit.kelvin),
                composable_states=[alch]))

    # Initial sampler state.
    sampler_state = states.SamplerState(
        positions=positions,
        box_vectors=(
            alch_system.getDefaultPeriodicBoxVectors()
            if alch_system.usesPeriodicBoundaryConditions()
            else None))

    n_samples_per = max(1, n_production_steps // sample_stride)
    n_states = n_windows
    u_kn = np.zeros((n_states, n_states * n_samples_per))
    N_k = np.array([n_samples_per] * n_states)

    # LangevinDynamicsMove.apply creates a fresh context each
    # call (that's the whole point — no context reuse, no stale
    # integrator state, no NaN from swap-induced kicks).
    equil_move = mcmc.LangevinDynamicsMove(
        timestep=timestep_fs * ommunit.femtosecond,
        collision_rate=friction_ps / ommunit.picosecond,
        n_steps=n_equilibration_steps,
        reassign_velocities=True)
    prod_move = mcmc.LangevinDynamicsMove(
        timestep=timestep_fs * ommunit.femtosecond,
        collision_rate=friction_ps / ommunit.picosecond,
        n_steps=sample_stride,
        reassign_velocities=False)

    col = 0
    for k, lam in enumerate(schedule):
        comp_state = compound_states[k]
        # Equilibrate at this λ. reassign_velocities=True at
        # the start of each window decouples windows.
        equil_move.apply(comp_state, sampler_state)
        # Production + cross-λ energy evaluation.
        for _ in range(n_samples_per):
            prod_move.apply(comp_state, sampler_state)
            # Compute reduced potential at every target λ.
            # SamplerState.reduced_potential handles context
            # creation + energy eval; no state scrambling.
            for n in range(n_states):
                u = compound_states[n].reduced_potential(
                    sampler_state)
                u_kn[n, col] = u
            col += 1

    dG_kT, ddG_kT = estimate_free_energy(u_kn, N_k)
    kT = 0.0019872041 * temperature_K
    result.ok = True
    result.n_samples_per_window = n_samples_per
    result.dG_kcalmol = dG_kT * kT
    result.dG_uncertainty_kcalmol = ddG_kT * kT
    result.wall_seconds = time.time() - t0
    return result


def estimate_free_energy(u_kn, N_k) -> tuple[float, float]:
    """Run MBAR and return (ΔG_λ=0→λ=1, uncertainty) in kT."""
    from pymbar import MBAR

    mbar = MBAR(u_kn, N_k)
    try:
        r = mbar.compute_free_energy_differences()
        Delta_f = r["Delta_f"]
        dDelta_f = r["dDelta_f"]
    except AttributeError:
        Delta_f, dDelta_f = mbar.getFreeEnergyDifferences()
    K = len(N_k)
    return float(Delta_f[K - 1, 0]), float(dDelta_f[K - 1, 0])
