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
    # GHMC acceptance per window — professor's acceptance-rate
    # requirement. < 0.70 means the timestep is too large and
    # the ΔG cannot be trusted.
    ghmc_acceptance: list = field(default_factory=list)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] FEP sampling  {self.reason}"
        accept_tag = ""
        if self.ghmc_acceptance:
            mean_acc = sum(self.ghmc_acceptance) / len(
                self.ghmc_acceptance)
            min_acc = min(self.ghmc_acceptance)
            accept_tag = (
                f"  GHMC accept mean={mean_acc:.0%} "
                f"min={min_acc:.0%}")
        return (
            f"[OK]   FEP sampling  "
            f"ΔG = {self.dG_kcalmol:+.2f} ± "
            f"{self.dG_uncertainty_kcalmol:.2f} kcal/mol  "
            f"({self.n_windows} windows × "
            f"{self.n_samples_per_window} samples,  "
            f"wall {self.wall_seconds:.1f} s){accept_tag}")


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
    timestep_fs: float = 1.0,       # prof-checklist: 1 fs for
                                    # debug; drop from 2.0
    temperature_K: float = 298.15,
    friction_ps: float = 20.0,      # prof-checklist: match GHMC
                                    # canonical friction
    seed: int = 1,
    sample_stride: int = 100,
) -> AlchemicalSamplingResult:
    """Run Langevin MD at each λ via openmmtools MCMC moves;
    return ΔG between λ=0 and λ=1 estimated via MBAR."""
    import numpy as np
    from openmm import (
        unit as ommunit,
        Context as _Context,
        VerletIntegrator as Verlet,
        LocalEnergyMinimizer as _Min,
    )
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

    # GHMCMove (Generalised Hybrid Monte Carlo) is Metropolised
    # Langevin — bad configurations get REJECTED instead of
    # blown up to NaN. The acceptance fraction is also a live
    # diagnostic: < 0.70 means the timestep is too large and
    # we cannot trust the ΔG. This is the professor's
    # "monitor the acceptance" requirement.
    equil_move = mcmc.GHMCMove(
        timestep=timestep_fs * ommunit.femtosecond,
        collision_rate=friction_ps / ommunit.picosecond,
        n_steps=n_equilibration_steps)
    prod_move = mcmc.GHMCMove(
        timestep=timestep_fs * ommunit.femtosecond,
        collision_rate=friction_ps / ommunit.picosecond,
        n_steps=sample_stride)

    # Minimise the initial state at λ=1 (fully coupled) before
    # starting any MD — packmol placement has close contacts
    # that send the first integrator step to NaN. Done under a
    # throwaway context so the MCMC move machinery gets a
    # relaxed sampler_state.
    try:
        _min_int = Verlet(1.0 * ommunit.femtosecond)
        _min_ctx = _Context(alch_system, _min_int)
        compound_states[-1].apply_to_context(_min_ctx)
        sampler_state.apply_to_context(_min_ctx)
        _Min.minimize(_min_ctx, tolerance=1.0, maxIterations=2000)
        sampler_state.update_from_context(_min_ctx)
        del _min_ctx, _min_int
    except Exception as e:
        logger.debug("initial minimisation skipped: %s", e)

    # Iterate the schedule in REVERSE: start fully coupled
    # (λ=1), which is the physical initial state from packmol
    # (ligand embedded in a water box, all interactions on), and
    # decouple as λ → 0. Running forward would start at λ=0 with
    # the ligand already phantom, waters collapse into its space
    # before we can equilibrate at any coupled state, and the
    # integrator NaNs on the first step because waters are sitting
    # on top of each other where the ligand used to be. This is
    # the canonical FEP order used by perses/yank.
    schedule_order = list(reversed(list(enumerate(schedule))))
    for k, lam in schedule_order:
        comp_state = compound_states[k]

        # Professor's checklist item 3: minimise before each λ
        # window. A "warm" config from the previous λ carries
        # overlaps the new Hamiltonian won't tolerate.
        try:
            _mi = Verlet(1.0 * ommunit.femtosecond)
            _mctx = _Context(alch_system, _mi)
            comp_state.apply_to_context(_mctx)
            sampler_state.apply_to_context(_mctx)
            _Min.minimize(_mctx, tolerance=1.0,
                          maxIterations=2000)
            sampler_state.update_from_context(_mctx)
            del _mctx, _mi
        except Exception as e:
            logger.debug(
                "per-λ minimise skipped at λ=%.2f: %s", lam, e)

        # Reset the equilibration move's acceptance counter
        # before each window so the per-window statistic is clean.
        equil_move.reset_statistics()
        equil_move.apply(comp_state, sampler_state)

        # Reset production-move counter for the window.
        prod_move.reset_statistics()
        for s in range(n_samples_per):
            col = k * n_samples_per + s
            prod_move.apply(comp_state, sampler_state)
            # Reduced potential at every target λ.
            for n in range(n_states):
                u = compound_states[n].reduced_potential(
                    sampler_state)
                u_kn[n, col] = u
        # Record acceptance for this window.
        try:
            result.ghmc_acceptance.append(
                float(prod_move.fraction_accepted))
        except Exception:
            result.ghmc_acceptance.append(0.0)

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
