"""Alchemical MD sampling + MBAR for Layer 1.3 FEP.

Minimal non-AI FEP engine: run short Langevin MD at each λ state,
collect reduced potentials on a common grid, and BAR/MBAR-combine
into a free-energy difference.

No replica exchange, no REST2, no neural surrogate — this is the
baseline that every non-AI FEP pipeline shares. Accuracy is set
by (n_windows, n_steps_per_window); neither parameter affects the
shape of the code.

Entry points:
  sample_alchemical_windows(...)  → u_kn matrix + N_k
  estimate_free_energy(u_kn, N_k) → (ΔG_kcalmol, dΔG_kcalmol)
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
    """Linear schedule in λ ∈ [0, 1]. Turns both electrostatics
    and sterics off simultaneously. Real production FEP decouples
    electrostatics first (λ_elec 1→0) then sterics (λ_ster 1→0)
    to avoid the Coulomb singularity; we collapse to one λ for
    simplicity since the ligand is small and softcore handles the
    sterics endpoint — good enough for hydration of neutral
    drug-like molecules."""
    if n_windows < 2:
        n_windows = 2
    return [i / (n_windows - 1) for i in range(n_windows)]


def sample_alchemical_windows(
    alch_system,
    topology,
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
    """Run Langevin MD at each λ; return the u_kn matrix for MBAR.

    The returned result carries ΔG already computed from MBAR.
    For a FreeSolv-grade absolute hydration number set
    n_production_steps ≥ 50_000 (=100 ps at 2 fs) and
    n_windows = 12. Defaults here are a fast-smoke tuning
    (~10 s/compound) that verifies the pipeline runs end-to-end.
    """
    import numpy as np
    from openmm import (
        LangevinMiddleIntegrator,
        LocalEnergyMinimizer,
        unit as ommunit,
        Context,
        Platform,
    )
    from openmmtools.alchemy import AlchemicalState

    t0 = time.time()
    schedule = _default_lambda_schedule(n_windows)
    result = AlchemicalSamplingResult(
        ok=False, n_windows=n_windows,
        lambda_schedule=schedule)

    # Requires the alchemical region to have been built with
    # name='' (default) so the lambda_sterics / lambda_electrostatics
    # globals are unsuffixed. src/fep/__init__.py does this already.
    try:
        alch_state = AlchemicalState.from_system(alch_system)
    except Exception as e:
        result.reason = (
            f"AlchemicalState.from_system failed: "
            f"{str(e)[:200]}")
        return result

    # Pick the fastest platform we can get without caring which —
    # Reference is always present, CUDA/OpenCL/Metal speeds are
    # a bonus. This is a small system (dozens to hundreds of
    # atoms), so we don't need a GPU for the pipeline smoke.
    try:
        platform = Platform.getPlatformByName("CPU")
    except Exception:
        platform = None

    T = temperature_K * ommunit.kelvin
    integrator = LangevinMiddleIntegrator(
        T, friction_ps / ommunit.picosecond,
        timestep_fs * ommunit.femtosecond)
    integrator.setRandomNumberSeed(int(seed))

    ctx = Context(alch_system, integrator, platform) \
        if platform else Context(alch_system, integrator)
    ctx.setPositions(positions)

    kT = (0.0019872041 * temperature_K)  # kcal/mol (R*T at T_K)

    n_samples_per = max(1, n_production_steps // sample_stride)
    n_states = n_windows
    u_kn = np.zeros((n_states, n_states * n_samples_per))
    N_k = np.array([n_samples_per] * n_states)

    col = 0
    for k, lam in enumerate(schedule):
        # Update both electrostatics and sterics to same λ.
        alch_state.lambda_electrostatics = lam
        alch_state.lambda_sterics = lam
        alch_state.apply_to_context(ctx)

        # Minimise after each λ change: packmol placement +
        # alchemical coupling together produce close contacts
        # that send Langevin to NaN on the first integrator
        # step. Tight tolerance + many iterations is needed for
        # solvated systems with polar/aromatic solutes.
        try:
            LocalEnergyMinimizer.minimize(
                ctx, tolerance=1.0, maxIterations=2000)
        except Exception as e:
            logger.debug(
                "minimisation skipped at λ=%.2f: %s", lam, e)

        # Equilibrate.
        integrator.step(n_equilibration_steps)

        # Production: at each stride, evaluate the reduced
        # potential at EVERY λ state (to populate u_kn).
        for _ in range(n_samples_per):
            integrator.step(sample_stride)
            # Snapshot full state — positions AND velocities —
            # so the inner cross-λ energy-evaluation loop can
            # restore them cleanly. Without saving velocities
            # the Langevin integrator's internal velocity
            # cache goes stale after setPositions calls inside
            # the inner loop and the next integrator.step()
            # blows up with a NaN. This is the load-bearing
            # fix for solvated-system sampling.
            state = ctx.getState(
                getPositions=True, getVelocities=True)
            pos = state.getPositions(asNumpy=True)
            vel = state.getVelocities(asNumpy=True)
            # Evaluate this sample at each target state n.
            for n, lam_n in enumerate(schedule):
                alch_state.lambda_electrostatics = lam_n
                alch_state.lambda_sterics = lam_n
                alch_state.apply_to_context(ctx)
                ctx.setPositions(pos)
                U = ctx.getState(
                    getEnergy=True
                ).getPotentialEnergy().value_in_unit(
                    ommunit.kilocalorie_per_mole)
                u_kn[n, col] = U / kT
            # Restore λ_k, positions AND velocities for the
            # next production step.
            alch_state.lambda_electrostatics = lam
            alch_state.lambda_sterics = lam
            alch_state.apply_to_context(ctx)
            ctx.setPositions(pos)
            ctx.setVelocities(vel)
            col += 1

    # Free-energy estimate via MBAR.
    dG_kT, ddG_kT = estimate_free_energy(u_kn, N_k)
    result.ok = True
    result.n_samples_per_window = n_samples_per
    result.dG_kcalmol = dG_kT * kT
    result.dG_uncertainty_kcalmol = ddG_kT * kT
    result.wall_seconds = time.time() - t0
    return result


def estimate_free_energy(u_kn, N_k) -> tuple[float, float]:
    """Run MBAR on a (K × total_samples) reduced-potential matrix
    and return (ΔG_λ=0→λ=1, uncertainty) in kT units."""
    from pymbar import MBAR

    mbar = MBAR(u_kn, N_k)
    try:
        # pymbar 4 API: compute_free_energy_differences returns
        # dict with 'Delta_f' (matrix) and 'dDelta_f' keys.
        r = mbar.compute_free_energy_differences()
        Delta_f = r["Delta_f"]
        dDelta_f = r["dDelta_f"]
    except AttributeError:
        Delta_f, dDelta_f = mbar.getFreeEnergyDifferences()
    # ΔG from state 0 (λ=0 = decoupled) to state K-1 (λ=1 = fully
    # interacting). Sign convention: we report ΔG(decouple) =
    # G(λ=0) − G(λ=1) = Delta_f[K-1, 0] which is positive when
    # the interacting state is more stable.
    K = len(N_k)
    return float(Delta_f[K - 1, 0]), float(dDelta_f[K - 1, 0])
