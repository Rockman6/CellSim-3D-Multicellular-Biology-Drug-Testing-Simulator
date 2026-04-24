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
    """Coupled single-λ schedule (legacy). Deprecated in favour
    of `_split_lambda_schedule` — see Milestone A post-mortem in
    BENCHMARKS.md for the water-penetration artifact that motivates
    the switch."""
    if n_windows < 2:
        n_windows = 2
    return [i / (n_windows - 1) for i in range(n_windows)]


def _split_lambda_schedule(
    n_windows: int,
) -> list[tuple[float, float]]:
    """Canonical absolute-hydration FEP schedule: decouple
    electrostatics first, then sterics. Returns a list of
    (λ_electrostatics, λ_sterics) tuples indexed 0..n_windows-1
    where index 0 is the decoupled state (both 0) and index
    n_windows-1 is the coupled state (both 1).

    Going from coupled (index K-1) to decoupled (index 0):
      stage A: λ_sterics = 1, λ_electrostatics: 1 → 0
      stage B: λ_electrostatics = 0, λ_sterics: 1 → 0

    Why split: the openmmtools default of coupling
    `lambda_electrostatics = lambda_sterics = λ` creates a
    partially-attractive well at intermediate λ that waters
    penetrate. This inflates F(decoupled) on polar solutes by
    4-6 kcal/mol — the root cause of the Milestone A pilot-1
    ethanol overshoot. Decoupling electrostatics first while
    keeping full LJ repulsion excludes water from the ligand
    volume throughout, eliminating the penetration artifact.

    Layout with n_windows=11 (the Milestone A default): 6 elec
    windows including both endpoints, 5 sterics windows sharing
    the (elec=0, sterics=1) midpoint so no duplicate state.
    """
    if n_windows < 3:
        n_windows = 3
    # Split ~half-and-half. With n=11: 6 elec + 5 sterics (the
    # (elec=0, sterics=1) midpoint belongs to the elec stage).
    n_elec = (n_windows + 1) // 2     # 6 for n=11
    n_st = n_windows - n_elec + 1     # 6 for n=11 — but the (0,1)
                                      # point is shared, so we only
                                      # emit n_st - 1 = 5 new tuples.

    schedule: list[tuple[float, float]] = []
    # Sterics stage (index 0 = decoupled both): sterics = 0 → 1
    # while elec = 0. Yields (0,0), (0,1/(n_st-1)), ..., (0, 1).
    for i in range(n_st):
        lam_s = i / (n_st - 1)
        schedule.append((0.0, lam_s))
    # Elec stage starts from the shared (elec=0, sterics=1) which
    # is already at schedule[-1]. Continue from elec=1/(n_elec-1)
    # through elec=1, all at sterics=1.
    for i in range(1, n_elec):
        lam_e = i / (n_elec - 1)
        schedule.append((lam_e, 1.0))
    # Sanity: endpoints
    assert schedule[0] == (0.0, 0.0), schedule[0]
    assert schedule[-1] == (1.0, 1.0), schedule[-1]
    assert len(schedule) == n_windows, (len(schedule), n_windows)
    return schedule


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
        Platform,
    )
    from openmmtools import mcmc, states, cache
    from openmmtools.alchemy import AlchemicalState

    # Prefer the fastest available platform: Metal (Apple silicon
    # OpenMM 8.2+) → OpenCL (works on Apple silicon via Metal shim,
    # Intel, AMD, NVIDIA) → CUDA → CPU. On an M5 Max this gives a
    # ~5-10× speedup over CPU. openmmtools' context cache uses the
    # first platform available in its list unless we pin one.
    _preferred = ("Metal", "OpenCL", "CUDA", "CPU")
    _chosen_platform = None
    for _name in _preferred:
        try:
            _p = Platform.getPlatformByName(_name)
            cache.global_context_cache.platform = _p
            _chosen_platform = _name
            break
        except Exception:
            continue
    logger.info(
        "FEP sampling platform: %s", _chosen_platform or "default")

    t0 = time.time()
    # Canonical split schedule: decouple electrostatics first, then
    # sterics. See Milestone A post-mortem in BENCHMARKS.md for why
    # the legacy coupled-λ schedule was a water-penetration bug.
    schedule = _split_lambda_schedule(n_windows)
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
    for lam_e, lam_s in schedule:
        alch = AlchemicalState(
            lambda_electrostatics=lam_e, lambda_sterics=lam_s)
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
    # starting any MD. Direct setPositions + LocalEnergyMinimizer,
    # bypassing CompoundThermodynamicState.apply_to_context (which
    # expects a thermostat present in the System and raises a
    # silent AttributeError if the OpenMM system from openff
    # Interchange doesn't carry one). LocalEnergyMinimizer doesn't
    # need a thermostat; it only needs the bare System + positions.
    # After minimise, call applyConstraints(1e-6) to enforce H-bond
    # and water-SETTLE constraints that minimisation alone can
    # leave slightly off by ~0.1 Å — GHMC rejects 100% of proposals
    # when given positions with those violations.
    def _apply_alchemical_params(ctx, comp_state):
        """Set lambda_* globals on the context without touching a
        thermostat. Extracted from CompoundThermodynamicState."""
        for cs in comp_state._composable_states:
            cs.apply_to_context(ctx)
    try:
        _min_int = Verlet(1.0 * ommunit.femtosecond)
        _min_ctx = _Context(alch_system, _min_int)
        _apply_alchemical_params(_min_ctx, compound_states[-1])
        _min_ctx.setPositions(sampler_state.positions)
        if sampler_state.box_vectors is not None:
            _min_ctx.setPeriodicBoxVectors(
                *sampler_state.box_vectors)
        _Min.minimize(_min_ctx, tolerance=1.0, maxIterations=2000)
        _min_ctx.applyConstraints(1e-6)
        sampler_state.update_from_context(_min_ctx)
        del _min_ctx, _min_int
    except Exception as e:
        logger.warning("initial minimisation failed: %s", e)

    # Iterate the schedule in REVERSE: start fully coupled
    # (λ=1), which is the physical initial state from packmol
    # (ligand embedded in a water box, all interactions on), and
    # decouple as λ → 0. Running forward would start at λ=0 with
    # the ligand already phantom, waters collapse into its space
    # before we can equilibrate at any coupled state, and the
    # integrator NaNs on the first step because waters are sitting
    # on top of each other where the ligand used to be. This is
    # the canonical FEP order used by perses/yank.
    # Persistent eval-context: used ONLY for cross-λ energy
    # evaluation via direct AlchemicalState control. openmmtools'
    # CompoundThermodynamicState.reduced_potential fails on
    # thermostat-less systems (openff Interchange doesn't add
    # one); direct control is the canonical fallback used by
    # perses when the system lacks a thermostat. We build this
    # once and reuse across all sample collections.
    _eval_int = Verlet(1.0 * ommunit.femtosecond)
    _eval_ctx = _Context(alch_system, _eval_int)
    _eval_state = AlchemicalState.from_system(alch_system)
    beta = 1.0 / (0.0019872041 * temperature_K)  # 1/(kT) kcal/mol

    schedule_order = list(reversed(list(enumerate(schedule))))
    for k, (lam_e, lam_s) in schedule_order:
        comp_state = compound_states[k]

        # Professor's checklist item 3: minimise before each λ
        # window + applyConstraints. Bypass CompoundTS.apply_to_
        # context (no thermostat in the system).
        try:
            _mi = Verlet(1.0 * ommunit.femtosecond)
            _mctx = _Context(alch_system, _mi)
            _apply_alchemical_params(_mctx, comp_state)
            _mctx.setPositions(sampler_state.positions)
            if sampler_state.box_vectors is not None:
                _mctx.setPeriodicBoxVectors(
                    *sampler_state.box_vectors)
            _Min.minimize(_mctx, tolerance=1.0,
                          maxIterations=2000)
            _mctx.applyConstraints(1e-6)
            sampler_state.update_from_context(_mctx)
            del _mctx, _mi
        except Exception as e:
            logger.warning(
                "per-λ minimise FAILED at λ_e=%.2f λ_s=%.2f: %s",
                lam_e, lam_s, e)

        # Reset the equilibration move's acceptance counter
        # before each window so the per-window statistic is clean.
        equil_move.reset_statistics()
        equil_move.apply(comp_state, sampler_state)

        # Reset production-move counter for the window.
        prod_move.reset_statistics()
        for s in range(n_samples_per):
            col = k * n_samples_per + s
            prod_move.apply(comp_state, sampler_state)
            # Eval reduced potential at every target λ via direct
            # control: setPositions → set λ via AlchemicalState →
            # read potential energy → u = β·U.
            _eval_ctx.setPositions(sampler_state.positions)
            if sampler_state.box_vectors is not None:
                _eval_ctx.setPeriodicBoxVectors(
                    *sampler_state.box_vectors)
            for n, (lam_en, lam_sn) in enumerate(schedule):
                _eval_state.lambda_electrostatics = lam_en
                _eval_state.lambda_sterics = lam_sn
                _eval_state.apply_to_context(_eval_ctx)
                U = _eval_ctx.getState(
                    getEnergy=True).getPotentialEnergy(
                    ).value_in_unit(ommunit.kilocalorie_per_mole)
                u_kn[n, col] = beta * U
        # Record acceptance for this window.
        try:
            result.ghmc_acceptance.append(
                float(prod_move.fraction_accepted))
        except Exception:
            result.ghmc_acceptance.append(0.0)
    del _eval_ctx, _eval_int

    try:
        dG_kT, ddG_kT = estimate_free_energy(u_kn, N_k)
    except Exception as e:
        result.reason = _biologist_reason_for_mbar_error(
            e, n_windows=n_windows,
            n_production_steps=n_production_steps)
        result.n_samples_per_window = n_samples_per
        result.wall_seconds = time.time() - t0
        return result
    kT = 0.0019872041 * temperature_K
    result.ok = True
    result.n_samples_per_window = n_samples_per
    result.dG_kcalmol = dG_kT * kT
    result.dG_uncertainty_kcalmol = ddG_kT * kT
    result.wall_seconds = time.time() - t0
    return result


def _biologist_reason_for_mbar_error(
    exc: BaseException,
    *,
    n_windows: int,
    n_production_steps: int,
) -> str:
    """Translate pymbar/MBAR internal exceptions into biologist-
    actionable guidance. The raw pymbar text says things like
    "Column sum W_nk = 0" which is not meaningful to anyone who
    hasn't read the MBAR paper. The translation names the physical
    cause (adjacent λ-windows don't overlap, sampling too short)
    and suggests exact parameter bumps.
    """
    msg = str(exc).lower()
    # Suggested bumps: 2× windows, 3× prod steps, capped so we
    # don't print absurdly large numbers that scare off users on
    # laptops. These are the Milestone-A-tier numbers that
    # reliably converge streptavidin-class binders on GPU.
    suggest_windows = min(max(n_windows * 2, 21), 31)
    suggest_prod = min(max(n_production_steps * 3, 25000), 100000)
    hint = (
        f" → try --n-windows {suggest_windows} "
        f"--n-production-steps {suggest_prod} "
        "(Milestone-A-tier sampling; GPU recommended)")

    # Pattern-match the distinctive pymbar failure strings.
    if ("column sum" in msg) or ("w_nk" in msg) or (
            "overlap" in msg):
        return (
            "adjacent λ-windows don't overlap — sampling was "
            "insufficient for MBAR to reweight between states. "
            "This is common for tight binders (|ΔG| > 10 kcal/mol) "
            "at smoke-tier parameters." + hint)
    if ("singular" in msg) or ("convergence" in msg):
        return (
            "MBAR self-consistent iteration did not converge — "
            "typically caused by states with near-identical "
            "reduced potentials (redundant λ schedule) or NaN/Inf "
            "samples." + hint)
    if ("bounds" in msg) or ("negative uncertainty" in msg):
        return (
            "MBAR returned non-physical uncertainty — too few "
            "independent samples per window." + hint)
    if ("nan" in msg) or ("inf" in msg):
        return (
            "reduced-potential matrix contains NaN/Inf — the "
            "integrator blew up on at least one window. Try a "
            "smaller timestep (--timestep-fs 0.5) or verify the "
            "input structure has no steric clashes." + hint)
    # Fallback: preserve the original text but frame it clearly.
    return (
        f"MBAR failed with: {str(exc)[:200]}." + hint)


def estimate_free_energy(u_kn, N_k) -> tuple[float, float]:
    """Run MBAR and return (ΔG_λ=0→λ=1, uncertainty) in kT.

    Raises on convergence / overlap failure. Callers should wrap
    with _biologist_reason_for_mbar_error to translate."""
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
