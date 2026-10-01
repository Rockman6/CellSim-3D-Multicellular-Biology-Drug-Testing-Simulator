"""Full cycle transient — cohort structure, synchronisation, dose timing.

`cycle.py` gives the asymptotic growth rate λ from a characteristic
equation, and `combination.py` builds schedules by adding λ·Δt across
blocks. That is a QUASI-STEADY-STATE picture: it assumes the phase
distribution has relaxed to its stable shape within every block. It cannot
represent the thing that makes real scheduling interesting.

Kill the S-phase cells and the survivors are, by construction, NOT
asynchronous — they are all in G1 and G2/M. That cohort then marches
through the cycle together. A second dose landing when the cohort is in S
kills far more than the same dose landing when the cohort is in G1. The
population has a MEMORY of the first dose, and λ alone has no way to
express it.

This module integrates the phase compartments directly:

    dN₁/dt = 2·k_n·N_n − (k₁ + d₁(t))·N₁
    dN_i/dt = k_{i−1}·N_{i−1} − (k_i + d_i(t))·N_i

block by block (exact at the discontinuities, rather than fighting them
inside one solve), so the kill rates can switch at dose boundaries.

What it buys, beyond nicer curves:

  * **It validates — and bounds — the QSS approximation.** For blocks long
    compared with the cycle the transient agrees with
    `combination.evaluate_schedule`; for short blocks they diverge, and
    the transient is the correct one. The disagreement is now MEASURABLE
    rather than merely warned about.
  * **Dose timing becomes optimisable.** `best_second_dose_delay` scans
    the gap between two pulses and reports the one that maximises kill —
    a question the QSS model is structurally unable to answer, since every
    gap gives it the same answer.

Synchrony is quantified as the total-variation distance between the
current phase distribution and the stable (asynchronous) one: 0 means
fully asynchronous, larger means the population is marching in step.

Non-AI: direct integration of the linear compartment ODE. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from src.cell.cycle import CellCycle, stable_phase_fractions
from src.cell.combination import ScheduleBlock, apply_drugs


@dataclass
class PhaseTrajectory:
    """Time course of the phase compartments under a schedule."""

    t_h: List[float]
    counts: Dict[str, List[float]]
    total: List[float]
    synchrony: List[float]
    quiescent_fraction: float
    phases: List[str]

    def total_with_quiescent(self, i: int) -> float:
        q = self.quiescent_fraction
        return (1.0 - q) * self.total[i] + q

    @property
    def log10_kill(self) -> float:
        """Net log-kill at the end of the schedule."""
        return -math.log10(max(self.total_with_quiescent(-1), 1e-300))

    @property
    def peak_synchrony(self) -> float:
        return max(self.synchrony) if self.synchrony else 0.0

    def phase_fractions_at(self, i: int) -> Dict[str, float]:
        tot = self.total[i]
        if tot <= 0:
            return {p: 0.0 for p in self.phases}
        return {p: self.counts[p][i] / tot for p in self.phases}

    def summary(self) -> str:
        return (f"{self.log10_kill:+.2f} log10 kill over {self.t_h[-1]:.0f} h; "
                f"peak synchrony {self.peak_synchrony:.2f} "
                f"(0 = asynchronous)")


def _synchrony_index(fracs: Dict[str, float],
                     stable: Dict[str, float]) -> float:
    """Total-variation distance from the asynchronous distribution."""
    return 0.5 * sum(abs(fracs[p] - stable[p]) for p in stable)


def simulate_schedule_transient(
    cycle: CellCycle,
    blocks: Sequence[ScheduleBlock],
    *,
    points_per_block: int = 200,
    initial_fractions: Optional[Dict[str, float]] = None,
    rtol: float = 1e-9,
    atol: float = 1e-14,
) -> PhaseTrajectory:
    """Integrate the phase compartments through a treatment schedule.

    Starts from the stable (asynchronous) distribution unless
    `initial_fractions` says otherwise, normalised to a total of 1.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    if not blocks:
        raise ValueError("need at least one block")

    phases = cycle.phases
    stable = stable_phase_fractions(cycle)
    if initial_fractions is None:
        y = np.array([stable[p] for p in phases], dtype=float)
    else:
        missing = set(phases) - set(initial_fractions)
        if missing:
            raise ValueError(f"initial_fractions missing {sorted(missing)}")
        y = np.array([float(initial_fractions[p]) for p in phases])
        if y.min() < 0:
            raise ValueError("initial_fractions must be ≥ 0")
        s = y.sum()
        if s <= 0:
            raise ValueError("initial_fractions must not be all zero")
        y = y / s

    t_all: List[float] = [0.0]
    counts: Dict[str, List[float]] = {p: [float(v)] for p, v in zip(phases, y)}
    totals: List[float] = [float(y.sum())]
    sync: List[float] = [_synchrony_index(
        {p: float(v) / y.sum() for p, v in zip(phases, y)}, stable)]

    t_offset = 0.0
    for blk in blocks:
        eff, kills = apply_drugs(cycle, blk.drugs)
        ks = eff.rates_per_h()
        ds = [float(kills.get(p, 0.0)) for p in phases]
        n = len(phases)

        def rhs(t, v):
            out = np.empty(n)
            out[0] = 2.0 * ks[-1] * v[-1] - (ks[0] + ds[0]) * v[0]
            for i in range(1, n):
                out[i] = ks[i - 1] * v[i - 1] - (ks[i] + ds[i]) * v[i]
            return out

        t_eval = np.linspace(0.0, blk.hours, points_per_block)
        sol = solve_ivp(rhs, (0.0, blk.hours), y, method="LSODA",
                        t_eval=t_eval, rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(f"cycle transient failed: {sol.message}")

        for j in range(1, len(t_eval)):        # skip duplicate start point
            col = np.maximum(sol.y[:, j], 0.0)
            tot = float(col.sum())
            t_all.append(t_offset + float(t_eval[j]))
            for p, v in zip(phases, col):
                counts[p].append(float(v))
            totals.append(tot)
            fr = ({p: float(v) / tot for p, v in zip(phases, col)}
                  if tot > 0 else {p: 0.0 for p in phases})
            sync.append(_synchrony_index(fr, stable))

        y = np.maximum(sol.y[:, -1], 0.0)
        t_offset += blk.hours

    return PhaseTrajectory(t_h=t_all, counts=counts, total=totals,
                           synchrony=sync,
                           quiescent_fraction=cycle.quiescent_fraction,
                           phases=phases)


@dataclass
class TimingScan:
    """Kill achieved by a second pulse as a function of the gap after the first."""

    delays_h: List[float]
    log10_kill: List[float]
    best_delay_h: float
    best_log10_kill: float
    worst_delay_h: float
    worst_log10_kill: float

    @property
    def timing_benefit_log10(self) -> float:
        """How much the best timing beats the worst — the value of scheduling."""
        return self.best_log10_kill - self.worst_log10_kill

    def summary(self) -> str:
        return (f"best gap {self.best_delay_h:.1f} h "
                f"({self.best_log10_kill:+.2f} log), worst {self.worst_delay_h:.1f} h "
                f"({self.worst_log10_kill:+.2f} log) — timing is worth "
                f"{self.timing_benefit_log10:.2f} log10")


def best_second_dose_delay(
    cycle: CellCycle,
    drug,
    *,
    pulse_hours: float = 2.0,
    delays_h: Optional[Sequence[float]] = None,
    tail_hours: float = 0.0,
    points_per_block: int = 120,
) -> TimingScan:
    """Scan the gap between two identical pulses for the best total kill.

    A question the quasi-steady-state model cannot answer: with λ·Δt
    bookkeeping every gap gives the same total, because the drug-free gap
    contributes the same regrowth regardless of WHEN it happens. Only the
    cohort structure makes timing matter.
    """
    if pulse_hours <= 0:
        raise ValueError("pulse_hours must be > 0")
    if delays_h is None:
        # Scan a bit more than one cycle so the optimum is inside the window.
        n = 40
        span = 1.3 * cycle.cycle_time_h
        delays_h = [span * i / (n - 1) for i in range(n)]

    results: List[float] = []
    for gap in delays_h:
        blocks = [ScheduleBlock([drug], pulse_hours)]
        if gap > 0:
            blocks.append(ScheduleBlock([], gap))
        blocks.append(ScheduleBlock([drug], pulse_hours))
        if tail_hours > 0:
            blocks.append(ScheduleBlock([], tail_hours))
        tr = simulate_schedule_transient(cycle, blocks,
                                         points_per_block=points_per_block)
        results.append(tr.log10_kill)

    best_i = max(range(len(results)), key=lambda i: results[i])
    worst_i = min(range(len(results)), key=lambda i: results[i])
    return TimingScan(delays_h=list(delays_h), log10_kill=results,
                      best_delay_h=float(delays_h[best_i]),
                      best_log10_kill=results[best_i],
                      worst_delay_h=float(delays_h[worst_i]),
                      worst_log10_kill=results[worst_i])
