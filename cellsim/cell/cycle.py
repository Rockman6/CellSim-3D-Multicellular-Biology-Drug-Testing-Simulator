"""Cell-cycle phases — why a phase-specific drug has a kill ceiling.

The fate layer treats a cell as a single state with one net growth rate.
Real cells run a cycle — G1 → S → G2/M → divide — and many important drugs
only act in ONE phase: antimetabolites (5-FU, methotrexate, gemcitabine)
kill cells replicating DNA in S; taxanes and vinca alkaloids kill cells in
mitosis. A cell sitting in G1 simply is not a target.

Model: a linear chain of phases with exponentially-distributed dwell
times, ending in division into two G1 cells. With phase rate k_i = 1/T_i
and a phase-specific kill rate d_i:

    dN₁/dt = 2·k_n·N_n − (k₁ + d₁)·N₁
    dN_i/dt = k_{i−1}·N_{i−1} − (k_i + d_i)·N_i        (i ≥ 2)

Exponential growth N_i(t) = n_i·e^{λt} gives a CHARACTERISTIC EQUATION for
the population growth rate λ:

    2·Π_i [ k_i / (k_i + d_i + λ) ]  =  1

and the stable phase distribution follows from
n_i = n_{i−1}·k_{i−1}/(k_i + d_i + λ).

Two results fall out that a single-rate model cannot express:

  * **The kill ceiling.** Push an S-specific drug's kill rate to infinity
    and λ does NOT go to −∞: it approaches −k_G1, the rate at which cells
    ENTER S. You cannot kill cells faster than they arrive in the phase
    where the drug works. Cytotoxicity beyond that point buys nothing —
    only more exposure TIME does. (A test pins this limit.)
  * **Exponential dwell times grow faster than a fixed cycle length.**
    Variability lets some cells divide early, so λ > ln2/ΣT_i. The model
    says so rather than assuming doubling time = cycle time.

QUIESCENCE (G0). Cells that exit the cycle are refractory to every
cycle-specific drug — a sanctuary in TIME, exactly analogous to the
spatial sanctuary in `tissue.py`. Modelled here as a refractory fraction
that neither divides nor dies, so it sets a FLOOR on what any
cycle-specific therapy can achieve.

Cited: typical mammalian cycle ~24 h with G1 ~11 h, S ~8 h, G2/M ~5 h
(Alberts, Molecular Biology of the Cell; BioNumbers BNID 108374/112260).

Non-AI: closed-form linear compartment kinetics. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# Typical mammalian cycle phase durations (hours) — Alberts, MBoC;
# BioNumbers BNID 108374 / 112260.
_DEFAULT_PHASE_HOURS = {"G1": 11.0, "S": 8.0, "G2M": 5.0}
_LN2 = math.log(2.0)


@dataclass
class CellCycle:
    """Phase durations of a proliferating cell.

    `quiescent_fraction` is the share of the population sitting in G0 —
    refractory to cycle-specific drugs.
    """

    phase_hours: Dict[str, float] = field(
        default_factory=lambda: dict(_DEFAULT_PHASE_HOURS))
    quiescent_fraction: float = 0.0
    source: str = "default:typical-mammalian"
    trust: str = "uncalibrated"

    def __post_init__(self):
        if not self.phase_hours:
            raise ValueError("need at least one phase")
        for name, h in self.phase_hours.items():
            if h <= 0:
                raise ValueError(f"phase {name} duration must be > 0")
        if not (0.0 <= self.quiescent_fraction < 1.0):
            raise ValueError("quiescent_fraction must be in [0,1)")

    @property
    def phases(self) -> List[str]:
        return list(self.phase_hours.keys())

    @property
    def cycle_time_h(self) -> float:
        return sum(self.phase_hours.values())

    def rates_per_h(self) -> List[float]:
        return [1.0 / self.phase_hours[p] for p in self.phases]


def growth_rate_per_h(cycle: CellCycle,
                      kill_per_h: Optional[Dict[str, float]] = None
                      ) -> float:
    """Population growth rate λ solving 2·Π k_i/(k_i + d_i + λ) = 1.

    Solved by bisection on the strictly decreasing function
    f(λ) = 2·Π k_i/(k_i + d_i + λ). λ is bounded below by −min(k_i + d_i)
    (where the product diverges) and grows without bound above.
    """
    ks = cycle.rates_per_h()
    ds = [float((kill_per_h or {}).get(p, 0.0)) for p in cycle.phases]
    if any(d < 0 for d in ds):
        raise ValueError("kill rates must be ≥ 0")

    def f(lam: float) -> float:
        prod = 1.0
        for k, d in zip(ks, ds):
            denom = k + d + lam
            if denom <= 0:
                return math.inf
            prod *= k / denom
        return 2.0 * prod

    # Lower bound: just above the pole at −min(k_i + d_i), where f → +∞.
    pole = -min(k + d for k, d in zip(ks, ds))
    lo = pole + 1e-12
    hi = 1.0
    while f(hi) > 1.0:                 # push up until f drops below 1
        hi *= 2.0
        if hi > 1e9:
            raise RuntimeError("growth rate did not bracket")
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) > 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def stable_phase_fractions(cycle: CellCycle,
                           kill_per_h: Optional[Dict[str, float]] = None
                           ) -> Dict[str, float]:
    """Fraction of the cycling population in each phase at steady growth.

    This is what makes phase-specific drugs finite: only this share of the
    population is a target at any instant.
    """
    lam = growth_rate_per_h(cycle, kill_per_h)
    ks = cycle.rates_per_h()
    ds = [float((kill_per_h or {}).get(p, 0.0)) for p in cycle.phases]
    n = [1.0]
    for i in range(1, len(ks)):
        n.append(n[i - 1] * ks[i - 1] / (ks[i] + ds[i] + lam))
    total = sum(n)
    return {p: v / total for p, v in zip(cycle.phases, n)}


@dataclass
class CycleOutcome:
    """Effect of a (possibly phase-specific) drug on a cycling population."""

    growth_rate_per_h: float
    untreated_growth_rate_per_h: float
    doubling_time_h: Optional[float]
    phase_fractions: Dict[str, float]
    targeted_phases: List[str]
    targeted_fraction: float           # share of cells vulnerable at any instant
    max_kill_rate_per_h: float         # ceiling set by entry into the phase
    quiescent_fraction: float

    @property
    def is_killing(self) -> bool:
        return self.growth_rate_per_h < 0.0

    def cycling_population_at(self, t_h: float) -> float:
        return math.exp(self.growth_rate_per_h * t_h)

    def total_population_at(self, t_h: float) -> float:
        """Cycling cells decay/grow; the quiescent fraction persists."""
        q = self.quiescent_fraction
        return (1.0 - q) * self.cycling_population_at(t_h) + q

    def log10_kill_at(self, t_h: float) -> float:
        return -math.log10(max(self.total_population_at(t_h), 1e-300))

    def summary(self) -> str:
        dt = (f"doubling {self.doubling_time_h:.1f} h"
              if self.doubling_time_h else
              f"decaying (λ={self.growth_rate_per_h:.4f}/h)")
        tgt = ("+".join(self.targeted_phases) if self.targeted_phases
               else "none")
        return (f"target={tgt} ({100*self.targeted_fraction:.0f}% of cells "
                f"vulnerable at any instant); {dt}; "
                f"ceiling λ_min={-self.max_kill_rate_per_h:.4f}/h")


def evaluate_cycle_drug(cycle: CellCycle,
                        kill_per_h: Optional[Dict[str, float]] = None
                        ) -> CycleOutcome:
    """Effect of a phase-specific (or non-specific) drug on the cycle."""
    kills = {k: float(v) for k, v in (kill_per_h or {}).items()}
    unknown = set(kills) - set(cycle.phases)
    if unknown:
        raise ValueError(f"unknown phases in kill_per_h: {sorted(unknown)}")

    lam = growth_rate_per_h(cycle, kills)
    lam0 = growth_rate_per_h(cycle, None)
    fracs = stable_phase_fractions(cycle, kills)
    targeted = [p for p, v in kills.items() if v > 0]

    # Kill ceiling: with an infinitely cytotoxic drug in the targeted
    # phases, the population can decay no faster than cells ENTER the
    # first targeted phase — i.e. the rate of the phase feeding it.
    if targeted:
        idx = [cycle.phases.index(p) for p in targeted]
        ks = cycle.rates_per_h()
        # Entry into targeted phase i comes from phase i−1 (or division).
        feeders = [ks[(i - 1) % len(ks)] for i in idx]
        ceiling = min(feeders)
    else:
        ceiling = 0.0

    return CycleOutcome(
        growth_rate_per_h=lam,
        untreated_growth_rate_per_h=lam0,
        doubling_time_h=(_LN2 / lam) if lam > 0 else None,
        phase_fractions=fracs,
        targeted_phases=targeted,
        targeted_fraction=sum(fracs[p] for p in targeted) if targeted else 0.0,
        max_kill_rate_per_h=ceiling,
        quiescent_fraction=cycle.quiescent_fraction,
    )
