"""Combination therapy — why two drugs, and why the ORDER matters.

Single agents fail for two reasons this layer can now express: resistance
(`resistance.py`) and the transit-limited kill ceiling of phase-specific
drugs (`cycle.py`). Combinations are the standard answer to both, but they
are not automatically additive — a badly scheduled pair can be WORSE than
one drug alone.

Two mechanisms, both derived rather than asserted.

1. **Resistance multiplies.** If resistance to A and to B arise
   independently, a cell must acquire both to survive the pair, so the
   double-resistant fraction is the PRODUCT f_A·f_B — two 10⁻⁴ clones give
   10⁻⁸, which a tumour of realistic size simply may not contain. That
   product is the whole quantitative case for combination therapy. Where
   the mechanisms overlap (one efflux pump exporting both drugs), the
   fractions do NOT multiply, and `cross_resistance` interpolates to the
   pessimistic limit min(f_A, f_B).

2. **Cytostatic drugs ANTAGONISE phase-specific ones.** A drug that
   arrests cells in G1 stops them entering S — so an S-phase drug given
   at the same time has nothing to kill. In `cycle.py` terms the arrest
   lowers k_G1, and the phase-specific kill ceiling IS −k_G1. Blocking the
   cycle therefore lowers the ceiling of its partner. The model predicts
   the antagonism from the ceiling result; it is not a special case bolted
   on.

Scheduling is evaluated as piecewise-constant blocks: within a block the
drug set is fixed, the population grows at that block's λ, and log-kills
add. HONESTY LIMIT — this is a QUASI-STEADY-STATE treatment: it assumes
the phase distribution re-equilibrates fast compared with the block
length. Genuine synchronisation effects (giving a block short enough that
the cohort structure itself is the point) need the full transient and are
NOT captured here. Blocks much shorter than a cycle time should not be
trusted.

Non-AI: closed-form compartment kinetics. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from src.cell.cycle import CellCycle, growth_rate_per_h, evaluate_cycle_drug


@dataclass
class DrugAction:
    """What one drug does to a cycling cell.

    `kill_per_h` is a per-phase cytotoxic rate (empty = not cytotoxic).
    `arrest_phase` + `arrest_strength` describe a CYTOSTATIC block: exit
    from that phase is slowed by (1 − arrest_strength).
    """

    name: str
    kill_per_h: Dict[str, float] = field(default_factory=dict)
    arrest_phase: Optional[str] = None
    arrest_strength: float = 0.0
    trust: str = "uncalibrated"

    def __post_init__(self):
        if any(v < 0 for v in self.kill_per_h.values()):
            raise ValueError("kill rates must be ≥ 0")
        if not (0.0 <= self.arrest_strength < 1.0):
            raise ValueError("arrest_strength must be in [0,1)")
        if self.arrest_strength > 0 and not self.arrest_phase:
            raise ValueError("arrest_strength needs an arrest_phase")


def apply_drugs(cycle: CellCycle, drugs: Sequence[DrugAction]):
    """Combine drug actions into an effective cycle + per-phase kill map.

    Arrests multiply (two G1 blocks are more complete than one); kills add
    (independent cytotoxic mechanisms in the same phase).
    """
    hours = dict(cycle.phase_hours)
    kills: Dict[str, float] = {}
    for d in drugs:
        for phase, rate in d.kill_per_h.items():
            if phase not in hours:
                raise ValueError(f"{d.name}: unknown phase {phase!r}")
            kills[phase] = kills.get(phase, 0.0) + rate
        if d.arrest_phase:
            if d.arrest_phase not in hours:
                raise ValueError(f"{d.name}: unknown phase {d.arrest_phase!r}")
            # Slowing exit = lengthening the phase.
            hours[d.arrest_phase] = hours[d.arrest_phase] / \
                (1.0 - d.arrest_strength)
    eff = CellCycle(phase_hours=hours,
                    quiescent_fraction=cycle.quiescent_fraction,
                    source=cycle.source + "+drugs", trust=cycle.trust)
    return eff, kills


def evaluate_combination(cycle: CellCycle, drugs: Sequence[DrugAction]):
    """Growth rate and phase structure under a set of simultaneous drugs."""
    eff, kills = apply_drugs(cycle, drugs)
    return evaluate_cycle_drug(eff, kills)


@dataclass
class ScheduleBlock:
    """A stretch of time during which a fixed drug set is present."""

    drugs: List[DrugAction]
    hours: float

    def __post_init__(self):
        if self.hours <= 0:
            raise ValueError("block hours must be > 0")


@dataclass
class ScheduleOutcome:
    """Result of a piecewise-constant treatment schedule."""

    total_hours: float
    log10_kill: float
    block_rates_per_h: List[float]
    block_labels: List[str]
    quiescent_fraction: float
    quasi_steady_state_ok: bool
    warning: str = ""

    @property
    def surviving_fraction(self) -> float:
        return 10.0 ** (-self.log10_kill)

    def summary(self) -> str:
        seq = " → ".join(f"{lab}({r:+.4f}/h)"
                         for lab, r in zip(self.block_labels,
                                           self.block_rates_per_h))
        w = f"  ⚠ {self.warning}" if self.warning else ""
        return (f"{seq} over {self.total_hours:.0f} h → "
                f"{self.log10_kill:+.2f} log10 kill "
                f"({100*self.surviving_fraction:.3g}% surviving){w}")


def evaluate_schedule(cycle: CellCycle,
                      blocks: Sequence[ScheduleBlock]) -> ScheduleOutcome:
    """Total log-kill of a sequence of treatment blocks.

    Within each block the population grows at that block's λ, so the
    log-kills add. See the module docstring for the quasi-steady-state
    limitation.
    """
    if not blocks:
        raise ValueError("need at least one block")

    rates: List[float] = []
    labels: List[str] = []
    integral = 0.0
    for b in blocks:
        eff, kills = apply_drugs(cycle, b.drugs)
        lam = growth_rate_per_h(eff, kills)
        rates.append(lam)
        labels.append("+".join(d.name for d in b.drugs) or "no drug")
        integral += lam * b.hours

    total_h = sum(b.hours for b in blocks)
    q = cycle.quiescent_fraction
    cycling = math.exp(integral)
    surviving = (1.0 - q) * cycling + q
    log10_kill = -math.log10(max(surviving, 1e-300))

    # Quasi-steady-state check: blocks should be long vs the cycle time.
    shortest = min(b.hours for b in blocks)
    ok = shortest >= cycle.cycle_time_h
    warn = ("" if ok else
            f"block of {shortest:.0f} h is shorter than the {cycle.cycle_time_h:.0f} h "
            "cycle — quasi-steady-state assumption is not reliable; "
            "synchronisation effects are not modelled")

    return ScheduleOutcome(
        total_hours=total_h, log10_kill=log10_kill,
        block_rates_per_h=rates, block_labels=labels,
        quiescent_fraction=q, quasi_steady_state_ok=ok, warning=warn)


def double_resistant_fraction(f_a: float, f_b: float,
                              cross_resistance: float = 0.0) -> float:
    """Fraction of cells resistant to BOTH drugs.

    `cross_resistance` = 0 → independent mechanisms, so the fractions
    MULTIPLY (f_a·f_b) — the quantitative case for combination therapy.
    = 1 → one shared mechanism (e.g. a pump exporting both), so a cell
    resistant to the rarer drug is resistant to both: min(f_a, f_b).
    Intermediate values interpolate.
    """
    for f in (f_a, f_b):
        if not (0.0 <= f <= 1.0):
            raise ValueError("fractions must be in [0,1]")
    if not (0.0 <= cross_resistance <= 1.0):
        raise ValueError("cross_resistance must be in [0,1]")
    independent = f_a * f_b
    shared = min(f_a, f_b)
    return independent + cross_resistance * (shared - independent)
