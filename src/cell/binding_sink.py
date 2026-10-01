"""Intracellular binding sink — free vs total drug inside the cell.

Occupancy, efflux and permeation are all driven by the FREE drug
concentration, but most of the drug inside a cell is often BOUND — to
abundant nonspecific partners (serum/cytosolic protein, membrane lipid,
high-copy off-targets). That reservoir buffers the free concentration: a
cell can hold a large TOTAL amount while the free level that actually
reaches the target stays low, and the reservoir also slows apparent efflux
(it keeps refilling the free pool).

Physics (one saturable sink at capacity S_tot, dissociation K_s):

    total = free + S_tot · free / (K_s + free)

Given total, free is the positive root of  f² + (K_s + S_tot − total)·f
− total·K_s = 0. The free fraction f_u = free/total falls as the sink
capacity rises or its affinity tightens.

Limits (tested): S_tot = 0 → free = total (no sink); a large, tight sink
→ free ≪ total. Saturable: once free ≫ K_s the sink is full and the free
fraction climbs back toward 1.

Provenance: sink capacity + affinity are inputs carrying source + trust
(nonspecific-binding parameters would come from measured fu or a
partition model). Non-AI: closed-form mass-action buffering.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class BindingSink:
    """A saturable intracellular reservoir that binds free drug."""

    capacity_M: float          # S_tot — total sink sites
    Kd_M: float                # K_s — sink dissociation constant
    name: str = "nonspecific binding"
    source: str = "input"
    trust: str = "uncalibrated"

    def __post_init__(self):
        if self.capacity_M < 0 or self.Kd_M <= 0:
            raise ValueError("capacity ≥ 0 and Kd > 0 required")

    def bound_at(self, free_M: float) -> float:
        """Sink-bound drug concentration at a given free concentration."""
        return self.capacity_M * free_M / (self.Kd_M + free_M)


@dataclass
class SinkResult:
    free_M: float
    total_M: float
    bound_M: float
    free_fraction: float           # free / total (the "fu" the target sees)
    sink_saturated_fraction: float  # bound / capacity
    sink_trust: str

    def summary(self) -> str:
        return (f"free {self.free_M:.2e} M of total {self.total_M:.2e} M "
                f"(fu={self.free_fraction:.3f}); sink "
                f"{100*self.sink_saturated_fraction:.0f}% full "
                f"[{self.sink_trust}]")


def total_from_free(free_M: float, sink: BindingSink) -> SinkResult:
    """Total intracellular drug given the free concentration (explicit)."""
    if free_M < 0:
        raise ValueError("free_M must be ≥ 0")
    bound = sink.bound_at(free_M)
    total = free_M + bound
    fu = (free_M / total) if total > 0 else 1.0
    sat = (bound / sink.capacity_M) if sink.capacity_M > 0 else 0.0
    return SinkResult(free_M=free_M, total_M=total, bound_M=bound,
                      free_fraction=fu, sink_saturated_fraction=sat,
                      sink_trust=sink.trust)


def free_from_total(total_M: float, sink: BindingSink) -> SinkResult:
    """Free drug given total intracellular drug (positive root of quadratic)."""
    if total_M < 0:
        raise ValueError("total_M must be ≥ 0")
    if total_M == 0.0:
        return SinkResult(0.0, 0.0, 0.0, 1.0, 0.0, sink.trust)
    Ks, S = sink.Kd_M, sink.capacity_M
    b = Ks + S - total_M
    disc = b * b + 4.0 * total_M * Ks        # > 0
    free = (-b + math.sqrt(disc)) / 2.0
    free = max(0.0, min(free, total_M))
    return total_from_free(free, sink)
