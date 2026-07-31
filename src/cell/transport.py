"""Active efflux transport (e.g. P-glycoprotein) — pumping the drug back out.

Passive permeation and pH trapping both drive the drug toward (or past)
C_out. An efflux pump does the opposite: it spends energy to move the drug
OUT of the cell, holding the intracellular concentration BELOW C_out so it
can't reach an intracellular target. This is the central mechanism of
multidrug resistance (P-gp / ABCB1, BCRP, MRPs).

The pump is SATURABLE (Michaelis-Menten), which is why it matters clinically:
at low drug it clamps C_in near zero, but a high enough dose overwhelms it.

Steady state balances passive net influx against the pump efflux:

    k_perm·(C_out − C_in)  =  V_max · C_in / (K_m + C_in)

with k_perm = P·A/V = 1/τ (from `compartments`). Rearranged this is a
quadratic in C_in with one physical (positive) root:

    k_perm·C_in² + (k_perm·K_m + V_max − k_perm·C_out)·C_in
                 − k_perm·C_out·K_m = 0

Limits (checked in tests): V_max→0 gives C_in→C_out (no pump); a strong
pump at low dose gives C_in ≪ C_out; a saturating dose (C_out ≫ K_m)
lets C_in recover toward C_out.

Provenance: the pump kinetics (V_max, K_m) are inputs carrying source +
trust. K_m is the pump's affinity for the drug and can come from the
`src/bridge.affinity_to_michaelis` prior; V_max = turnover × pump copies
/ cell volume. Not from our own physics yet → `uncalibrated` by default.

Non-AI: closed-form Michaelis-Menten balance. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


@dataclass
class EffluxPump:
    """Saturable efflux-pump kinetics (concentration-rate units).

    V_max_M_per_s is the maximal efflux expressed as intracellular
    concentration cleared per second (turnover × copies / cell volume).
    """

    Vmax_M_per_s: float
    Km_M: float
    name: str = "efflux pump"
    source: str = "input"
    trust: str = "uncalibrated"        # 'calibrated' only if measured/derived

    def __post_init__(self):
        if self.Vmax_M_per_s < 0 or self.Km_M <= 0:
            raise ValueError("Vmax must be ≥ 0 and Km must be > 0")

    def efflux_rate(self, C_in_M: float) -> float:
        """Michaelis-Menten efflux rate at intracellular concentration C_in."""
        return self.Vmax_M_per_s * C_in_M / (self.Km_M + C_in_M)


@dataclass
class EffluxSteadyState:
    """Intracellular steady state under passive influx + active efflux."""

    C_in_M: float
    C_out_M: float
    accumulation_ratio: float          # C_in / C_out (< 1 when pump wins)
    fold_reduction_vs_passive: float   # C_out / C_in (how much the pump lowers it)
    pump_saturated_fraction: float     # efflux / Vmax at steady state
    tau_perm_s: float
    pump_trust: str
    residual: float                    # |influx − efflux| at the solution (~0)

    def summary(self) -> str:
        return (f"C_in/C_out = {self.accumulation_ratio:.3f} "
                f"({self.fold_reduction_vs_passive:.1f}× lower than passive); "
                f"pump {100*self.pump_saturated_fraction:.0f}% saturated "
                f"[{self.pump_trust}]")


def efflux_steady_state(C_out_M: float, permeability, geometry,
                        pump: EffluxPump) -> EffluxSteadyState:
    """Steady-state intracellular concentration under passive influx + efflux.

    Args:
        C_out_M: extracellular concentration (buffered reservoir).
        permeability, geometry: as in `compartments` (give k_perm = 1/τ).
        pump: the efflux-pump kinetics.

    Solves the influx = efflux balance in closed form (positive root of the
    quadratic).
    """
    from src.cell.compartments import equilibration_tau_s

    if C_out_M < 0:
        raise ValueError("C_out_M must be ≥ 0")
    tau = equilibration_tau_s(permeability, geometry)
    k = 1.0 / tau
    Vmax, Km = pump.Vmax_M_per_s, pump.Km_M

    if C_out_M == 0.0:
        C_in = 0.0
    else:
        a = k
        b = k * Km + Vmax - k * C_out_M
        c = -k * C_out_M * Km
        disc = b * b - 4.0 * a * c          # > 0 since c < 0
        C_in = (-b + math.sqrt(disc)) / (2.0 * a)
        C_in = max(0.0, min(C_in, C_out_M))  # physical bounds

    influx = k * (C_out_M - C_in)
    efflux = pump.efflux_rate(C_in)
    ratio = (C_in / C_out_M) if C_out_M > 0 else 1.0
    fold = (C_out_M / C_in) if C_in > 0 else math.inf
    sat = (efflux / Vmax) if Vmax > 0 else 0.0

    return EffluxSteadyState(
        C_in_M=C_in,
        C_out_M=C_out_M,
        accumulation_ratio=ratio,
        fold_reduction_vs_passive=fold,
        pump_saturated_fraction=sat,
        tau_perm_s=tau,
        pump_trust=pump.trust,
        residual=abs(influx - efflux),
    )


def pump_from_michaelis_prior(prior, pump_copies: float,
                              cell_volume_L: float,
                              name: str = "efflux pump") -> EffluxPump:
    """Build an `EffluxPump` from a `bridge.affinity_to_michaelis` prior.

    V_max = kcat · pump_copies / (N_A · cell_volume) expressed as M/s; K_m
    is the pump's substrate affinity from the prior. Trust follows the
    prior (uncalibrated unless the prior itself is calibrated).
    """
    if getattr(prior, "type", None) != "michaelis":
        raise ValueError("pump_from_michaelis_prior needs a michaelis prior")
    _AVOGADRO = 6.02214076e23
    kcat = float(prior.parameters["kcat_per_s"])
    Km = float(prior.parameters["KM_M"])
    if pump_copies < 0 or cell_volume_L <= 0:
        raise ValueError("pump_copies ≥ 0 and cell_volume_L > 0 required")
    Vmax = kcat * pump_copies / (_AVOGADRO * cell_volume_L)   # mol/L/s = M/s
    trust = "calibrated" if getattr(prior, "accuracy_known", False) else \
        "uncalibrated"
    return EffluxPump(Vmax_M_per_s=Vmax, Km_M=Km, name=name,
                      source=f"michaelis:{getattr(prior, 'method', '')}",
                      trust=trust)
