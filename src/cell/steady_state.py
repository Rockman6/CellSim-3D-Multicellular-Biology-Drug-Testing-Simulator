"""Unified intracellular steady state — all transport effects composed.

Each transport module is exact alone. This composes them into ONE
self-consistent steady state: passive permeation of the neutral form (pH
partitioning), active efflux (saturable pump), and an intracellular
binding sink — then reports the FREE concentration that drives target
occupancy and the TOTAL concentration stored.

The membrane-permeant species is the neutral, unbound drug. At steady
state the passive neutral flux balances the pump efflux:

    k_perm·(C_out·f_neu(pH_out) − f·f_neu(pH_in))  =  V_max·f/(K_m + f)

where f is the free intracellular concentration (all protonation states,
unbound to the sink). This is a quadratic in f with one physical root.

Physics insight the composition makes explicit: **the binding sink does
NOT change the steady-state free concentration.** At steady state the sink
is in equilibrium with f, so there is no net flux into it; it only sets
how much TOTAL drug is stored (free + bound) and how long the approach
takes. So f is fixed by permeation + partitioning + efflux; the sink acts
on the total and the transient.

Cross-checks (tests) reduce this to each standalone module in its limit:
  * no pump, no pH gradient          → f = C_out            (permeation)
  * no pump                          → f = C_out·R          (partitioning)
  * no pH gradient, no sink          → matches efflux_steady_state
  * f then feeds occupancy / sink exactly as the standalone modules.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from src.cell.compartments import Permeability, CellGeometry, equilibration_tau_s
from src.cell.partitioning import neutral_fraction
from src.cell.transport import EffluxPump
from src.cell.binding_sink import BindingSink
from src.cell.occupancy import hill_occupancy

_TRUST_RANK = {
    "trustworthy_absolute": 3, "calibrated": 3,
    "rank_order_only": 2, "do_not_trust_absolute": 1, "uncalibrated": 0,
}


@dataclass
class IntracellularSteadyState:
    """Composed steady state: free + total intracellular concentration."""

    free_in_M: float               # drives occupancy, efflux, permeation
    total_in_M: float              # free + sink-bound
    sink_bound_M: float
    C_out_M: float
    free_accumulation: float       # free_in / C_out
    total_accumulation: float      # total_in / C_out
    tau_perm_s: float
    residual: float                # |neutral influx − efflux| at solution
    inputs_calibrated: bool
    limiting_trust: str

    def occupancy(self, prior) -> float:
        """Target occupancy at the FREE intracellular concentration."""
        Kd = float(prior.parameters["Kd_M"])
        n = float(prior.parameters.get("n_hill", 1.0))
        return hill_occupancy(self.free_in_M, Kd, n)

    def summary(self) -> str:
        return (f"free {self.free_in_M:.2e} M ({self.free_accumulation:.2f}× C_out), "
                f"total {self.total_in_M:.2e} M ({self.total_accumulation:.1f}× C_out) "
                f"— inputs {'calibrated' if self.inputs_calibrated else 'UNCALIBRATED'}")


def solve_steady_state(
    C_out_M: float, *,
    permeability: Permeability,
    geometry: CellGeometry,
    pKa: Optional[float] = None,
    ion_type: str = "neutral",
    pH_in: float = 7.2,
    pH_out: float = 7.4,
    pump: Optional[EffluxPump] = None,
    sink: Optional[BindingSink] = None,
) -> IntracellularSteadyState:
    """Composed intracellular steady state (permeation + pH + efflux + sink).

    Args:
        C_out_M: extracellular (free) concentration.
        permeability, geometry: give k_perm = 1/τ.
        pKa, ion_type, pH_in, pH_out: pH partitioning of the neutral form.
            ion_type='neutral' disables partitioning (f_neu ≡ 1).
        pump: optional efflux pump (default: none).
        sink: optional intracellular binding sink (affects TOTAL only).
    """
    if C_out_M < 0:
        raise ValueError("C_out_M must be ≥ 0")
    tau = equilibration_tau_s(permeability, geometry)
    k = 1.0 / tau

    if ion_type == "neutral" or pKa is None:
        fn_in = fn_out = 1.0
    else:
        fn_in = neutral_fraction(pH_in, pKa, ion_type)
        fn_out = neutral_fraction(pH_out, pKa, ion_type)

    Vmax = pump.Vmax_M_per_s if pump else 0.0
    Km = pump.Km_M if pump else 1.0

    # Balance: k·(C_out·fn_out − f·fn_in) = Vmax·f/(Km+f). Quadratic in f.
    A = k * C_out_M * fn_out
    B = k * fn_in
    if C_out_M == 0.0:
        free = 0.0
    elif Vmax == 0.0:
        free = A / B if B > 0 else 0.0            # f = C_out·fn_out/fn_in
    else:
        a = B
        b = B * Km + Vmax - A
        c = -A * Km
        disc = b * b - 4.0 * a * c                 # > 0 since c < 0
        free = (-b + math.sqrt(disc)) / (2.0 * a)
        free = max(0.0, free)

    bound = sink.bound_at(free) if sink else 0.0
    total = free + bound

    influx = k * (C_out_M * fn_out - free * fn_in)
    efflux = Vmax * free / (Km + free) if Vmax > 0 else 0.0

    trusts = ["calibrated" if permeability.trust == "calibrated" else "uncalibrated"]
    if pump:
        trusts.append(pump.trust)
    if sink:
        trusts.append(sink.trust)
    if ion_type != "neutral" and pKa is not None:
        trusts.append("uncalibrated")   # pKa/pH provenance not asserted here
    limiting = min(trusts, key=lambda t: _TRUST_RANK.get(t, 0))
    calibrated = all(_TRUST_RANK.get(t, 0) >= 3 for t in trusts)

    return IntracellularSteadyState(
        free_in_M=free,
        total_in_M=total,
        sink_bound_M=bound,
        C_out_M=C_out_M,
        free_accumulation=(free / C_out_M) if C_out_M > 0 else 1.0,
        total_accumulation=(total / C_out_M) if C_out_M > 0 else 1.0,
        tau_perm_s=tau,
        residual=abs(influx - efflux),
        inputs_calibrated=calibrated,
        limiting_trust=limiting,
    )
