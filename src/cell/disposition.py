"""Full-system transient: the whole composed disposition over time.

`dynamics.py` integrates permeation + binding. `steady_state.py` composes
permeation + pH partitioning + efflux + sink but only at t→∞. This closes
the gap: the transient of the FULL composed system, integrated with a
stiff-aware adaptive solver, whose steady state PROVABLY equals
`solve_steady_state` (the cross-check test).

We track the TOTAL intracellular drug T and recover the free concentration
f = free_from_total(T, sink) each step (the nonspecific sink is fast, so
it is in instantaneous equilibrium with the free pool). The neutral free
form permeates and the pump effluxes the free form:

    dT/dt = k_perm·(C_out·f_neu(pH_out) − f·f_neu(pH_in)) − V_max·f/(K_m+f)

At steady state dT/dt = 0 reproduces exactly the `solve_steady_state`
balance, so the ODE must converge to it — the anchoring check.

Physics the transient makes visible: the binding sink is a CAPACITOR. It
does not change the steady-state free concentration (see `steady_state`),
but it stores drug, so it SLOWS the approach — a bigger sink means a
longer time to reach the same steady state.

Non-AI: closed-form fluxes + LSODA. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

from src.cell.compartments import Permeability, CellGeometry, equilibration_tau_s
from src.cell.partitioning import neutral_fraction
from src.cell.transport import EffluxPump
from src.cell.binding_sink import BindingSink, free_from_total
from src.cell.occupancy import hill_occupancy
from src.cell.steady_state import solve_steady_state


@dataclass
class DispositionTimeCourse:
    """Trajectory of the full composed intracellular disposition."""

    t_s: List[float]
    free_M: List[float]
    total_M: List[float]
    occupancy: Optional[List[float]]
    C_out_M: float
    tau_perm_s: float
    steady_state_free_M: float          # analytic (solve_steady_state)
    converged_free_M: float             # ODE final
    inputs_calibrated: bool
    limiting_trust: str

    def time_to_fraction_of_steady_state(self, frac: float = 0.5
                                         ) -> Optional[float]:
        """First time free reaches `frac`·steady-state free (None if never)."""
        thr = frac * self.steady_state_free_M
        for t, f in zip(self.t_s, self.free_M):
            if f >= thr:
                return t
        return None

    def summary(self) -> str:
        t50 = self.time_to_fraction_of_steady_state(0.5)
        t50s = f"{t50:.0f}s" if t50 is not None else "n/a"
        cal = "calibrated" if self.inputs_calibrated else \
            f"UNCALIBRATED [{self.limiting_trust}]"
        return (f"steady free {self.steady_state_free_M:.2e} M "
                f"(ODE {self.converged_free_M:.2e}), t½≈{t50s}, "
                f"τ_perm={self.tau_perm_s:.0f}s — {cal}")


def simulate_disposition(
    C_out_M: float, *,
    permeability: Permeability,
    geometry: CellGeometry,
    prior=None,
    pKa: Optional[float] = None,
    ion_type: str = "neutral",
    pH_in: float = 7.2,
    pH_out: float = 7.4,
    pump: Optional[EffluxPump] = None,
    sink: Optional[BindingSink] = None,
    t_end_s: Optional[float] = None,
    n_points: int = 200,
    rtol: float = 1e-8,
    atol: float = 1e-18,
) -> DispositionTimeCourse:
    """Integrate the full composed disposition; converges to solve_steady_state.

    `prior` (optional Hill prior) adds a target-occupancy trace at the free
    concentration. All other args match `solve_steady_state`.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

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

    # Analytic steady state (the ODE must converge to this).
    ss = solve_steady_state(
        C_out_M, permeability=permeability, geometry=geometry, pKa=pKa,
        ion_type=ion_type, pH_in=pH_in, pH_out=pH_out, pump=pump, sink=sink)

    # Horizon: permeation τ inflated by the sink's buffering capacity β
    # (a bigger sink slows the approach). Generous multiple ensures we
    # reach steady state.
    beta = 1.0 + (sink.capacity_M / sink.Kd_M if sink else 0.0)
    if t_end_s is None:
        t_end_s = 8.0 * tau * max(1.0, beta)

    def free_of_total(T):
        if sink is None:
            return T
        return free_from_total(T, sink).free_M

    def rhs(t, y):
        T = max(0.0, y[0])
        f = free_of_total(T)
        influx = k * (C_out_M * fn_out - f * fn_in)
        efflux = Vmax * f / (Km + f) if Vmax > 0 else 0.0
        return (influx - efflux,)

    t_eval = np.linspace(0.0, t_end_s, n_points)
    sol = solve_ivp(rhs, (0.0, t_end_s), (0.0,), method="LSODA",
                    t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")

    total = np.maximum(sol.y[0], 0.0)
    free = np.array([free_of_total(float(T)) for T in total])

    occ = None
    if prior is not None:
        Kd = float(prior.parameters["Kd_M"])
        n = float(prior.parameters.get("n_hill", 1.0))
        occ = [float(hill_occupancy(float(f), Kd, n)) for f in free]

    return DispositionTimeCourse(
        t_s=[float(x) for x in t_eval],
        free_M=[float(x) for x in free],
        total_M=[float(x) for x in total],
        occupancy=occ,
        C_out_M=C_out_M,
        tau_perm_s=tau,
        steady_state_free_M=ss.free_in_M,
        converged_free_M=float(free[-1]),
        inputs_calibrated=ss.inputs_calibrated,
        limiting_trust=ss.limiting_trust,
    )
