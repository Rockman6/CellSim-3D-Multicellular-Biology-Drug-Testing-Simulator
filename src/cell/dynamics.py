"""Time-resolved kinetics: drug crosses the membrane, then binds a target.

The other cell modules give equilibrium answers. This one integrates the
full coupled ODE so you can watch the approach — and it is built so its
steady state PROVABLY reduces to those equilibria (a consistency check,
`tests/cell/test_dynamics_smoke.py`).

Kinetics from thermodynamics
----------------------------
A binding ΔG (→ K_d prior from `src/bridge`) fixes only the RATIO of the
rate constants:  K_d = k_off / k_on.  We take the association rate k_on
as a physics input (diffusion-limited, capped below the Berg–Purcell /
Smoluchowski limit by an orientational-steric factor — see
`docs/cellsim_full_reference.tex`) and DERIVE k_off = k_on · K_d.

Crucial consequence for honesty: the EQUILIBRIUM (t→∞ occupancy) depends
only on K_d and is therefore exactly as trustworthy as the K_d prior. k_on
sets only the *timescale* to reach it, so its uncertainty is flagged
separately (`kinetics_calibrated`) and never contaminates the equilibrium
trust.

Coupled system (buffered extracellular reservoir L_out held constant):

    dL_in/dt = k_perm·(L_out − L_in) − (k_on·L_in·T − k_off·C)
    dT/dt    = −(k_on·L_in·T − k_off·C)
    dC/dt    = +(k_on·L_in·T − k_off·C)

with k_perm = P·A/V = 1/τ (from `compartments`), T + C = T_tot conserved.

The system is STIFF (fast binding, slow permeation), so it is integrated
with LSODA (adaptive, stiff-aware), not a fixed-step hand-rolled Euler.

Non-AI: closed-form rate law + a standard ODE solver. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

from src.cell.compartments import (
    Permeability, CellGeometry, spherical_cell_geometry, equilibration_tau_s)
from src.cell.occupancy import hill_occupancy

# Representative protein–ligand association rate (M⁻¹ s⁻¹). Labelled
# project convention: measured small-molecule on-rates cluster around
# 1e6–1e7, well below the ~1e9 Smoluchowski diffusion limit because of
# the orientational/steric factor. Overridable per call. The equilibrium
# is INDEPENDENT of this value; it only sets the kinetic timescale.
_DEFAULT_KON_PER_M_PER_S = 1.0e6


@dataclass
class RateConstants:
    """k_on / k_off / K_d with provenance. Equilibrium uses only K_d."""

    kon_per_M_per_s: float
    koff_per_s: float
    Kd_M: float
    kon_source: str = "convention:representative-1e6"
    kinetics_calibrated: bool = False   # True only if k_on is measured
    equilibrium_trust: str = "uncalibrated"
    accuracy_known: bool = False


def rate_constants_from_prior(prior,
                              kon_per_M_per_s: float = _DEFAULT_KON_PER_M_PER_S,
                              kon_source: Optional[str] = None,
                              kinetics_calibrated: bool = False
                              ) -> RateConstants:
    """Derive (k_on, k_off, K_d) from a Hill `RateLawPrior`.

    k_off = k_on · K_d. The equilibrium trust comes from the prior; the
    kinetic calibration is a separate flag (default False — the default
    k_on is a convention, not a measurement).
    """
    if getattr(prior, "type", None) != "hill":
        raise ValueError("rate_constants_from_prior needs a Hill prior")
    if kon_per_M_per_s <= 0:
        raise ValueError("kon must be > 0")
    Kd = float(prior.parameters["Kd_M"])
    return RateConstants(
        kon_per_M_per_s=kon_per_M_per_s,
        koff_per_s=kon_per_M_per_s * Kd,
        Kd_M=Kd,
        kon_source=kon_source or "convention:representative-1e6",
        kinetics_calibrated=kinetics_calibrated,
        equilibrium_trust=getattr(prior, "trust", "uncalibrated"),
        accuracy_known=getattr(prior, "accuracy_known", False),
    )


@dataclass
class TimeCourse:
    """Integrated trajectory of intracellular drug binding.

    Arrays are aligned with `t_s`. `equilibrium_occupancy` is the analytic
    t→∞ value (Hill at L_out); `converged_occupancy` is the ODE's final
    value — the two agreeing is the module's correctness check.
    """

    t_s: List[float]
    L_in_M: List[float]
    free_target_M: List[float]
    complex_M: List[float]
    occupancy: List[float]                 # C / T_tot over time
    T_total_M: float
    L_out_M: float
    rates: RateConstants
    tau_perm_s: float
    equilibrium_occupancy: float
    converged_occupancy: float
    max_mass_imbalance: float              # sup|T+C−T_tot|/T_tot (should be ~0)
    notes: str = ""

    @property
    def equilibrium_decision_grade(self) -> bool:
        """The EQUILIBRIUM is decision-grade on K_d alone (k_on-independent)."""
        return self.rates.accuracy_known and \
            self.rates.equilibrium_trust == "trustworthy_absolute"

    def occupancy_at(self, t_query_s: float) -> float:
        """Linear-interpolated occupancy at an arbitrary time."""
        ts = self.t_s
        if t_query_s <= ts[0]:
            return self.occupancy[0]
        if t_query_s >= ts[-1]:
            return self.occupancy[-1]
        # ts is monotonically increasing.
        lo, hi = 0, len(ts) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if ts[mid] <= t_query_s:
                lo = mid
            else:
                hi = mid
        frac = (t_query_s - ts[lo]) / (ts[hi] - ts[lo])
        return self.occupancy[lo] + frac * (self.occupancy[hi] - self.occupancy[lo])

    def time_to_occupancy(self, target_fraction_of_eq: float = 0.5
                          ) -> Optional[float]:
        """First time occupancy reaches `frac`·equilibrium (None if never)."""
        thr = target_fraction_of_eq * self.equilibrium_occupancy
        for i, occ in enumerate(self.occupancy):
            if occ >= thr:
                return self.t_s[i]
        return None

    def summary(self) -> str:
        grade = ("equilibrium decision-grade"
                 if self.equilibrium_decision_grade else
                 f"equilibrium NOT decision-grade [{self.rates.equilibrium_trust}]")
        kin = ("k_on measured" if self.rates.kinetics_calibrated
               else "k_on=convention → TIMESCALE uncalibrated")
        t50 = self.time_to_occupancy(0.5)
        t50s = f"{t50:.0f}s" if t50 is not None else "n/a"
        return (f"eq occupancy {100*self.equilibrium_occupancy:.1f}% "
                f"(ODE {100*self.converged_occupancy:.1f}%), t½≈{t50s}, "
                f"τ_perm={self.tau_perm_s:.0f}s — {grade}; {kin}")


def simulate_binding_in_cell(
    prior,
    *,
    target_conc_M: float,
    drug_out_M: float,
    permeability: Permeability,
    geometry: Optional[CellGeometry] = None,
    t_end_s: Optional[float] = None,
    n_points: int = 200,
    kon_per_M_per_s: float = _DEFAULT_KON_PER_M_PER_S,
    kinetics_calibrated: bool = False,
    rtol: float = 1e-8,
    atol: float = 1e-16,
) -> TimeCourse:
    """Integrate drug permeation + intracellular target binding over time.

    Drug is applied OUTSIDE at `drug_out_M` (buffered reservoir); it
    permeates in and reversibly binds an intracellular target present at
    `target_conc_M`. Returns the full trajectory plus the analytic
    equilibrium the ODE must converge to.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    if target_conc_M < 0 or drug_out_M < 0:
        raise ValueError("concentrations must be ≥ 0")
    if geometry is None:
        geometry = spherical_cell_geometry()

    rates = rate_constants_from_prior(
        prior, kon_per_M_per_s=kon_per_M_per_s,
        kinetics_calibrated=kinetics_calibrated)
    tau = equilibration_tau_s(permeability, geometry)
    k_perm = 1.0 / tau
    kon, koff = rates.kon_per_M_per_s, rates.koff_per_s
    L_out, T_tot = drug_out_M, target_conc_M

    # Default horizon: several of the SLOWEST timescale so we reach steady
    # state. Slowest of permeation (τ) and binding relaxation (~1/koff).
    if t_end_s is None:
        slow = max(tau, 1.0 / koff if koff > 0 else tau)
        t_end_s = 8.0 * slow

    def rhs(t, y):
        L_in, T, C = y
        binding = kon * L_in * T - koff * C
        return (k_perm * (L_out - L_in) - binding, -binding, binding)

    t_eval = np.linspace(0.0, t_end_s, n_points)
    sol = solve_ivp(rhs, (0.0, t_end_s), (0.0, T_tot, 0.0),
                    method="LSODA", t_eval=t_eval,
                    rtol=rtol, atol=atol, dense_output=False)
    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")

    L_in = sol.y[0]
    T = sol.y[1]
    C = sol.y[2]
    occ = (C / T_tot) if T_tot > 0 else np.zeros_like(C)

    # Analytic equilibrium: at steady state L_in → L_out and C/T_tot is the
    # Hill occupancy at L_out (n=1). This is what the ODE must converge to.
    eq_occ = hill_occupancy(L_out, rates.Kd_M, 1.0) if T_tot > 0 else 0.0
    mass_imbalance = float(np.max(np.abs((T + C) - T_tot)) / T_tot) \
        if T_tot > 0 else 0.0

    return TimeCourse(
        t_s=[float(x) for x in t_eval],
        L_in_M=[float(x) for x in L_in],
        free_target_M=[float(x) for x in T],
        complex_M=[float(x) for x in C],
        occupancy=[float(x) for x in occ],
        T_total_M=T_tot,
        L_out_M=L_out,
        rates=rates,
        tau_perm_s=tau,
        equilibrium_occupancy=float(eq_occ),
        converged_occupancy=float(occ[-1]),
        max_mass_imbalance=mass_imbalance,
    )
