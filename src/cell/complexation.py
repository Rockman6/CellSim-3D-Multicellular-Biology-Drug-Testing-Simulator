"""Drug–drug complexation: two drugs bind EACH OTHER to form a new species.

Answers the literal question "if two drugs stick together and form a new
thing, how much of it is there, and how strongly do they stick?". This is
distinct from competition (`competition.py`, where two drugs contend for a
target site). Here the two drugs are the binding partners:

    A + B  ⇌  AB,     K_d = [A][B] / [AB]

The "how strongly they stick" ruler is K_d for the A·B pair — a binding
free energy like any other, so it comes from the same `src/bridge` prior
(dock/FEP A against B). Given total (formulated) concentrations A_tot and
B_tot, the free monomer and complex concentrations solve the standard
bimolecular quadratic:

    [AB]² − (A_tot + B_tot + K_d)·[AB] + A_tot·B_tot = 0
    [AB] = ½[(A_tot+B_tot+K_d) − √((A_tot+B_tot+K_d)² − 4·A_tot·B_tot)]

(the root ≤ min(A_tot, B_tot); the other root is unphysical). [AB] is
monotonically DECREASING in K_d, so its CI comes exactly from the prior's
K_d CI endpoints — same rigorous interval approach as the other cell
modules. Trust rides through.

Not modelled: higher-order aggregates (A₂B, ABₙ), covalent adduct
formation (irreversible), or the new species' OWN target binding (that is
just another prior fed back into the occupancy model).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple


def _complex_conc(A_tot: float, B_tot: float, Kd: float) -> float:
    """Physical root of the bimolecular binding quadratic for [AB]."""
    b = A_tot + B_tot + Kd
    disc = b * b - 4.0 * A_tot * B_tot
    disc = max(disc, 0.0)                     # guard tiny negative round-off
    ab = 0.5 * (b - math.sqrt(disc))
    # Clamp to the physical range [0, min(A_tot, B_tot)].
    return max(0.0, min(ab, A_tot, B_tot))


@dataclass
class ComplexationResult:
    """Equilibrium speciation of A + B ⇌ AB."""

    Kd_M: float
    A_total_M: float
    B_total_M: float
    complex_M: float                          # [AB]
    free_A_M: float
    free_B_M: float
    complexed_fraction_A: float               # [AB]/A_tot
    complex_ci95: Optional[Tuple[float, float]] = None
    trust: str = "uncalibrated"
    accuracy_known: bool = False

    @property
    def decision_grade(self) -> bool:
        return self.accuracy_known and self.trust == "trustworthy_absolute"

    def summary(self) -> str:
        fa = 100.0 * self.complexed_fraction_A
        grade = "decision-grade" if self.decision_grade else \
            f"NOT decision-grade [{self.trust}]"
        ci = ""
        if self.complex_ci95 is not None:
            lo, hi = self.complex_ci95
            ci = f" (95% CI {lo:.2e}–{hi:.2e} M)"
        return (f"[AB]={self.complex_M:.2e} M{ci}; "
                f"{fa:.1f}% of A is complexed; "
                f"free A={self.free_A_M:.2e}, free B={self.free_B_M:.2e} M "
                f"— {grade}")


def complex_equilibrium(prior_AB, A_total_M: float,
                        B_total_M: float) -> ComplexationResult:
    """Equilibrium of A + B ⇌ AB given the pair's binding prior.

    Args:
        prior_AB: a Hill `RateLawPrior` whose K_d is the A·B affinity.
        A_total_M, B_total_M: total (formulated) concentrations.

    The complex CI (when the prior carries a K_d CI) is exact: [AB] is
    monotonically decreasing in K_d, so evaluate at the K_d CI endpoints.
    """
    if getattr(prior_AB, "type", None) != "hill":
        raise ValueError("complex_equilibrium needs a Hill prior")
    if A_total_M < 0 or B_total_M < 0:
        raise ValueError("total concentrations must be ≥ 0")

    params = prior_AB.parameters
    Kd = float(params["Kd_M"])

    ab = _complex_conc(A_total_M, B_total_M, Kd)
    free_A = A_total_M - ab
    free_B = B_total_M - ab
    frac_A = (ab / A_total_M) if A_total_M > 0 else 0.0

    complex_ci = None
    kd_ci = (getattr(prior_AB, "parameter_ci95", None) or {}).get("Kd_M")
    if kd_ci is not None:
        kd_lo, kd_hi = kd_ci
        # [AB] decreasing in K_d: low K_d → MORE complex.
        ab_hi = _complex_conc(A_total_M, B_total_M, kd_lo)
        ab_lo = _complex_conc(A_total_M, B_total_M, kd_hi)
        complex_ci = (ab_lo, ab_hi)

    return ComplexationResult(
        Kd_M=Kd,
        A_total_M=A_total_M,
        B_total_M=B_total_M,
        complex_M=ab,
        free_A_M=free_A,
        free_B_M=free_B,
        complexed_fraction_A=frac_A,
        complex_ci95=complex_ci,
        trust=getattr(prior_AB, "trust", "uncalibrated"),
        accuracy_known=getattr(prior_AB, "accuracy_known", False),
    )
