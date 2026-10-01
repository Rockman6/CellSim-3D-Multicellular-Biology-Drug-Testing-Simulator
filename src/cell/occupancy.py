"""Single-site target occupancy at a given drug concentration.

The simplest cell-level readout: given a drug concentration [L] around a
target and a binding K_d (as a `RateLawPrior` from `src/bridge`), what
fraction of the target is bound?

    θ([L]) = [L]^n / (K_d^n + [L]^n)          (Hill occupancy)

θ is what a downstream signalling model actually needs — an ion channel
that is 90 % blocked behaves very differently from one 10 % blocked. This
module computes θ AND propagates the prior's uncertainty and trust so the
occupancy is never more confident than the K_d behind it.

Uncertainty propagation is exact, not linearised: θ is a monotonically
DECREASING function of K_d, so the occupancy CI is obtained by evaluating
θ at the K_d CI endpoints (a wider K_d CI → a wider θ CI), which is
correct for any n > 0.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple


def hill_occupancy(ligand_conc_M: float, Kd_M: float,
                   n_hill: float = 1.0) -> float:
    """Fraction of target bound at [L] for dissociation constant K_d.

    θ = L^n / (K_d^n + L^n). Returns a value in [0, 1]. At [L] = K_d and
    n = 1 this is exactly 0.5, the definition of K_d.
    """
    if ligand_conc_M < 0 or Kd_M <= 0:
        raise ValueError("ligand_conc_M must be ≥ 0 and Kd_M must be > 0")
    if ligand_conc_M == 0:
        return 0.0
    ln = ligand_conc_M ** n_hill
    kn = Kd_M ** n_hill
    return ln / (kn + ln)


@dataclass
class OccupancyResult:
    """Cell-level occupancy readout with the prior's uncertainty + trust.

    `theta_ci95` is populated only when the source prior carried a K_d CI.
    `trust` and `accuracy_known` are inherited from the prior so a
    consumer can refuse to act on an untrustworthy readout.
    """

    theta: float                                   # fraction bound [0,1]
    ligand_conc_M: float
    Kd_M: float
    n_hill: float
    theta_ci95: Optional[Tuple[float, float]] = None
    trust: str = "uncalibrated"
    accuracy_known: bool = False
    method: str = ""
    notes: str = ""

    @property
    def decision_grade(self) -> bool:
        """True only if the underlying K_d is trustworthy in absolute terms.

        An occupancy computed from a `rank_order_only`, `do_not_trust_
        absolute`, or `uncalibrated` K_d is informative for ordering but
        must not be read as a calibrated absolute fraction.
        """
        return self.accuracy_known and self.trust == "trustworthy_absolute"

    def summary(self) -> str:
        pct = 100.0 * self.theta
        if self.theta_ci95 is not None:
            lo, hi = self.theta_ci95
            ci = f" (95% CI {100*lo:.1f}–{100*hi:.1f} %)"
        else:
            ci = " (no CI — uncalibrated)"
        grade = "decision-grade" if self.decision_grade else \
            f"NOT decision-grade [{self.trust}]"
        return (f"occupancy {pct:.1f} %{ci} at "
                f"[L]={self.ligand_conc_M:.2e} M — {grade}")


def occupancy_from_prior(prior, ligand_conc_M: float) -> OccupancyResult:
    """Compute target occupancy from a `src/bridge` Hill `RateLawPrior`.

    Carries the prior's K_d CI through to an occupancy CI (exact, via the
    monotonic K_d→θ relationship) and inherits its `trust` / accuracy
    flags so the readout is never more confident than the K_d behind it.
    """
    if getattr(prior, "type", None) != "hill":
        raise ValueError(
            f"occupancy_from_prior needs a Hill prior, got '{getattr(prior, 'type', None)}'")
    params = prior.parameters
    Kd_M = float(params["Kd_M"])
    n_hill = float(params.get("n_hill", 1.0))

    theta = hill_occupancy(ligand_conc_M, Kd_M, n_hill)

    theta_ci95 = None
    kd_ci = (prior.parameter_ci95 or {}).get("Kd_M")
    if kd_ci is not None:
        kd_lo, kd_hi = kd_ci
        # θ decreases in K_d: the LOW K_d gives the HIGH occupancy.
        theta_hi = hill_occupancy(ligand_conc_M, kd_lo, n_hill)
        theta_lo = hill_occupancy(ligand_conc_M, kd_hi, n_hill)
        theta_ci95 = (theta_lo, theta_hi)

    return OccupancyResult(
        theta=theta,
        ligand_conc_M=ligand_conc_M,
        Kd_M=Kd_M,
        n_hill=n_hill,
        theta_ci95=theta_ci95,
        trust=getattr(prior, "trust", "uncalibrated"),
        accuracy_known=getattr(prior, "accuracy_known", False),
        method=getattr(prior, "method", ""),
        notes=getattr(prior, "notes", ""),
    )
