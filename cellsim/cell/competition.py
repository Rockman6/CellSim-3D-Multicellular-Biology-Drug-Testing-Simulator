"""Competitive binding: several drugs contending for ONE target site.

Answers "if I put two drugs in, do they interfere, and how does that
change the occupancy 'ruler'?". When ligands share a single binding site
they compete: each one present lowers the others' occupancy.

Physics (single-site competitive equilibrium, the standard result):

    xᵢ = [Lᵢ] / K_d,ᵢ
    θᵢ = xᵢ / (1 + Σⱼ xⱼ)              (fraction of site bound by ligand i)
    free = 1 / (1 + Σⱼ xⱼ)              (fraction of site unbound)

This is the n_hill = 1 single-site model; it does not describe
cooperativity or two ligands binding simultaneously at distinct
sub-sites. Ligands forming a NEW species by binding each other
(drug–drug complexation A + B ⇌ AB) is a different edge — a bimolecular
reaction, not competition — and is not modelled here (see
docs/campaign1_closeout.md Gap C).

Uncertainty: θᵢ is monotonically DECREASING in K_d,ᵢ and monotonically
INCREASING in every other K_d,ⱼ. So a rigorous 95 % interval on θᵢ comes
from interval arithmetic on the priors' K_d CIs — take K_d,ᵢ at one CI
end and all competitors at the opposite end. This is an exact bound for
the interval, not a linearisation. It is reported only when every
competing prior carries a K_d CI.

Trust rides through: the competitive readout is decision-grade only if
every participating K_d is trustworthy in absolute terms.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple


@dataclass
class LigandOccupancy:
    """One ligand's share of the shared site."""

    label: str
    theta: float
    Kd_M: float
    conc_M: float
    theta_ci95: Optional[Tuple[float, float]] = None
    trust: str = "uncalibrated"
    accuracy_known: bool = False


@dataclass
class CompetitiveResult:
    ligands: List[LigandOccupancy]
    free_fraction: float
    free_fraction_ci95: Optional[Tuple[float, float]] = None

    @property
    def decision_grade(self) -> bool:
        """Only if EVERY competing K_d is trustworthy in absolute terms."""
        return all(l.accuracy_known and l.trust == "trustworthy_absolute"
                   for l in self.ligands) and bool(self.ligands)

    def summary(self) -> str:
        parts = [f"{l.label} {100*l.theta:.1f}%" for l in self.ligands]
        grade = "decision-grade" if self.decision_grade else "NOT decision-grade"
        return (f"shared site: {', '.join(parts)}, "
                f"free {100*self.free_fraction:.1f}% — {grade}")


def _kd_and_ci(prior):
    params = getattr(prior, "parameters", {})
    Kd = float(params["Kd_M"])
    n_hill = float(params.get("n_hill", 1.0))
    if abs(n_hill - 1.0) > 1e-9:
        raise ValueError(
            "competitive_occupancy is a single-site (n_hill=1) model; "
            f"got n_hill={n_hill}. Cooperative competition is not supported.")
    ci = (getattr(prior, "parameter_ci95", None) or {}).get("Kd_M")
    return Kd, ci


def competitive_occupancy(
    entries: Sequence[Tuple[str, object, float]],
) -> CompetitiveResult:
    """Occupancy of a single shared site by several competing ligands.

    Args:
        entries: sequence of (label, RateLawPrior, concentration_M). Each
            prior must be a Hill prior with n_hill = 1.

    Returns: per-ligand occupancy (with rigorous CI when all priors carry
    a K_d CI) plus the free-site fraction.
    """
    if not entries:
        raise ValueError("need at least one (label, prior, conc) entry")

    labels, priors, concs = [], [], []
    kds, cis = [], []
    for label, prior, conc_M in entries:
        if conc_M < 0:
            raise ValueError(f"{label}: concentration must be ≥ 0")
        Kd, ci = _kd_and_ci(prior)
        labels.append(label)
        priors.append(prior)
        concs.append(float(conc_M))
        kds.append(Kd)
        cis.append(ci)

    n = len(entries)
    x = [concs[i] / kds[i] for i in range(n)]      # [L]/Kd per ligand
    denom = 1.0 + sum(x)
    thetas = [x[i] / denom for i in range(n)]
    free = 1.0 / denom

    all_have_ci = all(ci is not None for ci in cis)

    ligands: List[LigandOccupancy] = []
    for i in range(n):
        theta_ci = None
        if all_have_ci:
            # θᵢ max: Kd_i at its LOW end (bigger x_i), competitors at
            # their HIGH end (smaller x_j). θᵢ min: the reverse.
            def _theta_i(kd_i, other_end):
                xi = concs[i] / kd_i
                s = xi
                for j in range(n):
                    if j == i:
                        continue
                    kd_j = cis[j][other_end]
                    s_add = concs[j] / kd_j
                    s += s_add
                return xi / (1.0 + s)
            theta_hi = _theta_i(cis[i][0], 1)   # Kd_i low, others high
            theta_lo = _theta_i(cis[i][1], 0)   # Kd_i high, others low
            theta_ci = (theta_lo, theta_hi)
        ligands.append(LigandOccupancy(
            label=labels[i], theta=thetas[i], Kd_M=kds[i], conc_M=concs[i],
            theta_ci95=theta_ci,
            trust=getattr(priors[i], "trust", "uncalibrated"),
            accuracy_known=getattr(priors[i], "accuracy_known", False),
        ))

    free_ci = None
    if all_have_ci:
        # free is DECREASING in every x, i.e. INCREASING in every Kd.
        free_hi = 1.0 / (1.0 + sum(concs[j] / cis[j][1] for j in range(n)))
        free_lo = 1.0 / (1.0 + sum(concs[j] / cis[j][0] for j in range(n)))
        free_ci = (free_lo, free_hi)

    return CompetitiveResult(ligands=ligands, free_fraction=free,
                             free_fraction_ci95=free_ci)
