"""Two-compartment passive membrane permeation (extracellular ↔ cytoplasm).

Answers "the drug is outside the cell — how much gets *inside*, and how
fast?", so target occupancy can be computed at the LOCAL intracellular
concentration instead of a single bulk number.

Physics (passive diffusion, inhomogeneous-solubility-diffusion reduced to
a single permeability coefficient P):

    dC_in/dt = (P · A / V) · (C_out − C_in)

with a large, buffered extracellular reservoir (C_out fixed). Closed form:

    C_in(t) = C_out · (1 − exp(−t/τ)),     τ = V / (P · A)

For a sphere of radius r this collapses to the clean τ = r / (3P): the
permeability sets the *timescale*, and passive diffusion with no trapping
drives C_in → C_out at steady state.

Provenance discipline (MISSION.md — never a magic number): the
permeability P is an INPUT carrying its own source + trust. It would
ultimately come from a Martini bilayer PMF + position-dependent diffusion
profile (Layer 1.5); that calculation is not wired yet, so a P not from
our own physics is flagged `uncalibrated`. Cell geometry uses a cited
typical mammalian cell; override for a specific cell type.

Honest limitations of THIS version (documented, not hidden):
  * passive permeation only — no active transport / efflux pumps;
  * neutral-species model — no pH / ion trapping (weak-base lysosomal
    accumulation needs pKa + compartment pH; deferred);
  * no intracellular binding sink buffering the free concentration;
  * P uncertainty is not yet propagated into the occupancy CI — only the
    K_d CI is. The permeability `trust` flag rides through so the readout
    is not called decision-grade when P is uncalibrated.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from src.cell.occupancy import occupancy_from_prior, OccupancyResult

# Typical mammalian cell radius ~10 µm (Alberts, Molecular Biology of the
# Cell; BioNumbers BNID 100434). Used only when the caller doesn't supply
# a specific geometry.
_DEFAULT_CELL_RADIUS_UM = 10.0


@dataclass
class Permeability:
    """Passive membrane permeability coefficient with provenance.

    `trust`: 'calibrated' only when P came from our own physics (a Martini
    PMF/diffusion calculation). Literature or placeholder values are
    'uncalibrated' — used, but flagged so a readout built on them is not
    called decision-grade.
    """

    P_cm_per_s: float
    source: str = "placeholder"          # 'martini' | 'literature' | 'placeholder'
    trust: str = "uncalibrated"           # 'calibrated' | 'uncalibrated'
    note: str = ""

    def __post_init__(self):
        if self.P_cm_per_s <= 0:
            raise ValueError("P_cm_per_s must be > 0")


@dataclass
class CellGeometry:
    """Compartment geometry in CGS (cm, cm², cm³) so τ = V/(P·A) is in s."""

    volume_cm3: float
    membrane_area_cm2: float
    label: str = "spherical cell"

    def __post_init__(self):
        if self.volume_cm3 <= 0 or self.membrane_area_cm2 <= 0:
            raise ValueError("volume and area must be > 0")


def spherical_cell_geometry(radius_um: float = _DEFAULT_CELL_RADIUS_UM
                            ) -> CellGeometry:
    """Geometry of a spherical cell of the given radius (µm → CGS)."""
    if radius_um <= 0:
        raise ValueError("radius_um must be > 0")
    r_cm = radius_um * 1e-4                     # µm → cm
    v = (4.0 / 3.0) * math.pi * r_cm ** 3
    a = 4.0 * math.pi * r_cm ** 2
    return CellGeometry(volume_cm3=v, membrane_area_cm2=a,
                        label=f"sphere r={radius_um:g}µm")


def equilibration_tau_s(perm: Permeability, geom: CellGeometry) -> float:
    """Time constant τ = V / (P · A) for C_in to approach C_out (seconds)."""
    return geom.volume_cm3 / (perm.P_cm_per_s * geom.membrane_area_cm2)


def intracellular_concentration(C_out_M: float, t_s: float,
                                perm: Permeability,
                                geom: CellGeometry) -> float:
    """Intracellular concentration at time t under passive permeation.

    C_in(t) = C_out · (1 − exp(−t/τ)). t_s = math.inf gives the steady
    state (→ C_out for this passive, non-trapping model).
    """
    if C_out_M < 0 or t_s < 0:
        raise ValueError("C_out_M and t_s must be ≥ 0")
    tau = equilibration_tau_s(perm, geom)
    if math.isinf(t_s):
        frac = 1.0
    else:
        frac = 1.0 - math.exp(-t_s / tau)
    return C_out_M * frac


@dataclass
class CellOccupancyResult:
    """Target occupancy computed at the LOCAL intracellular concentration.

    Distinguishes the naive bulk readout (occupancy at C_out) from the
    membrane-aware one (occupancy at C_in), and carries BOTH the K_d
    trust and the permeability trust.
    """

    occupancy: OccupancyResult            # computed at C_in
    C_out_M: float
    C_in_M: float
    t_s: float
    tau_s: float
    fraction_equilibrated: float
    permeability_trust: str
    naive_occupancy_theta: float          # occupancy if you ignored the membrane

    @property
    def decision_grade(self) -> bool:
        """Decision-grade only if the K_d AND the permeability are both
        calibrated — a membrane-aware readout is only as trustworthy as
        its weakest physics input."""
        return (self.occupancy.decision_grade
                and self.permeability_trust == "calibrated")

    def summary(self) -> str:
        eq = 100.0 * self.fraction_equilibrated
        grade = "decision-grade" if self.decision_grade else "NOT decision-grade"
        return (
            f"at t={self.t_s:.0f}s: C_in={self.C_in_M:.2e} M "
            f"({eq:.0f}% equilibrated, τ={self.tau_s:.0f}s) → "
            f"target occupancy {100*self.occupancy.theta:.1f}% "
            f"(bulk-C_out would read {100*self.naive_occupancy_theta:.1f}%) "
            f"— {grade}")


def occupancy_in_cell(prior, C_out_M: float, perm: Permeability,
                      geom: Optional[CellGeometry] = None,
                      t_s: float = math.inf) -> CellOccupancyResult:
    """Target occupancy at the intracellular concentration reached by time t.

    Ties the layers together: a molecular ΔG (→ K_d prior) plus a membrane
    permeability give the fraction of an *intracellular* target bound when
    the drug is applied *outside* the cell. Contrast `naive_occupancy_theta`
    (occupancy at the bulk C_out) with the membrane-aware value to see how
    much the barrier matters at time t.
    """
    if geom is None:
        geom = spherical_cell_geometry()
    tau = equilibration_tau_s(perm, geom)
    C_in = intracellular_concentration(C_out_M, t_s, perm, geom)
    frac_eq = 1.0 if math.isinf(t_s) else (1.0 - math.exp(-t_s / tau))

    occ_in = occupancy_from_prior(prior, C_in)
    occ_bulk = occupancy_from_prior(prior, C_out_M)

    return CellOccupancyResult(
        occupancy=occ_in,
        C_out_M=C_out_M,
        C_in_M=C_in,
        t_s=t_s,
        tau_s=tau,
        fraction_equilibrated=frac_eq,
        permeability_trust=perm.trust,
        naive_occupancy_theta=occ_bulk.theta,
    )
