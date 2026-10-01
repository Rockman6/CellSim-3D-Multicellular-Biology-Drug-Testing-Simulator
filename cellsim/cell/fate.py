"""Cell fate — does the cell proliferate, arrest, or die?

Everything upstream computes how much drug reaches a target and how much
of that target is bound. This module closes the loop to the thing the
project is actually about: what the CELL does about it.

Two steps, both standard and cited.

1. Occupancy → effect (Black–Leff operational model, 1983, Proc R Soc
   Lond B 220:141):

       E(θ) = τ·θ / (1 + τ·θ)

   τ is the transduction coefficient ("receptor reserve"). τ ≫ 1 means a
   small fraction of targets bound already produces a near-maximal
   response — spare receptors, common for enzymes and GPCRs. τ = 1 means
   effect tracks occupancy weakly. This is why "% target engagement" and
   "% effect" are not the same number, and the model says so explicitly
   rather than assuming E = θ.

2. Effect → net growth rate:

       k_net = k_prolif·(1 − c·E)  −  k_death_basal  −  k_maxkill·E

   where c ∈ [0,1] is the cytostatic fraction (how much of the drug's
   action is stopping division rather than killing) and k_maxkill is the
   maximal cytotoxic rate. Fate is read off the sign of k_net:

       k_net > 0  proliferating       k_net ≈ 0  arrested
       k_net < 0  dying

   Population follows N(t) = N₀·exp(k_net·t) for a constant exposure.

The clinically meaningful derived quantity is the CRITICAL OCCUPANCY θ*
at which k_net crosses zero — the target engagement a drug must sustain
to stop the population growing. It has a closed form (see
`critical_occupancy`), and if the required effect exceeds 1 the drug
CANNOT kill these cells at any dose — a property of the pharmacology, not
of the concentration.

Cited defaults: mammalian doubling time ~24 h (Alberts, MBoC; BioNumbers
BNID 100685) → k_prolif = ln2/24 h. Basal apoptosis is slow compared with
division in a growing culture.

Non-AI: closed-form pharmacodynamics. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

# Mammalian cell doubling time ~24 h (Alberts, Molecular Biology of the
# Cell; BioNumbers BNID 100685).
_DEFAULT_DOUBLING_TIME_S = 24.0 * 3600.0
_LN2 = math.log(2.0)

# Below this |k_net| (relative to the basal proliferation rate) the
# population is neither meaningfully growing nor shrinking: arrested.
_ARREST_BAND_FRACTION = 0.05


@dataclass
class CellFateParams:
    """Growth / death pharmacodynamics of a cell type.

    All rates in s⁻¹. Provenance: these are cell-line properties, inputs
    carrying source + trust (measured doubling time, drug Emax from a
    viability assay). Defaults describe a generic dividing mammalian cell.
    """

    k_prolif_per_s: float = _LN2 / _DEFAULT_DOUBLING_TIME_S
    k_death_basal_per_s: float = 0.0
    k_maxkill_per_s: float = 2.0 * _LN2 / _DEFAULT_DOUBLING_TIME_S
    cytostatic_fraction: float = 1.0     # c ∈ [0,1]
    transduction_tau: float = 1.0        # τ — receptor reserve
    source: str = "default:generic-mammalian"
    trust: str = "uncalibrated"

    def __post_init__(self):
        if self.k_prolif_per_s < 0 or self.k_death_basal_per_s < 0 \
                or self.k_maxkill_per_s < 0:
            raise ValueError("rates must be ≥ 0")
        if not (0.0 <= self.cytostatic_fraction <= 1.0):
            raise ValueError("cytostatic_fraction must be in [0,1]")
        if self.transduction_tau <= 0:
            raise ValueError("transduction_tau must be > 0")

    @property
    def doubling_time_s(self) -> Optional[float]:
        """Untreated doubling time, or None if the cells don't grow."""
        k = self.k_prolif_per_s - self.k_death_basal_per_s
        return (_LN2 / k) if k > 0 else None


def effect_from_occupancy(occupancy: float, tau: float = 1.0) -> float:
    """Black–Leff transduction: fraction of maximal effect at occupancy θ."""
    if not (0.0 <= occupancy <= 1.0):
        raise ValueError("occupancy must be in [0,1]")
    if tau <= 0:
        raise ValueError("tau must be > 0")
    return tau * occupancy / (1.0 + tau * occupancy)


@dataclass
class CellFate:
    """The cell's fate under a given sustained target occupancy."""

    occupancy: float
    effect: float
    k_net_per_s: float
    fate: str                       # 'proliferating' | 'arrested' | 'dying'
    doubling_time_s: Optional[float]      # None unless proliferating
    half_life_s: Optional[float]          # None unless dying
    trust: str = "uncalibrated"

    def population_at(self, t_s: float, N0: float = 1.0) -> float:
        """N(t) = N₀·exp(k_net·t) for sustained exposure."""
        if t_s < 0 or N0 < 0:
            raise ValueError("t_s and N0 must be ≥ 0")
        return N0 * math.exp(self.k_net_per_s * t_s)

    def summary(self) -> str:
        if self.fate == "proliferating":
            detail = f"doubling {self.doubling_time_s / 3600:.1f} h"
        elif self.fate == "dying":
            detail = f"half-life {self.half_life_s / 3600:.1f} h"
        else:
            detail = "growth arrested"
        return (f"occupancy {100*self.occupancy:.1f}% → effect "
                f"{100*self.effect:.1f}% → {self.fate.upper()} ({detail}) "
                f"[{self.trust}]")


def fate_from_occupancy(occupancy: float,
                        params: Optional[CellFateParams] = None,
                        trust: str = "uncalibrated") -> CellFate:
    """Cell fate under a sustained target occupancy."""
    p = params or CellFateParams()
    E = effect_from_occupancy(occupancy, p.transduction_tau)
    k_net = (p.k_prolif_per_s * (1.0 - p.cytostatic_fraction * E)
             - p.k_death_basal_per_s - p.k_maxkill_per_s * E)

    band = _ARREST_BAND_FRACTION * p.k_prolif_per_s
    if k_net > band:
        fate, dt, hl = "proliferating", _LN2 / k_net, None
    elif k_net < -band:
        fate, dt, hl = "dying", None, _LN2 / abs(k_net)
    else:
        fate, dt, hl = "arrested", None, None

    return CellFate(occupancy=occupancy, effect=E, k_net_per_s=k_net,
                    fate=fate, doubling_time_s=dt, half_life_s=hl,
                    trust=trust)


def critical_occupancy(params: Optional[CellFateParams] = None
                       ) -> Optional[float]:
    """Occupancy θ* at which k_net = 0 — the bar a drug must clear.

    Closed form: setting k_net = 0 gives the required effect

        E* = (k_prolif − k_death_basal) / (k_prolif·c + k_maxkill)

    and inverting the Black–Leff relation gives θ* = E*/(τ·(1 − E*)).
    Returns None if E* ≥ 1 (the drug can never stop these cells, at any
    occupancy) or if the cells are not growing to begin with.
    """
    p = params or CellFateParams()
    numer = p.k_prolif_per_s - p.k_death_basal_per_s
    if numer <= 0:
        return 0.0                        # already not growing
    denom = p.k_prolif_per_s * p.cytostatic_fraction + p.k_maxkill_per_s
    if denom <= 0:
        return None
    E_star = numer / denom
    if E_star >= 1.0:
        return None                        # unreachable — drug can't win
    theta = E_star / (p.transduction_tau * (1.0 - E_star))
    return theta if theta <= 1.0 else None


@dataclass
class TissueFateProfile:
    """Cell fate as a function of depth into tissue."""

    x_um: List[float]
    occupancy: List[float]
    k_net_per_s: List[float]
    fates: List[str]
    critical_occupancy: Optional[float]
    killed_fraction: float               # fraction of depth that is dying
    trust: str = "uncalibrated"

    def summary(self) -> str:
        crit = (f"θ*={100*self.critical_occupancy:.1f}%"
                if self.critical_occupancy is not None else "θ* unreachable")
        return (f"{crit}; {100*self.killed_fraction:.0f}% of tissue depth "
                f"dying, {100*(1-self.killed_fraction):.0f}% surviving "
                f"[{self.trust}]")


def tissue_fate_profile(profile, prior,
                        params: Optional[CellFateParams] = None
                        ) -> TissueFateProfile:
    """Fate of the cells at each depth of a `tissue.TissueProfile`.

    The payoff of the whole stack: cells near the vessel see enough drug to
    die while cells deeper in survive and keep proliferating — the spatial
    picture of why a tumour regrows from its interior.
    """
    p = params or CellFateParams()
    occs = profile.occupancy_profile(prior)
    fates = [fate_from_occupancy(o, p) for o in occs]
    dying = sum(1 for f in fates if f.fate == "dying")
    return TissueFateProfile(
        x_um=list(profile.x_um),
        occupancy=occs,
        k_net_per_s=[f.k_net_per_s for f in fates],
        fates=[f.fate for f in fates],
        critical_occupancy=critical_occupancy(p),
        killed_fraction=dying / len(fates) if fates else 0.0,
        trust=profile.trust,
    )
