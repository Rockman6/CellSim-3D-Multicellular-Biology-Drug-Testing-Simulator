"""Resistance evolution — why treatment works, then stops working.

The fate layer says which cells die under a given exposure. But a tumour
is not one clone: a small resistant subpopulation is usually present
before treatment starts (or arises by mutation). Treatment does not
create resistance — it SELECTS for it, by killing the sensitive majority
and leaving the resistant minority to repopulate.

Two competing subpopulations, each growing at its own net rate:

    dN_s/dt = k_s·N_s        dN_r/dt = k_r·N_r
    N_i(t)  = N_i(0)·exp(k_i·t)

where k_s, k_r come from `fate.fate_from_occupancy` at each clone's OWN
effective occupancy. Resistance is modelled mechanistically as a shift in
the resistant clone's occupancy at the same drug concentration — a
mutated target binds the drug more weakly (K_d up), or an upregulated
pump lowers the intracellular concentration. Both routes are just "the
resistant clone sees less effective engagement".

The selection coefficient s = k_r − k_s drives everything:

    resistant fraction   f(t) = f₀·e^{st} / (1 − f₀ + f₀·e^{st})

which is a logistic sweep. RELAPSE is when the total population returns
to its starting size.

The clinically important asymmetry: resistance usually carries a FITNESS
COST (the mutated target works less well, the pump costs ATP), so without
drug s < 0 and the resistant clone is OUT-competed. That is the rationale
for drug holidays and adaptive therapy — and it falls out of this model
rather than being asserted.

Non-AI: closed-form two-clone competition. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

from src.cell.fate import CellFateParams, fate_from_occupancy


@dataclass
class ResistantClone:
    """A resistant subpopulation defined by how much engagement it escapes.

    `occupancy_scale` ∈ (0,1] multiplies the occupancy the sensitive clone
    would experience at the same drug concentration: 0.1 means the
    resistant cells achieve only a tenth of the target engagement (weaker
    binding, or a pump keeping the drug out).

    `fitness_cost` ∈ [0,1) reduces the resistant clone's DRUG-FREE
    proliferation rate — the price of the resistance mechanism.
    """

    occupancy_scale: float = 0.1
    fitness_cost: float = 0.1
    initial_fraction: float = 1e-4
    mechanism: str = "target mutation"
    source: str = "input"
    trust: str = "uncalibrated"

    def __post_init__(self):
        if not (0.0 < self.occupancy_scale <= 1.0):
            raise ValueError("occupancy_scale must be in (0,1]")
        if not (0.0 <= self.fitness_cost < 1.0):
            raise ValueError("fitness_cost must be in [0,1)")
        if not (0.0 <= self.initial_fraction < 1.0):
            raise ValueError("initial_fraction must be in [0,1)")


@dataclass
class ResistanceOutcome:
    """Two-clone competition under a sustained exposure."""

    occupancy: float
    k_sensitive_per_s: float
    k_resistant_per_s: float
    selection_coefficient_per_s: float     # s = k_r − k_s
    initial_resistant_fraction: float
    fate_sensitive: str
    fate_resistant: str
    trust: str = "uncalibrated"

    def total_population_at(self, t_s: float, N0: float = 1.0) -> float:
        f0 = self.initial_resistant_fraction
        return N0 * ((1.0 - f0) * math.exp(self.k_sensitive_per_s * t_s)
                     + f0 * math.exp(self.k_resistant_per_s * t_s))

    def resistant_fraction_at(self, t_s: float) -> float:
        """Logistic sweep of the resistant clone."""
        f0 = self.initial_resistant_fraction
        if f0 <= 0.0:
            return 0.0
        num = f0 * math.exp(self.selection_coefficient_per_s * t_s)
        return num / (1.0 - f0 + num)

    def time_to_relapse_s(self, horizon_s: float = 365 * 24 * 3600.0
                          ) -> Optional[float]:
        """When the total population regrows to its starting size.

        None if it never does within `horizon_s` (durable response) or if
        it never shrinks in the first place (no response).
        """
        if self.total_population_at(horizon_s) < 1.0:
            return None                       # still below baseline
        # Must dip below 1 first, else there was no response to relapse from.
        n_probe = 400
        dipped = any(self.total_population_at(horizon_s * i / n_probe) < 0.999
                     for i in range(1, n_probe + 1))
        if not dipped:
            return None
        lo, hi = 0.0, horizon_s
        # Find the last crossing back up through N0.
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.total_population_at(mid) < 1.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def summary(self) -> str:
        s_day = self.selection_coefficient_per_s * 86400.0
        t_rel = self.time_to_relapse_s()
        rel = (f"relapse at {t_rel / 86400.0:.0f} d" if t_rel is not None
               else "no relapse in 1 y")
        return (f"occ {100*self.occupancy:.0f}%: sensitive {self.fate_sensitive}, "
                f"resistant {self.fate_resistant}; s={s_day:+.3f}/day; {rel} "
                f"[{self.trust}]")


def resistance_outcome(occupancy: float,
                       clone: Optional[ResistantClone] = None,
                       params: Optional[CellFateParams] = None
                       ) -> ResistanceOutcome:
    """Compete a sensitive and a resistant clone at a sustained occupancy."""
    c = clone or ResistantClone()
    p = params or CellFateParams()

    f_sens = fate_from_occupancy(occupancy, p)

    # The resistant clone: less engagement, and a fitness cost on its
    # drug-free proliferation.
    p_res = CellFateParams(
        k_prolif_per_s=p.k_prolif_per_s * (1.0 - c.fitness_cost),
        k_death_basal_per_s=p.k_death_basal_per_s,
        k_maxkill_per_s=p.k_maxkill_per_s,
        cytostatic_fraction=p.cytostatic_fraction,
        transduction_tau=p.transduction_tau,
        source=f"{p.source}+resistant:{c.mechanism}",
        trust=p.trust)
    f_res = fate_from_occupancy(occupancy * c.occupancy_scale, p_res)

    return ResistanceOutcome(
        occupancy=occupancy,
        k_sensitive_per_s=f_sens.k_net_per_s,
        k_resistant_per_s=f_res.k_net_per_s,
        selection_coefficient_per_s=f_res.k_net_per_s - f_sens.k_net_per_s,
        initial_resistant_fraction=c.initial_fraction,
        fate_sensitive=f_sens.fate,
        fate_resistant=f_res.fate,
        trust=c.trust,
    )
