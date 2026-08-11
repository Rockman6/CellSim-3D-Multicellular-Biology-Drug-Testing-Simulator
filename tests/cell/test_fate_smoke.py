"""Cell fate — proliferate, arrest, or die.

Rigor anchors:
  * untreated cells recover the cited 24 h doubling time;
  * the CLOSED-FORM critical occupancy θ* agrees with the numerically
    located zero-crossing of k_net (independent derivation vs search);
  * population growth matches N₀·exp(k_net·t) exactly;
  * receptor reserve (τ) makes low occupancy sufficient — the reason
    "% engagement" ≠ "% effect";
  * a drug whose maximal effect cannot outrun proliferation reports
    θ* = None instead of a fake threshold.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    CellFateParams,
    effect_from_occupancy,
    fate_from_occupancy,
    critical_occupancy,
    tissue_fate_profile,
    penetration_profile_first_order,
)

_DAY = 24.0 * 3600.0


def test_untreated_cells_double_on_schedule():
    f = fate_from_occupancy(0.0)
    assert f.fate == "proliferating"
    assert abs(f.doubling_time_s - _DAY) / _DAY < 1e-9
    # Population doubles in exactly one doubling time.
    assert abs(f.population_at(_DAY, N0=1.0) - 2.0) < 1e-9


def test_full_occupancy_kills():
    f = fate_from_occupancy(1.0)
    assert f.fate == "dying"
    assert f.k_net_per_s < 0
    assert f.half_life_s is not None
    assert f.population_at(_DAY) < 1.0


def test_closed_form_critical_occupancy_matches_numeric_zero_crossing():
    """THE anchor: the analytic θ* equals the numerically bisected root."""
    for params in (
        CellFateParams(),
        CellFateParams(transduction_tau=10.0),
        CellFateParams(cytostatic_fraction=0.0),
        CellFateParams(k_death_basal_per_s=1e-6),
    ):
        theta_star = critical_occupancy(params)
        assert theta_star is not None
        # k_net must vanish there.
        assert abs(fate_from_occupancy(theta_star, params).k_net_per_s) < 1e-12
        # Bisect independently and compare.
        lo, hi = 0.0, 1.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if fate_from_occupancy(mid, params).k_net_per_s > 0:
                lo = mid
            else:
                hi = mid
        assert abs(theta_star - 0.5 * (lo + hi)) < 1e-9


def test_receptor_reserve_lowers_the_bar():
    """Larger τ (spare receptors) → less occupancy needed to stop growth."""
    low = critical_occupancy(CellFateParams(transduction_tau=1.0))
    high = critical_occupancy(CellFateParams(transduction_tau=20.0))
    assert high < low
    # And effect always exceeds occupancy when τ > 1.
    assert effect_from_occupancy(0.1, tau=20.0) > 0.1


def test_drug_that_cannot_win_reports_none():
    """If max effect can't outrun proliferation, θ* is None, not a number."""
    weak = CellFateParams(k_maxkill_per_s=0.0, cytostatic_fraction=0.5)
    assert critical_occupancy(weak) is None
    # Even at full occupancy the cells still grow.
    assert fate_from_occupancy(1.0, weak).fate == "proliferating"


def test_arrest_band_between_growth_and_death():
    """Just at θ* the population is arrested, not growing or dying."""
    p = CellFateParams()
    theta = critical_occupancy(p)
    assert fate_from_occupancy(theta, p).fate == "arrested"


def test_tissue_fate_shows_survivors_in_the_depths():
    """The payoff: cells die near the vessel and survive deeper in."""
    prior = binding_to_hill(-9.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    prof = penetration_profile_first_order(1e-5, thickness_um=150.0,
                                           k_per_s=1e-1)
    tf = tissue_fate_profile(prof, prior)
    assert tf.fates[0] == "dying"            # vessel-side cells die
    assert tf.fates[-1] == "proliferating"   # deep cells keep growing
    assert 0.0 < tf.killed_fraction < 1.0
    assert tf.critical_occupancy is not None
