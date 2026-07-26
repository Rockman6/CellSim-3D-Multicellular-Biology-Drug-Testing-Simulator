"""Competitive binding: two drugs contend for one site — does the ruler move?

Pins the single-site competitive equilibrium, the conservation identity
(all occupancies + free = 1), the direction of interference (adding a
competitor lowers the other's occupancy), rigorous interval-arithmetic
CI propagation over the K_d CIs, and trust ride-through.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import competitive_occupancy  # noqa: E402


def _prior(dG, rec="benchmarks/dock/3ptb.pdb", sigma=0.2):
    return binding_to_hill(dG, uncertainty_kcalmol=sigma, receptor=rec,
                           method="vina")


def test_occupancies_and_free_sum_to_one():
    a = _prior(-9.0)
    b = _prior(-7.0)
    res = competitive_occupancy([("A", a, 1e-7), ("B", b, 1e-7)])
    total = sum(l.theta for l in res.ligands) + res.free_fraction
    assert abs(total - 1.0) < 1e-12


def test_adding_a_competitor_lowers_occupancy():
    """A's occupancy must drop when a competitor B is added."""
    a = _prior(-9.0)
    b = _prior(-9.0)  # equally strong competitor
    alone = competitive_occupancy([("A", a, 1e-7)])
    with_b = competitive_occupancy([("A", a, 1e-7), ("B", b, 1e-7)])
    theta_alone = alone.ligands[0].theta
    theta_with_b = with_b.ligands[0].theta
    assert theta_with_b < theta_alone
    # Two equal competitors at equal conc → equal shares.
    assert abs(with_b.ligands[0].theta - with_b.ligands[1].theta) < 1e-12


def test_stronger_binder_wins_the_site():
    strong = _prior(-11.0)   # tighter K_d
    weak = _prior(-6.0)
    res = competitive_occupancy([("strong", strong, 1e-7),
                                 ("weak", weak, 1e-7)])
    theta = {l.label: l.theta for l in res.ligands}
    assert theta["strong"] > theta["weak"]


def test_ci_brackets_point_estimate_and_widens_with_uncertainty():
    a = _prior(-9.0, sigma=0.2)
    b = _prior(-7.0, sigma=0.2)
    res = competitive_occupancy([("A", a, 1e-7), ("B", b, 1e-7)])
    for l in res.ligands:
        assert l.theta_ci95 is not None
        lo, hi = l.theta_ci95
        assert lo <= l.theta <= hi
    # free-fraction CI present and brackets the point value.
    flo, fhi = res.free_fraction_ci95
    assert flo <= res.free_fraction <= fhi


def test_trust_rides_through():
    """Untrustworthy competitor → whole readout not decision-grade."""
    good = _prior(-9.0, rec="benchmarks/dock/3ptb.pdb")   # trustworthy
    bad = _prior(-9.0, rec="benchmarks/dock/1stp.pdb")     # do_not_trust
    res = competitive_occupancy([("good", good, 1e-7), ("bad", bad, 1e-7)])
    assert not res.decision_grade

    res2 = competitive_occupancy([("g1", good, 1e-7),
                                  ("g2", _prior(-8.0), 1e-7)])
    assert res2.decision_grade


def test_rejects_cooperative_prior():
    coop = binding_to_hill(-9.0, n_hill=2.0, uncertainty_kcalmol=0.2,
                           receptor="benchmarks/dock/3ptb.pdb", method="vina")
    try:
        competitive_occupancy([("coop", coop, 1e-7)])
        assert False, "should reject n_hill != 1"
    except ValueError:
        pass


def test_rejects_empty():
    try:
        competitive_occupancy([])
        assert False
    except ValueError:
        pass
