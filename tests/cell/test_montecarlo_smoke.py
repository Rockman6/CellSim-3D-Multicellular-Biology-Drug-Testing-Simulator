"""Joint Monte-Carlo uncertainty propagation — and its cross-check.

The headline rigor test: for a SINGLE prior feeding a MONOTONE readout,
the MC percentile CI must agree with the exact interval-arithmetic CI the
analytic occupancy module computes. Agreement validates BOTH methods
against each other; the MC then earns trust for the multi-prior /
non-monotone cases where the analytic interval no longer applies.

Also pins determinism (seed), that a wider ΔG σ gives a wider output
distribution, and trust ride-through.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    occupancy_from_prior,
    hill_occupancy,
    monte_carlo_propagate,
)


def _prior(dG, rec="benchmarks/dock/3ptb.pdb", sigma=0.4):
    return binding_to_hill(dG, uncertainty_kcalmol=sigma, receptor=rec,
                           method="vina")


def test_mc_ci_agrees_with_exact_interval_for_single_monotone_readout():
    """MC 95% CI ≈ analytic exact interval for one prior, monotone readout."""
    prior = _prior(-8.0, sigma=0.5)
    L = 1e-7

    analytic = occupancy_from_prior(prior, L)
    a_lo, a_hi = analytic.theta_ci95

    mc = monte_carlo_propagate(
        {"A": prior},
        lambda kd: hill_occupancy(L, kd["A"], 1.0),
        n=20000, seed=1)

    # Point estimate matches.
    assert abs(mc.median - analytic.theta) < 0.02
    # The two 95% CIs agree to a few percent of occupancy. (Analytic uses
    # ±1.96σ endpoints; MC uses 2.5/97.5 percentiles — same target.)
    assert abs(mc.ci95[0] - a_lo) < 0.03
    assert abs(mc.ci95[1] - a_hi) < 0.03


def test_deterministic_given_seed():
    prior = _prior(-9.0)
    fn = lambda kd: hill_occupancy(1e-7, kd["A"], 1.0)  # noqa: E731
    a = monte_carlo_propagate({"A": prior}, fn, n=3000, seed=7)
    b = monte_carlo_propagate({"A": prior}, fn, n=3000, seed=7)
    c = monte_carlo_propagate({"A": prior}, fn, n=3000, seed=8)
    assert a.mean == b.mean and a.ci95 == b.ci95
    assert a.mean != c.mean          # different seed → different draws


def test_wider_sigma_widens_the_output_distribution():
    fn = lambda kd: hill_occupancy(1e-7, kd["A"], 1.0)  # noqa: E731
    tight = monte_carlo_propagate({"A": _prior(-8.0, sigma=0.2)}, fn,
                                  n=8000, seed=3)
    wide = monte_carlo_propagate({"A": _prior(-8.0, sigma=1.5)}, fn,
                                 n=8000, seed=3)
    assert wide.std > tight.std
    assert (wide.ci95[1] - wide.ci95[0]) > (tight.ci95[1] - tight.ci95[0])


def test_joint_two_prior_competition_readout():
    """Two priors feeding a competitive-occupancy readout: the MC handles
    the joint uncertainty that per-prior intervals can't combine."""
    A = _prior(-9.0)
    B = _prior(-8.0)
    L_A = L_B = 1e-7

    def theta_A(kd):
        xA = L_A / kd["A"]
        xB = L_B / kd["B"]
        return xA / (1.0 + xA + xB)

    mc = monte_carlo_propagate({"A": A, "B": B}, theta_A, n=10000, seed=5)
    assert 0.0 <= mc.ci95[0] <= mc.median <= mc.ci95[1] <= 1.0
    assert mc.std > 0.0


def test_trust_is_worst_of_inputs():
    good = _prior(-9.0, rec="benchmarks/dock/3ptb.pdb")   # trustworthy
    bad = _prior(-9.0, rec="benchmarks/dock/1stp.pdb")     # do_not_trust
    fn = lambda kd: hill_occupancy(1e-7, kd["A"], 1.0)  # noqa: E731

    mc_good = monte_carlo_propagate({"A": good}, fn, n=2000, seed=0)
    assert mc_good.decision_grade

    def theta_A(kd):
        xA = 1e-7 / kd["A"]
        xB = 1e-7 / kd["B"]
        return xA / (1.0 + xA + xB)

    mc_mixed = monte_carlo_propagate({"A": good, "B": bad}, theta_A,
                                     n=2000, seed=0)
    assert mc_mixed.trust == "do_not_trust_absolute"   # worst wins
    assert not mc_mixed.decision_grade


def test_rejects_empty_and_tiny_n():
    fn = lambda kd: 0.0  # noqa: E731
    for bad in ([], {}):
        try:
            monte_carlo_propagate(bad, fn, n=100)
            assert False
        except ValueError:
            pass
    try:
        monte_carlo_propagate({"A": _prior(-9.0)}, fn, n=1)
        assert False
    except ValueError:
        pass
