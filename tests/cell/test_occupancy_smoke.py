"""First molecule→cell readout: target occupancy from a bridge prior.

Pins (a) the closed-form Hill occupancy, (b) exact CI propagation through
the monotonic K_d→θ relationship, and — the point of the whole bridge —
(c) that the prior's `trust` verdict rides all the way to the cell-level
readout, so an occupancy computed from an untrustworthy K_d is flagged
NOT decision-grade even though the number itself looks precise.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    hill_occupancy,
    occupancy_from_prior,
    OccupancyResult,
)


def test_occupancy_half_at_kd():
    """At [L] = K_d and n = 1, exactly half the target is bound."""
    assert abs(hill_occupancy(1e-8, 1e-8, 1.0) - 0.5) < 1e-12
    # Well above K_d → nearly saturated; well below → nearly empty.
    assert hill_occupancy(1e-6, 1e-8, 1.0) > 0.98
    assert hill_occupancy(1e-10, 1e-8, 1.0) < 0.02
    assert hill_occupancy(0.0, 1e-8, 1.0) == 0.0


def test_occupancy_ci_is_wider_when_kd_ci_is_wider():
    """A wider K_d CI must give a wider occupancy CI (exact, monotonic)."""
    # Same ΔG, but an untrustworthy target widens the K_d CI a lot.
    tight = binding_to_hill(-11.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    wide = binding_to_hill(-11.0, uncertainty_kcalmol=0.2,
                           receptor="benchmarks/dock/1stp.pdb", method="vina")
    L = 1e-8
    occ_tight = occupancy_from_prior(tight, L)
    occ_wide = occupancy_from_prior(wide, L)

    def width(occ):
        lo, hi = occ.theta_ci95
        return hi - lo

    assert width(occ_wide) > width(occ_tight)
    # CI must bracket the point estimate.
    for occ in (occ_tight, occ_wide):
        lo, hi = occ.theta_ci95
        assert lo <= occ.theta <= hi


def test_trust_rides_through_to_the_readout():
    """The whole point: an untrustworthy K_d yields a NOT-decision-grade
    occupancy, and a trustworthy one yields decision-grade."""
    untrust = binding_to_hill(-7.3, uncertainty_kcalmol=0.29,
                              receptor="benchmarks/dock/1stp.pdb",
                              method="vina")
    occ = occupancy_from_prior(untrust, 1e-8)
    assert occ.trust == "do_not_trust_absolute"
    assert not occ.decision_grade

    trust = binding_to_hill(-6.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb",
                            method="vina")
    occ2 = occupancy_from_prior(trust, 1e-6)
    assert occ2.trust == "trustworthy_absolute"
    assert occ2.decision_grade


def test_fep_readout_is_not_decision_grade_until_validated():
    """FEP-sourced K_d is uncalibrated for accuracy → readout not
    decision-grade, even though a CI (from precision) exists."""
    prior = binding_to_hill(-12.0, uncertainty_kcalmol=0.4,
                            method="FEP-DDM-MBAR")
    occ = occupancy_from_prior(prior, 1e-9)
    assert occ.theta_ci95 is not None       # precision CI present
    assert not occ.decision_grade           # but accuracy unknown
    assert occ.trust == "uncalibrated"


def test_no_ci_when_prior_has_none():
    prior = binding_to_hill(-10.0)          # no σ, no receptor
    occ = occupancy_from_prior(prior, 1e-8)
    assert occ.theta_ci95 is None
    assert isinstance(occ, OccupancyResult)
    assert "no CI" in occ.summary()


def test_rejects_non_hill_prior():
    from src.bridge import affinity_to_michaelis
    mm = affinity_to_michaelis(kcat_per_s=100.0, KM_M=1e-5)
    try:
        occupancy_from_prior(mm, 1e-8)
        assert False, "should reject a michaelis prior"
    except ValueError:
        pass
