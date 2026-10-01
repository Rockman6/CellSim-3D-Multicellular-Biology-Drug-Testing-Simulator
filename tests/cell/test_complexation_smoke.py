"""Drug-drug complexation A + B <=> AB: two drugs stick and form a new species.

Pins the bimolecular quadratic (mass balance + K_d identity), the limits
(tight K_d → nearly all complexed; loose K_d → little complex), monotonic
CI over K_d, and trust ride-through.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import complex_equilibrium  # noqa: E402


def _pair_prior(dG, rec="benchmarks/dock/3ptb.pdb", sigma=0.2):
    return binding_to_hill(dG, uncertainty_kcalmol=sigma, receptor=rec,
                           method="vina")


def test_mass_balance_and_kd_identity_hold():
    prior = _pair_prior(-9.0)
    res = complex_equilibrium(prior, A_total_M=1e-6, B_total_M=1e-6)
    # Mass balance: free + complex = total.
    assert abs((res.free_A_M + res.complex_M) - res.A_total_M) < 1e-18
    assert abs((res.free_B_M + res.complex_M) - res.B_total_M) < 1e-18
    # K_d identity: [A][B]/[AB] = K_d.
    kd_check = res.free_A_M * res.free_B_M / res.complex_M
    assert abs(kd_check - res.Kd_M) / res.Kd_M < 1e-6


def test_tight_binding_complexes_most_of_the_limiting_reagent():
    """Very tight K_d with excess B → almost all A is complexed."""
    tight = _pair_prior(-14.0)   # sub-nM
    res = complex_equilibrium(tight, A_total_M=1e-7, B_total_M=1e-5)
    assert res.complexed_fraction_A > 0.95


def test_weak_binding_leaves_little_complex():
    weak = _pair_prior(-4.0)     # ~mM K_d
    res = complex_equilibrium(weak, A_total_M=1e-8, B_total_M=1e-8)
    assert res.complexed_fraction_A < 0.05


def test_complex_ci_brackets_and_is_monotonic_in_kd():
    prior = _pair_prior(-9.0, sigma=0.3)
    res = complex_equilibrium(prior, 1e-6, 1e-6)
    assert res.complex_ci95 is not None
    lo, hi = res.complex_ci95
    assert lo <= res.complex_M <= hi
    assert lo >= 0.0


def test_trust_rides_through():
    good = _pair_prior(-9.0, rec="benchmarks/dock/3ptb.pdb")
    bad = _pair_prior(-9.0, rec="benchmarks/dock/1stp.pdb")
    assert complex_equilibrium(good, 1e-6, 1e-6).decision_grade
    assert not complex_equilibrium(bad, 1e-6, 1e-6).decision_grade


def test_rejects_non_hill_and_negative():
    from src.bridge import affinity_to_michaelis
    mm = affinity_to_michaelis(kcat_per_s=1.0, KM_M=1e-5)
    try:
        complex_equilibrium(mm, 1e-6, 1e-6)
        assert False
    except ValueError:
        pass
    try:
        complex_equilibrium(_pair_prior(-9.0), -1.0, 1e-6)
        assert False
    except ValueError:
        pass
