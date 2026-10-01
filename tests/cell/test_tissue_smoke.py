"""Spatial penetration into tissue — cells at depth see less drug.

Rigor anchor: the nonlinear saturable BVP solver must reproduce the
analytic first-order profile in the C ≪ K_m limit (where MM uptake IS
first-order with k = V_max/K_m). Also pins the boundary conditions, the
monotone decay, the penetration-depth scaling, and the therapeutic
consequence (deep cells untreated).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    penetration_depth_um,
    penetration_profile_first_order,
    penetration_profile_saturable,
)


def test_penetration_depth_scaling():
    """λ = √(D/k): faster uptake → shallower penetration."""
    lam_slow = penetration_depth_um(1e-6, 1e-3)
    lam_fast = penetration_depth_um(1e-6, 1e-1)
    assert lam_slow > lam_fast
    # √10× uptake ⇒ λ down by √10.
    assert abs(lam_slow / lam_fast - math.sqrt(100)) < 1e-6


def test_analytic_profile_boundary_conditions_and_decay():
    p = penetration_profile_first_order(1e-6, thickness_um=150.0, k_per_s=1e-2)
    # C(0) = C_vessel exactly.
    assert abs(p.C_M[0] - 1e-6) / 1e-6 < 1e-12
    # Monotone decreasing into the tissue.
    for a, b in zip(p.C_M, p.C_M[1:]):
        assert b <= a + 1e-18
    # Zero-flux at the far edge → the last two points are nearly equal.
    assert abs(p.C_M[-1] - p.C_M[-2]) / p.C_M[0] < 1e-3


def test_saturable_reduces_to_first_order_when_C_far_below_Km():
    """THE anchor: at C ≪ K_m, MM uptake is first-order with k=Vmax/Km."""
    Km = 1e-3            # far above the drug concentration
    Vmax = 1e-5
    k_eff = Vmax / Km    # = 1e-2 /s
    C0 = 1e-9            # C ≪ Km

    num = penetration_profile_saturable(
        C0, thickness_um=150.0, Vmax_M_per_s=Vmax, Km_M=Km)
    ana = penetration_profile_first_order(
        C0, thickness_um=150.0, k_per_s=k_eff)

    for cn, ca in zip(num.C_M, ana.C_M):
        assert abs(cn - ca) / C0 < 1e-4


def test_saturating_dose_penetrates_deeper_than_first_order():
    """When C ≫ Km the uptake saturates (zeroth-order), so a high dose
    reaches deeper than the linear model would predict."""
    Km, Vmax = 1e-9, 1e-11
    thin = penetration_profile_saturable(
        1e-6, thickness_um=150.0, Vmax_M_per_s=Vmax, Km_M=Km)  # C >> Km
    low = penetration_profile_saturable(
        1e-12, thickness_um=150.0, Vmax_M_per_s=Vmax, Km_M=Km)  # C << Km
    assert thin.far_fraction > low.far_fraction


def test_deep_cells_are_untreated_when_lambda_is_small():
    """Therapeutic consequence: λ ≪ L leaves the far tissue untreated."""
    prior = binding_to_hill(-9.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    # Fast uptake → shallow λ relative to a 150 µm Krogh half-distance.
    p = penetration_profile_first_order(1e-5, thickness_um=150.0, k_per_s=1.0)
    assert p.penetration_depth_um < 50.0
    assert p.far_fraction < 0.05
    treated = p.treated_fraction(prior, threshold=0.5)
    assert 0.0 < treated < 0.6        # only the vessel-side band is treated

    # Slow uptake → deep λ → nearly the whole slab treated.
    q = penetration_profile_first_order(1e-5, thickness_um=150.0, k_per_s=1e-4)
    assert q.treated_fraction(prior, threshold=0.5) > 0.95


def test_rejects_bad_inputs():
    for bad in (dict(thickness_um=0.0, k_per_s=1e-2),
                dict(thickness_um=100.0, k_per_s=0.0)):
        try:
            penetration_profile_first_order(1e-6, **bad)
            assert False
        except ValueError:
            pass
