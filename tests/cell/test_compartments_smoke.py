"""Two-compartment membrane permeation: drug outside → concentration inside.

Pins the passive-permeation physics (τ = r/3P for a sphere, C_in(t)
closed form, steady state → C_out) and the honesty carry-through: a
membrane-aware occupancy is decision-grade only if BOTH the K_d and the
permeability are calibrated, and the bulk-C_out reading is kept alongside
so the effect of the barrier is visible.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    Permeability,
    spherical_cell_geometry,
    equilibration_tau_s,
    intracellular_concentration,
    occupancy_in_cell,
)


def test_tau_matches_r_over_3P_for_a_sphere():
    """For a sphere, τ = V/(P·A) collapses to r/(3P)."""
    geom = spherical_cell_geometry(radius_um=10.0)
    perm = Permeability(P_cm_per_s=1e-5)
    tau = equilibration_tau_s(perm, geom)
    r_cm = 10.0 * 1e-4
    expected = r_cm / (3.0 * 1e-5)
    assert abs(tau - expected) / expected < 1e-9
    # ~33 s for these numbers — physically sensible.
    assert 30 < tau < 36


def test_concentration_rises_and_saturates_to_C_out():
    geom = spherical_cell_geometry(10.0)
    perm = Permeability(P_cm_per_s=1e-5)
    C_out = 1e-6
    tau = equilibration_tau_s(perm, geom)

    # At t = 0, nothing inside; at t = τ, ~63 %; at steady state, → C_out.
    assert intracellular_concentration(C_out, 0.0, perm, geom) == 0.0
    at_tau = intracellular_concentration(C_out, tau, perm, geom)
    assert abs(at_tau - C_out * (1 - math.exp(-1))) / C_out < 1e-9
    ss = intracellular_concentration(C_out, math.inf, perm, geom)
    assert abs(ss - C_out) < 1e-18


def test_lower_permeability_means_slower_equilibration():
    geom = spherical_cell_geometry(10.0)
    fast = Permeability(P_cm_per_s=1e-4)
    slow = Permeability(P_cm_per_s=1e-7)
    assert equilibration_tau_s(slow, geom) > equilibration_tau_s(fast, geom)


def test_membrane_matters_early_bulk_would_overstate_occupancy():
    """Early on, C_in << C_out, so a bulk reading overstates occupancy."""
    prior = binding_to_hill(-6.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    perm = Permeability(P_cm_per_s=1e-7, source="literature")  # slow
    geom = spherical_cell_geometry(10.0)
    tau = equilibration_tau_s(perm, geom)

    early = occupancy_in_cell(prior, 1e-4, perm, geom, t_s=tau * 0.01)
    assert early.C_in_M < early.C_out_M
    assert early.occupancy.theta < early.naive_occupancy_theta
    # At steady state the two converge (passive, no trapping).
    ss = occupancy_in_cell(prior, 1e-4, perm, geom, t_s=math.inf)
    assert abs(ss.occupancy.theta - ss.naive_occupancy_theta) < 1e-9


def test_uncalibrated_permeability_is_not_decision_grade():
    """Even a trustworthy K_d gives a non-decision-grade cell readout when
    the permeability is a placeholder (uncalibrated)."""
    prior = binding_to_hill(-6.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    assert prior.trust == "trustworthy_absolute"
    perm = Permeability(P_cm_per_s=1e-5, source="placeholder",
                        trust="uncalibrated")
    res = occupancy_in_cell(prior, 1e-6, perm, t_s=math.inf)
    assert not res.decision_grade      # permeability drags it down

    perm_ok = Permeability(P_cm_per_s=1e-5, source="martini",
                           trust="calibrated")
    res_ok = occupancy_in_cell(prior, 1e-6, perm_ok, t_s=math.inf)
    assert res_ok.decision_grade       # both K_d and P calibrated


def test_rejects_bad_inputs():
    for bad in (0.0, -1e-5):
        try:
            Permeability(P_cm_per_s=bad)
            assert False, "should reject non-positive P"
        except ValueError:
            pass
    try:
        spherical_cell_geometry(0.0)
        assert False, "should reject non-positive radius"
    except ValueError:
        pass
