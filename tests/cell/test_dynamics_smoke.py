"""Time-resolved binding kinetics — and the rigor checks that anchor it.

The ODE is only trustworthy if it reduces to the equilibria the analytic
modules compute. These tests pin:
  * k_off = k_on · K_d (kinetics from thermodynamics);
  * mass conservation T + C = T_tot throughout the trajectory;
  * the ODE steady state EQUALS the analytic Hill occupancy at L_out —
    the cross-check that ties dynamics to equilibrium;
  * monotone approach and permeability setting the timescale;
  * the equilibrium trust (K_d-only) is separate from kinetic calibration
    (k_on), so an uncalibrated k_on never taints the equilibrium verdict.
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
    occupancy_from_prior,
    rate_constants_from_prior,
    simulate_binding_in_cell,
)


def _prior(dG, rec="benchmarks/dock/3ptb.pdb", sigma=0.2):
    return binding_to_hill(dG, uncertainty_kcalmol=sigma, receptor=rec,
                           method="vina")


def test_koff_equals_kon_times_kd():
    p = _prior(-9.0)
    rc = rate_constants_from_prior(p, kon_per_M_per_s=1e6)
    assert abs(rc.koff_per_s - 1e6 * rc.Kd_M) / rc.koff_per_s < 1e-12
    assert not rc.kinetics_calibrated             # default k_on is a convention


def test_steady_state_matches_analytic_hill_occupancy():
    """THE consistency check: ODE t→∞ occupancy == analytic Hill at L_out."""
    p = _prior(-8.0)
    perm = Permeability(P_cm_per_s=1e-5)
    tc = simulate_binding_in_cell(
        p, target_conc_M=1e-9, drug_out_M=1e-7, permeability=perm)
    analytic = occupancy_from_prior(p, 1e-7).theta   # Hill at L_out
    assert abs(tc.equilibrium_occupancy - analytic) < 1e-12
    # ODE actually converged to it.
    assert abs(tc.converged_occupancy - tc.equilibrium_occupancy) < 1e-3


def test_mass_is_conserved_along_the_trajectory():
    p = _prior(-10.0)
    perm = Permeability(P_cm_per_s=1e-5)
    tc = simulate_binding_in_cell(
        p, target_conc_M=5e-9, drug_out_M=1e-7, permeability=perm)
    assert tc.max_mass_imbalance < 1e-6


def test_occupancy_starts_at_zero_and_rises_monotonically():
    p = _prior(-9.0)
    perm = Permeability(P_cm_per_s=1e-5)
    tc = simulate_binding_in_cell(
        p, target_conc_M=1e-9, drug_out_M=1e-7, permeability=perm)
    assert tc.occupancy[0] < 1e-9
    # Non-decreasing approach (allow tiny numerical wiggle).
    for a, b in zip(tc.occupancy, tc.occupancy[1:]):
        assert b >= a - 1e-9


def test_lower_permeability_slows_the_approach():
    p = _prior(-9.0)
    geom = spherical_cell_geometry(10.0)
    fast = simulate_binding_in_cell(
        p, target_conc_M=1e-9, drug_out_M=1e-7,
        permeability=Permeability(P_cm_per_s=1e-4), geometry=geom)
    slow = simulate_binding_in_cell(
        p, target_conc_M=1e-9, drug_out_M=1e-7,
        permeability=Permeability(P_cm_per_s=1e-7), geometry=geom)
    # Both reach the same equilibrium...
    assert abs(fast.equilibrium_occupancy - slow.equilibrium_occupancy) < 1e-12
    # ...but the low-permeability cell takes longer to get halfway.
    assert slow.time_to_occupancy(0.5) > fast.time_to_occupancy(0.5)


def test_equilibrium_trust_is_independent_of_kon_calibration():
    """A trustworthy K_d gives an equilibrium-decision-grade result even
    though the default k_on (timescale) is an uncalibrated convention."""
    good = _prior(-8.0, rec="benchmarks/dock/3ptb.pdb")
    perm = Permeability(P_cm_per_s=1e-5)
    tc = simulate_binding_in_cell(
        good, target_conc_M=1e-9, drug_out_M=1e-7, permeability=perm)
    assert tc.equilibrium_decision_grade          # K_d is trustworthy
    assert not tc.rates.kinetics_calibrated        # but the timescale isn't

    bad = _prior(-8.0, rec="benchmarks/dock/1stp.pdb")
    tc_bad = simulate_binding_in_cell(
        bad, target_conc_M=1e-9, drug_out_M=1e-7, permeability=perm)
    assert not tc_bad.equilibrium_decision_grade


def test_occupancy_at_interpolates():
    p = _prior(-9.0)
    perm = Permeability(P_cm_per_s=1e-5)
    tc = simulate_binding_in_cell(
        p, target_conc_M=1e-9, drug_out_M=1e-7, permeability=perm)
    mid_t = tc.t_s[len(tc.t_s) // 2]
    assert abs(tc.occupancy_at(mid_t) - tc.occupancy[len(tc.t_s) // 2]) < 1e-9
    # Clamps outside range.
    assert tc.occupancy_at(-1.0) == tc.occupancy[0]
    assert tc.occupancy_at(1e30) == tc.occupancy[-1]
