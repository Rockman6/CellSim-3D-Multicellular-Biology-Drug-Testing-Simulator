#!/usr/bin/env python3
"""Gates for the spatial dish (cellsim.cell.dish).

Two kinds. The field solvers are checked against closed forms, because
everything the dish says about oxygen and drug in tissue rests on them:
the zero-order sphere profile and diffusion limit used by Grimes et al.
2014 (J R Soc Interface 11:20131124) for DLD-1 spheroids, and the
first-order slab profile already in cellsim.cell.tissue. The lattice is
checked for the behaviours that make it a dish: exponential growth at
the line's doubling time while sparse, a plateau at confluence, oxygen
gradients and a necrotic core only where oxygen runs out, and clones
that inherit their mother's traits.
"""
from __future__ import annotations

import math
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import (APOPTOTIC, EMPTY, NECROTIC, Dish, DishParams,  # noqa: E402
                               run_dish, solve_steady, step_transient)
from cellsim.cell.engine import Params  # noqa: E402
from cellsim.cell.library import get_line  # noqa: E402
from cellsim.cell.tissue import penetration_profile_first_order  # noqa: E402

A549 = get_line("A549")


def test_oxygen_in_a_packed_sphere_matches_the_closed_form():
    """Zero-order consumption q in a sphere of radius R with surface
    pressure p0: p(r) = p0 - q (R^2 - r^2) / (6 D). The centre reaches
    zero at R = sqrt(6 D p0 / q), Grimes's diffusion limit: 233 um for
    their DLD-1 numbers."""
    dp = DishParams()
    D, p0, q = dp.o2_D_um2_per_s, dp.o2_medium_mmHg, dp.o2_uptake_mmHg_per_s
    r_limit = math.sqrt(6 * D * p0 / q)
    assert abs(r_limit - 233.0) < 1.0, r_limit
    R = 200.0
    edges = np.linspace(0, R, 401)
    p = solve_steady(edges, "sphere", D, p0, np.full(400, q), Km=1e-6)
    r = 0.5 * (edges[1:] + edges[:-1])
    exact = p0 - q * (R ** 2 - r ** 2) / (6 * D)
    assert np.max(np.abs(p - exact)) < 0.01 * p0, np.max(np.abs(p - exact))
    # just beyond the limit the centre is anoxic, just inside it is not
    for R_test, anoxic in ((0.95 * r_limit, False), (1.05 * r_limit, True)):
        e = np.linspace(0, R_test, 401)
        pc = solve_steady(e, "sphere", D, p0, np.full(400, q), Km=1e-3)[0]
        assert (pc < 1.0) == anoxic, (R_test, pc)


def test_first_order_slab_matches_the_tissue_module():
    """With Km far above the concentration, consumption is first order
    (k = uptake/Km) and the slab profile must be tissue.py's analytic
    cosh solution."""
    D_um2, k_per_s, L, c0 = 100.0, 1e-3, 300.0, 1.0
    edges = np.linspace(0, L, 601)                 # x' = distance from the no-flux face
    Km = 1e6
    c = solve_steady(edges, "plane", D_um2, c0, np.full(600, k_per_s * Km), Km=Km)
    ref = penetration_profile_first_order(c0 * 1e-6, thickness_um=L, k_per_s=k_per_s,
                                          D_cm2_per_s=D_um2 * 1e-8, n_points=601)
    x_from_vessel = L - 0.5 * (edges[1:] + edges[:-1])
    expect = np.interp(x_from_vessel, ref.x_um, np.array(ref.C_M) * 1e6)
    assert np.max(np.abs(c - expect)) < 0.005, np.max(np.abs(c - expect))


def test_transient_solver_relaxes_to_the_steady_profile():
    """Backward Euler with a linear sink must settle on the steady
    solution with the same sink, and with no sink on the boundary value."""
    edges = np.linspace(0, 250.0, 101)
    a = np.full(100, 2e-3)
    c = np.zeros(100)
    for _ in range(400):
        c = step_transient(edges, "sphere", 100.0, 1.0, c, 60.0, a, np.zeros(100), np.ones(100))
    steady = solve_steady(edges, "sphere", 100.0, 1.0, a * 1e6, Km=1e6)
    assert np.max(np.abs(c - steady)) < 1e-3, np.max(np.abs(c - steady))
    c = np.zeros(100)
    for _ in range(400):
        c = step_transient(edges, "sphere", 100.0, 1.0, c, 60.0, np.zeros(100), np.zeros(100),
                           np.full(100, 0.4))
    assert np.max(np.abs(c - 1.0)) < 1e-3


def test_sparse_monolayer_grows_at_the_doubling_time_then_stops_at_confluence():
    res = run_dish(A549, t_end_h=200.0, n_seed=60, grid_sites=40, record_every_h=12.0,
                   p=Params(cycle_cv=0.25), dp=DishParams(geometry="monolayer", push_sites=2))
    rec = res.record
    t, n = np.array(rec.t_h), np.array(rec.n_live, float)
    early = t <= 48.0
    doubling = math.log(2) / np.polyfit(t[early], np.log(n[early]), 1)[0]
    assert abs(doubling - A549.doubling_time_h) < 0.15 * A549.doubling_time_h, doubling
    assert rec.confluence[-1] > 0.9, rec.confluence[-1]
    assert n[-1] <= 40 * 40, "one cell per site"
    late_growth = n[-1] / n[-3]
    assert late_growth < 1.15, f"growth should stall at confluence, x{late_growth:.2f} in 24 h"
    assert rec.quiescent_frac[-1] > 0.5, rec.quiescent_frac[-1]


def test_lattice_bookkeeping_is_consistent_and_daughters_inherit():
    dish = Dish(A549, n_seed=30, grid_sites=24, seed=3, p=Params(cycle_cv=0.25),
                dp=DishParams(geometry="monolayer"))
    traits = {int(l): (b, u) for l, b, u in zip(dish.lineage, dish.cells.het_bcl2,
                                                dish.cells.het_uptake)}
    for _ in range(int(48 / dish.dp.field_dt_h)):
        dish.advance(None)
    n = len(dish.pos)
    assert n > 30, "the colony should have grown"
    # every live cell sits on its own site, and the site points back at it
    assert np.array_equal(dish.occ[tuple(dish.pos.T)], np.arange(n))
    assert len({tuple(x) for x in dish.pos}) == n
    assert np.sum(dish.occ >= 0) == n
    # every cell carries its founder's reserve and uptake exactly
    for l, b, u in zip(dish.lineage, dish.cells.het_bcl2, dish.cells.het_uptake):
        assert traits[int(l)] == (b, u)
    # complete cycles were logged, with plausible durations
    durations = np.array([d - b for b, d in dish.cycles])
    assert len(durations) > 5 and 10.0 < durations.mean() < 40.0, durations.mean()


def test_a_spheroid_goes_anoxic_and_necrotic_only_at_its_centre():
    """Consumption raised ninefold puts the diffusion limit at 78 um, so a
    2000-cell aggregate (R ~ 117 um) must have an anoxic centre, form a
    necrotic core there within a day, and keep its rim alive."""
    dp = replace(DishParams(geometry="spheroid", push_sites=float("inf")),
                 o2_uptake_mmHg_per_s=9 * 22.1)
    res = run_dish(get_line("DLD-1"), t_end_h=24.0, n_seed=2000, grid_sites=36,
                   record_every_h=6.0, seed=2, p=Params(cycle_cv=0.25), dp=dp)
    rec, dish = res.record, res.dish
    assert rec.min_o2_mmHg[0] < dp.o2_necrosis_full_mmHg, rec.min_o2_mmHg[0]
    assert rec.n_necrotic[-1] > 50, rec.n_necrotic[-1]
    assert rec.hypoxic_frac[-1] > 0.2, rec.hypoxic_frac[-1]
    c = dish.centroid()
    necrotic = np.argwhere(dish.occ == NECROTIC)
    d_nec = np.linalg.norm(necrotic - c, axis=1).mean()
    d_live = np.linalg.norm(dish.pos - c, axis=1).mean()
    assert d_nec < 0.7 * d_live, (d_nec, d_live)
    # the outer cells see medium-level oxygen
    outer = np.linalg.norm(dish.pos - c, axis=1) > np.percentile(np.linalg.norm(dish.pos - c, axis=1), 90)
    assert dish.o2[outer].mean() > 20.0, dish.o2[outer].mean()


def test_the_dish_is_deterministic_for_a_seed():
    a = run_dish(A549, t_end_h=24.0, n_seed=40, grid_sites=20, seed=5, record_every_h=24.0)
    b = run_dish(A549, t_end_h=24.0, n_seed=40, grid_sites=20, seed=5, record_every_h=24.0)
    assert a.record.n_live == b.record.n_live
    assert np.array_equal(a.dish.pos, b.dish.pos)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} dish gates pass")
    raise SystemExit(1 if failed else 0)
