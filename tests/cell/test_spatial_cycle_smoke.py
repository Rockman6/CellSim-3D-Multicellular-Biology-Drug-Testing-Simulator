"""Spatial × cycle coupling — cells on a lattice that carry a cycle phase.

Anchors:
  * THE tie-in: with a cycle and no drug, the STOCHASTIC SPATIAL model
    reproduces the analytic phase distribution from `cycle.py`
    (`stable_phase_fractions`) — two independent formalisms agreeing;
  * a phase-specific drug spares cells outside its target phase, so it
    leaves more survivors than a non-specific drug of the same strength;
  * the DOUBLE SANCTUARY: cells deep in tissue see less drug AND cycle
    more slowly, and a phase-specific agent fails against them twice over;
  * cycle=None reproduces the original behaviour exactly.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    CellCycle, simulate_colony, stable_phase_fractions,
)


def _prior():
    return binding_to_hill(-9.0, uncertainty_kcalmol=0.2,
                           receptor="benchmarks/dock/3ptb.pdb", method="vina")


def test_agent_phase_distribution_matches_analytic_cycle_theory():
    """The stochastic lattice reproduces cycle.py's stable distribution."""
    cycle = CellCycle()
    analytic = stable_phase_fractions(cycle)
    tr = simulate_colony(grid_size=(220, 220), n_seed_cells=3000,
                         cycle=cycle, dt_h=0.5, n_steps=120,
                         rng_seed=11, seed_scattered=True)
    got = tr.phase_fractions()
    assert tr.final_count > 3000                    # it grew
    for p in cycle.phases:
        assert abs(got[p] - analytic[p]) < 0.05, (p, got[p], analytic[p])


def test_cycle_none_is_unchanged():
    a = simulate_colony(grid_size=(30, 30), n_seed_cells=5, n_steps=40,
                        rng_seed=2)
    b = simulate_colony(grid_size=(30, 30), n_seed_cells=5, n_steps=40,
                        rng_seed=2, cycle=None)
    assert a.n_total == b.n_total
    assert a.phase_fractions() == {}


def test_phase_specific_drug_spares_cells_outside_the_phase():
    """Same kill strength, but only S-phase cells are targets."""
    cycle = CellCycle()
    prior = _prior()
    common = dict(grid_size=(70, 70), n_seed_cells=1500, cycle=cycle,
                  prior=prior, drug_field=lambda x, y: 1e-5,
                  dt_h=1.0, n_steps=60, rng_seed=4, seed_scattered=True)
    nonspecific = simulate_colony(**common)
    s_only = simulate_colony(phase_specific=["S"], **common)
    assert s_only.final_count > nonspecific.final_count
    # The survivors of the S-specific drug are depleted in S.
    fr = s_only.phase_fractions()
    assert fr["S"] < stable_phase_fractions(cycle)["S"]


def test_double_sanctuary_deep_cells_resist_a_phase_specific_drug():
    """Deep cells see less drug AND cycle slower — protected twice."""
    cycle = CellCycle()
    prior = _prior()
    nx = 60

    def drug(x, y):                     # decays away from the vessel
        return 1e-5 * math.exp(-x / 10.0)

    def speed(x, y):                    # so does nutrient supply
        return 0.15 + 0.85 * math.exp(-x / 15.0)

    tr = simulate_colony(grid_size=(nx, 40), n_seed_cells=1600,
                         cycle=cycle, phase_specific=["S"], prior=prior,
                         drug_field=drug, cycle_rate_field=speed,
                         dt_h=1.0, n_steps=80, rng_seed=6)
    assert not tr.extinct
    xs = [x for x, _, _ in tr.occupied]
    near = sum(1 for x in xs if x < 12)
    far = sum(1 for x in xs if x >= 36)
    assert far > near


def test_rejects_bad_cycle_inputs():
    cycle = CellCycle()
    try:
        simulate_colony(grid_size=(10, 10), phase_specific=["S"])
        assert False, "phase_specific without a cycle should raise"
    except ValueError:
        pass
    try:
        simulate_colony(grid_size=(10, 10), cycle=cycle,
                        phase_specific=["nonsense"])
        assert False
    except ValueError:
        pass
    try:
        simulate_colony(grid_size=(10, 10),
                        cycle_rate_field=lambda x, y: 1.0)
        assert False, "cycle_rate_field without a cycle should raise"
    except ValueError:
        pass
