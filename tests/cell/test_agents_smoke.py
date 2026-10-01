"""Agent-based cells — individuals dividing, dying, competing for space.

Anchors:
  * determinism given a seed;
  * the STOCHASTIC simulation reproduces the ANALYTIC exponential in
    expectation (dilute regime, so contact inhibition is not yet biting) —
    the bridge between the agent layer and every rate equation above it;
  * pure-death regime decays as exp(−k·t);
  * contact inhibition emerges: a colony saturates at grid capacity
    instead of growing exponentially forever;
  * a drug gradient produces spatial structure — survivors sit away from
    the vessel edge.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    Clone, CellFateParams, simulate_colony,
)

_LN2 = math.log(2.0)


def test_deterministic_given_seed():
    a = simulate_colony(grid_size=(40, 40), n_seed_cells=10, n_steps=40,
                        rng_seed=3)
    b = simulate_colony(grid_size=(40, 40), n_seed_cells=10, n_steps=40,
                        rng_seed=3)
    c = simulate_colony(grid_size=(40, 40), n_seed_cells=10, n_steps=40,
                        rng_seed=4)
    assert a.n_total == b.n_total
    assert sorted(a.occupied) == sorted(b.occupied)
    assert a.n_total != c.n_total


def test_stochastic_growth_matches_analytic_exponential():
    """Dilute colony (no crowding) must follow N₀·exp(k·t) on average."""
    p = CellFateParams()                       # 24 h doubling
    k_per_h = p.k_prolif_per_s * 3600.0
    n0, hours = 200, 24.0
    # Big grid so free neighbours are always available.
    runs = [simulate_colony(grid_size=(300, 300), n_seed_cells=n0,
                            params=p, dt_h=0.5, n_steps=int(hours / 0.5),
                            rng_seed=s, seed_scattered=True).final_count
            for s in range(6)]
    mean = sum(runs) / len(runs)
    expected = n0 * math.exp(k_per_h * hours)   # = 2·n0 after one doubling
    assert abs(mean - expected) / expected < 0.10


def test_pure_death_decays_exponentially():
    """No division, only death → N(t) = N₀·exp(−k·t)."""
    p = CellFateParams(k_prolif_per_s=0.0,
                       k_death_basal_per_s=_LN2 / (12 * 3600.0))
    n0 = 800
    tr = simulate_colony(grid_size=(120, 120), n_seed_cells=n0, params=p,
                         dt_h=1.0, n_steps=12, rng_seed=1)
    # One half-life (12 h) → about half remain.
    assert abs(tr.final_count / n0 - 0.5) < 0.08


def test_contact_inhibition_saturates_the_colony():
    """A small grid fills up and stops — surface growth, not exponential."""
    tr = simulate_colony(grid_size=(12, 12), n_seed_cells=1,
                         dt_h=2.0, n_steps=250, rng_seed=2)
    capacity = 12 * 12
    assert tr.final_count > 0.85 * capacity
    assert tr.final_count <= capacity
    # Growth stalled: the last quarter of the run added few cells.
    q = len(tr.n_total) // 4
    assert tr.n_total[-1] - tr.n_total[-q] < 0.10 * capacity


def test_drug_gradient_creates_spatial_survivors():
    """Vessel at x=0: cells near it die, cells far from it survive."""
    prior = binding_to_hill(-9.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    nx = 60

    def field(x, y):
        # Exponential decay away from the vessel edge.
        return 1e-5 * math.exp(-x / 8.0)

    tr = simulate_colony(grid_size=(nx, 40), n_seed_cells=1200,
                         prior=prior, drug_field=field,
                         dt_h=2.0, n_steps=60, rng_seed=5)
    assert not tr.extinct
    xs = [x for x, _, _ in tr.occupied]
    near = sum(1 for x in xs if x < 10)
    far = sum(1 for x in xs if x >= 30)
    assert far > near, "survivors should sit away from the vessel"


def test_resistant_clone_sweeps_under_drug():
    """Two clones on the grid: the resistant one takes over under drug."""
    prior = binding_to_hill(-9.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    clones = [Clone("sensitive", occupancy_scale=1.0),
              Clone("resistant", occupancy_scale=0.03, fitness_cost=0.1)]
    tr = simulate_colony(grid_size=(70, 70), n_seed_cells=900,
                         clones=clones, seed_fractions=[0.98, 0.02],
                         prior=prior, drug_field=lambda x, y: 1e-5,
                         dt_h=2.0, n_steps=60, rng_seed=7,
                         seed_scattered=True)
    assert tr.clone_fraction("resistant") > 0.5


def test_rejects_bad_inputs():
    try:
        simulate_colony(grid_size=(10, 10), drug_field=lambda x, y: 1e-6)
        assert False, "drug field without a prior should raise"
    except ValueError:
        pass
    try:
        Clone("x", occupancy_scale=0.0)
        assert False
    except ValueError:
        pass
