"""Cell-cycle phases and the kill ceiling of phase-specific drugs.

Anchors:
  * the returned λ actually satisfies the characteristic equation
    2·Π k/(k+d+λ) = 1 (independent check of the solver);
  * exponential dwell times grow FASTER than a fixed cycle length
    (λ > ln2/T_cycle) — variability lets some cells divide early;
  * phase fractions sum to 1 and S-phase share is sensible;
  * THE CEILING: an S-specific drug's kill rate saturates at the rate
    cells enter S — infinite cytotoxicity does not give infinite kill;
  * quiescent cells set a floor no cycle-specific drug can pass.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    CellCycle,
    growth_rate_per_h,
    stable_phase_fractions,
    evaluate_cycle_drug,
)

_LN2 = math.log(2.0)


def _char_eq(cycle, kills, lam):
    prod = 1.0
    for p, k in zip(cycle.phases, cycle.rates_per_h()):
        prod *= k / (k + kills.get(p, 0.0) + lam)
    return 2.0 * prod


def test_growth_rate_satisfies_characteristic_equation():
    cycle = CellCycle()
    for kills in ({}, {"S": 0.05}, {"G2M": 0.3}, {"G1": 0.02, "S": 0.1}):
        lam = growth_rate_per_h(cycle, kills)
        assert abs(_char_eq(cycle, kills, lam) - 1.0) < 1e-9


def test_exponential_dwell_beats_fixed_cycle_length():
    """Variability in phase duration makes the population grow faster
    than a deterministic cycle of the same mean length."""
    cycle = CellCycle()
    lam = growth_rate_per_h(cycle)
    naive = _LN2 / cycle.cycle_time_h
    assert lam > naive
    # ...but not absurdly so (same order of magnitude).
    assert lam < 3.0 * naive


def test_phase_fractions_are_a_distribution():
    cycle = CellCycle()
    fr = stable_phase_fractions(cycle)
    assert abs(sum(fr.values()) - 1.0) < 1e-12
    assert all(0.0 < v < 1.0 for v in fr.values())
    # G1 is the longest phase, so it holds the largest share.
    assert fr["G1"] == max(fr.values())


def test_phase_specific_kill_has_a_transit_limited_ceiling():
    """THE result: infinite S-phase cytotoxicity cannot beat the rate at
    which cells enter S (= 1/T_G1)."""
    cycle = CellCycle()
    k_G1 = 1.0 / cycle.phase_hours["G1"]
    lams = [growth_rate_per_h(cycle, {"S": d})
            for d in (1.0, 1e2, 1e4, 1e6, 1e9)]
    # Monotonically more negative, but converging.
    for a, b in zip(lams, lams[1:]):
        assert b < a
    assert lams[-1] > -k_G1                     # never passes the ceiling
    assert abs(lams[-1] + k_G1) < 1e-3          # converges to it
    # The reported ceiling matches.
    out = evaluate_cycle_drug(cycle, {"S": 1e6})
    assert abs(out.max_kill_rate_per_h - k_G1) < 1e-12


def test_only_a_fraction_of_cells_is_vulnerable_at_any_instant():
    cycle = CellCycle()
    out = evaluate_cycle_drug(cycle, {"S": 0.2})
    assert out.targeted_phases == ["S"]
    assert 0.05 < out.targeted_fraction < 0.6
    assert out.is_killing or out.growth_rate_per_h < \
        out.untreated_growth_rate_per_h


def test_non_specific_drug_has_no_such_ceiling():
    """Killing in EVERY phase escapes the transit limit."""
    cycle = CellCycle()
    k_G1 = 1.0 / cycle.phase_hours["G1"]
    everywhere = {p: 1e3 for p in cycle.phases}
    lam = growth_rate_per_h(cycle, everywhere)
    assert lam < -k_G1        # far past what an S-specific drug can reach


def test_quiescent_cells_set_a_floor():
    """G0 cells are refractory — a sanctuary in time."""
    cycle = CellCycle(quiescent_fraction=0.05)
    out = evaluate_cycle_drug(cycle, {"S": 5.0})
    assert out.is_killing
    # However long you treat, 5 % of the population survives.
    assert out.total_population_at(1e5) >= 0.05 - 1e-12
    assert out.log10_kill_at(1e5) < 1.31        # log10(1/0.05) = 1.30
    # With no quiescent pool the kill is unbounded.
    out2 = evaluate_cycle_drug(CellCycle(), {"S": 5.0})
    assert out2.log10_kill_at(1e5) > 100


def test_rejects_bad_inputs():
    try:
        CellCycle(phase_hours={"G1": 0.0})
        assert False
    except ValueError:
        pass
    try:
        CellCycle(quiescent_fraction=1.0)
        assert False
    except ValueError:
        pass
    try:
        evaluate_cycle_drug(CellCycle(), {"nonsense": 1.0})
        assert False
    except ValueError:
        pass
