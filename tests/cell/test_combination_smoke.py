"""Combination therapy — resistance multiplies, but scheduling can antagonise.

Anchors:
  * independent resistance MULTIPLIES (f_a·f_b) and shared resistance
    collapses to min(f_a, f_b) — the quantitative case for combinations;
  * a cytostatic G1 block ANTAGONISES an S-phase drug, because the arrest
    lowers k_G1 and the phase-specific kill ceiling IS −k_G1 (this is
    predicted by the cycle model, not bolted on);
  * sequencing the same two drugs beats giving them together, for exactly
    that reason;
  * log-kills add across blocks and the quasi-steady-state warning fires
    when a block is shorter than a cycle.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    CellCycle, DrugAction, ScheduleBlock,
    apply_drugs, evaluate_combination, evaluate_schedule,
    double_resistant_fraction, growth_rate_per_h,
)

# An S-phase cytotoxic (antimetabolite-like) and a G1 arrest (cytostatic).
_S_KILLER = DrugAction("Sdrug", kill_per_h={"S": 1.0})
_G1_BLOCK = DrugAction("G1block", arrest_phase="G1", arrest_strength=0.9)


def test_independent_resistance_multiplies():
    assert abs(double_resistant_fraction(1e-4, 1e-4) - 1e-8) < 1e-20
    # Shared mechanism → no benefit beyond the rarer clone.
    assert abs(double_resistant_fraction(1e-4, 1e-6, cross_resistance=1.0)
               - 1e-6) < 1e-20
    # Partial cross-resistance sits between the two.
    mid = double_resistant_fraction(1e-4, 1e-6, cross_resistance=0.5)
    assert 1e-10 < mid < 1e-6


def test_arrest_lengthens_the_phase_and_lowers_the_ceiling():
    cycle = CellCycle()
    eff, kills = apply_drugs(cycle, [_G1_BLOCK, _S_KILLER])
    # 90 % arrest → G1 takes 10× longer.
    assert abs(eff.phase_hours["G1"] - 10 * cycle.phase_hours["G1"]) < 1e-9
    assert kills == {"S": 1.0}
    # The kill ceiling is now the SLOWED G1 rate.
    out = evaluate_combination(cycle, [_G1_BLOCK, _S_KILLER])
    assert abs(out.max_kill_rate_per_h - 1.0 / eff.phase_hours["G1"]) < 1e-12


def test_cytostatic_antagonises_the_phase_specific_drug():
    """THE result: adding a G1 block makes the S-drug WORSE, not better."""
    cycle = CellCycle()
    lam_s_alone = growth_rate_per_h(cycle, {"S": 1.0})
    lam_both = evaluate_combination(cycle, [_G1_BLOCK, _S_KILLER]).growth_rate_per_h
    assert lam_s_alone < 0                     # S-drug alone kills
    assert lam_both > lam_s_alone              # adding the block kills LESS
    # And the block alone does not kill at all.
    lam_block = evaluate_combination(cycle, [_G1_BLOCK]).growth_rate_per_h
    assert lam_block > lam_s_alone


def test_sequential_beats_concurrent_for_this_pair():
    """Same two drugs, same total time — order changes the outcome."""
    cycle = CellCycle()
    hours = 120.0
    concurrent = evaluate_schedule(cycle, [
        ScheduleBlock([_G1_BLOCK, _S_KILLER], hours)])
    sequential = evaluate_schedule(cycle, [
        ScheduleBlock([_S_KILLER], hours / 2),
        ScheduleBlock([_G1_BLOCK], hours / 2)])
    assert sequential.log10_kill > concurrent.log10_kill
    # Both schedules ran the same wall-clock time.
    assert abs(sequential.total_hours - concurrent.total_hours) < 1e-9


def test_log_kills_add_across_blocks():
    cycle = CellCycle()
    a = ScheduleBlock([_S_KILLER], 48.0)
    b = ScheduleBlock([], 48.0)                # drug holiday: regrowth
    combined = evaluate_schedule(cycle, [a, b])
    lam_a, lam_b = combined.block_rates_per_h
    expected = -(lam_a * 48.0 + lam_b * 48.0) / math.log(10.0)
    assert abs(combined.log10_kill - expected) < 1e-9
    assert lam_b > 0                            # regrows during the holiday


def test_quiescent_floor_applies_to_schedules():
    cycle = CellCycle(quiescent_fraction=0.02)
    out = evaluate_schedule(cycle, [ScheduleBlock([_S_KILLER], 5000.0)])
    assert out.surviving_fraction >= 0.02 - 1e-12
    assert out.log10_kill < -math.log10(0.02) + 1e-6


def test_short_block_raises_the_quasi_steady_state_warning():
    cycle = CellCycle()
    out = evaluate_schedule(cycle, [ScheduleBlock([_S_KILLER], 2.0)])
    assert not out.quasi_steady_state_ok
    assert "synchronisation" in out.warning
    long = evaluate_schedule(cycle, [ScheduleBlock([_S_KILLER], 100.0)])
    assert long.quasi_steady_state_ok and long.warning == ""


def test_rejects_bad_inputs():
    try:
        DrugAction("x", arrest_strength=0.5)     # no phase given
        assert False
    except ValueError:
        pass
    try:
        apply_drugs(CellCycle(), [DrugAction("y", kill_per_h={"Z": 1.0})])
        assert False
    except ValueError:
        pass
    try:
        double_resistant_fraction(0.5, 1.5)
        assert False
    except ValueError:
        pass
