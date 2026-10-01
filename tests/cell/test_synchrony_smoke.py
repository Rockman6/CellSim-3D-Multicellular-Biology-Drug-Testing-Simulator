"""Full cycle transient — cohort structure and dose timing.

Anchors:
  * THE tie-in: under a constant drug the transient's asymptotic growth
    rate converges to the characteristic-equation λ from `cycle.py`, and
    its phase distribution converges to `stable_phase_fractions` —
    independent numerical confirmation of the analytic eigenvalue;
  * for LONG blocks the transient agrees with the quasi-steady-state
    `evaluate_schedule`; for SHORT blocks they diverge — so the QSS
    warning is now quantified, not just asserted;
  * an S-phase pulse SYNCHRONISES the survivors (synchrony rises from ~0);
  * dose timing matters — the best gap between two pulses beats the worst,
    which the QSS model structurally cannot represent.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    CellCycle, DrugAction, ScheduleBlock,
    growth_rate_per_h, stable_phase_fractions, evaluate_schedule,
    simulate_schedule_transient, best_second_dose_delay,
)

_S_KILLER = DrugAction("Sdrug", kill_per_h={"S": 1.0})


def test_transient_converges_to_characteristic_growth_rate():
    """THE anchor: numerical ODE reproduces the analytic eigenvalue."""
    cycle = CellCycle()
    for kills, drugs in (({}, []), ({"S": 0.3}, [DrugAction("d", {"S": 0.3})])):
        lam_analytic = growth_rate_per_h(cycle, kills)
        tr = simulate_schedule_transient(
            cycle, [ScheduleBlock(drugs, 400.0)], points_per_block=4000)
        # Slope of log(total) over the last stretch = asymptotic λ.
        i0 = int(0.8 * len(tr.t_h))
        lam_num = (math.log(tr.total[-1]) - math.log(tr.total[i0])) / \
            (tr.t_h[-1] - tr.t_h[i0])
        assert abs(lam_num - lam_analytic) < 1e-4


def test_phase_distribution_converges_to_stable():
    cycle = CellCycle()
    stable = stable_phase_fractions(cycle)
    tr = simulate_schedule_transient(
        cycle, [ScheduleBlock([], 500.0)], points_per_block=2000,
        initial_fractions={"G1": 1.0, "S": 0.0, "G2M": 0.0})
    final = tr.phase_fractions_at(-1)
    for p in cycle.phases:
        assert abs(final[p] - stable[p]) < 1e-3
    # Starting fully synchronised, synchrony decays away.
    assert tr.synchrony[0] > 0.4
    assert tr.synchrony[-1] < 1e-2


def test_long_blocks_agree_with_quasi_steady_state():
    """QSS is valid where it claims to be."""
    cycle = CellCycle()
    blocks = [ScheduleBlock([_S_KILLER], 240.0), ScheduleBlock([], 240.0)]
    qss = evaluate_schedule(cycle, blocks)
    tr = simulate_schedule_transient(cycle, blocks, points_per_block=1500)
    assert qss.quasi_steady_state_ok
    assert abs(tr.log10_kill - qss.log10_kill) < 0.25


def test_short_blocks_diverge_from_quasi_steady_state():
    """And it is wrong where it warns — now measurably so."""
    cycle = CellCycle()
    blocks = [ScheduleBlock([_S_KILLER], 3.0), ScheduleBlock([], 3.0)] * 12
    qss = evaluate_schedule(cycle, blocks)
    tr = simulate_schedule_transient(cycle, blocks, points_per_block=120)
    assert not qss.quasi_steady_state_ok
    assert abs(tr.log10_kill - qss.log10_kill) > 0.05


def test_s_phase_pulse_synchronises_the_survivors():
    """Killing S-phase cells leaves a cohort marching in step."""
    cycle = CellCycle()
    tr = simulate_schedule_transient(
        cycle, [ScheduleBlock([DrugAction("pulse", {"S": 5.0})], 6.0),
                ScheduleBlock([], 30.0)], points_per_block=400)
    assert tr.synchrony[0] < 1e-9        # starts asynchronous
    assert tr.peak_synchrony > 0.15      # pulse creates a cohort


def test_dose_timing_changes_the_kill():
    """A question QSS cannot answer: when should the second dose land?"""
    cycle = CellCycle()
    scan = best_second_dose_delay(cycle, DrugAction("pulse", {"S": 5.0}),
                                  pulse_hours=2.0, points_per_block=80)
    assert scan.timing_benefit_log10 > 0.02
    assert scan.best_delay_h != scan.worst_delay_h
    assert scan.best_log10_kill > scan.worst_log10_kill

    # The QSS model cannot see this. For it, a longer gap is ALWAYS worse
    # (more regrowth), so its best gap is trivially the shortest one and
    # its kill is monotonic in the gap. The transient is not monotonic —
    # that non-monotonicity IS the cohort structure.
    def qss_kill(gap):
        pulse = DrugAction("pulse", {"S": 5.0})
        blocks = [ScheduleBlock([pulse], 2.0)]
        if gap > 0:
            blocks.append(ScheduleBlock([], gap))
        blocks.append(ScheduleBlock([pulse], 2.0))
        return evaluate_schedule(cycle, blocks).log10_kill

    gaps = [g for g in scan.delays_h if g > 0]
    qss_vals = [qss_kill(g) for g in gaps]
    assert all(b <= a + 1e-12 for a, b in zip(qss_vals, qss_vals[1:])), \
        "QSS kill must decrease monotonically with the gap"
    tr_vals = [k for g, k in zip(scan.delays_h, scan.log10_kill) if g > 0]
    assert any(b > a + 1e-9 for a, b in zip(tr_vals, tr_vals[1:])), \
        "transient must be non-monotonic in the gap — the cohort effect"


def test_quiescent_floor_applies():
    cycle = CellCycle(quiescent_fraction=0.03)
    tr = simulate_schedule_transient(
        cycle, [ScheduleBlock([DrugAction("d", {"S": 5.0})], 1000.0)],
        points_per_block=800)
    assert tr.total_with_quiescent(-1) >= 0.03 - 1e-12
    assert tr.log10_kill < -math.log10(0.03) + 1e-6


def test_rejects_bad_inputs():
    try:
        simulate_schedule_transient(CellCycle(), [])
        assert False
    except ValueError:
        pass
    try:
        simulate_schedule_transient(
            CellCycle(), [ScheduleBlock([], 10.0)],
            initial_fractions={"G1": 1.0})
        assert False
    except ValueError:
        pass
