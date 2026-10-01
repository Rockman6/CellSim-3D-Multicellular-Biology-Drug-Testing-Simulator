#!/usr/bin/env python3
"""Phase-1 engine gates: the behaviours the cell model must reproduce.

Each test asserts a published phenotype or an internal consistency
property, not merely that the code runs. Seconds per test on CPU.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import (  # noqa: E402
    IX, Params, calibrate_cycle_scale, dose_response, ic50_from_curve, simulate,
)
from cellsim.cell.library import get_drug, get_line  # noqa: E402

A549 = get_line("A549")
K_A549 = calibrate_cycle_scale(A549)


def test_untreated_population_doubles_at_the_quoted_rate():
    """An untreated colony must grow at the line's doubling time, which is
    what the single cycle-rate scale is calibrated against."""
    res = simulate(A549, None, 0.0, t_end_h=72, n_cells=128, dt_h=0.02, k_cyc=K_A549)
    fold = res.alive_weight[-1, 0] / res.alive_weight[0, 0]
    expected = 2 ** (72 / A549.doubling_time_h)
    assert 0.75 * expected < fold < 1.35 * expected, (
        f"72 h growth x{fold:.2f}, expected ~x{expected:.2f}")
    assert res.deaths[-1, 0] == 0, "untreated cells must not die"


def test_untreated_phase_fractions_are_physiological():
    """Asynchronous culture: G1 should dominate, M should be a few percent
    (Casciari 1992 reports ~46/33/17/4 for HeLa; we gate loosely)."""
    res = simulate(A549, None, 0.0, t_end_h=48, n_cells=128, dt_h=0.02, k_cyc=K_A549)
    g1, s, g2, m = res.phase_frac[-1]
    assert g1 > 0.4, f"G1 fraction {g1:.2f} too low"
    assert m < 0.15, f"M fraction {m:.2f} too high"
    assert abs(g1 + s + g2 + m - 1.0) < 1e-6


def test_p53_pulses_with_a_purvis_like_period_under_sustained_damage():
    """Purvis & Lahav 2012 Science 336:1440 — sustained DNA damage drives
    p53 pulses with a ~5.5 h period. The frozen C++ prototype gave 3.0 h;
    the port is time-scaled to land in the published window."""
    res = simulate(A549, None, 0.0, t_end_h=48, n_cells=16, dt_h=0.02,
                   k_cyc=K_A549, forced_D=0.30, record_every_h=0.1)
    p53 = res.mean_p53
    peaks = [i for i in range(1, len(p53) - 1)
             if p53[i] > p53[i - 1] and p53[i] >= p53[i + 1] and p53[i] > 0.12]
    assert len(peaks) >= 2, f"expected pulses under sustained damage, got {len(peaks)}"
    period = float(np.mean(np.diff(res.t_h[peaks])))
    assert 3.5 <= period <= 8.0, f"pulse period {period:.2f} h outside the published window"


def test_undamaged_cells_show_no_pulsing():
    """Specificity: the oscillator must be stress-gated, not constitutive."""
    res = simulate(A549, None, 0.0, t_end_h=48, n_cells=16, dt_h=0.02,
                   k_cyc=K_A549, record_every_h=0.1)
    p53 = res.mean_p53
    assert np.nanmax(p53) < 0.12, f"basal p53 rose to {np.nanmax(p53):.3f} without damage"


def test_damage_executes_death_not_just_commitment():
    """The defect that made the C++ prototype unusable: MOMP saturated but
    caspase-3 stayed clamped at 0.04, so no cell ever died. Here heavy
    damage must drive caspase-3 past commitment and remove cells."""
    res = simulate(A549, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02,
                   k_cyc=K_A549, forced_D=0.60)
    assert res.deaths[-1, 0] > 0, "no cell died under heavy sustained damage"
    assert res.final.Y[:, IX["C3"]].max() > 0.5, "caspase-3 never crossed commitment"


def test_mutant_p53_resists_damage_induced_death():
    """TP53 loss-of-function removes transactivation, so the PUMA-driven
    intrinsic route is lost and the same damage kills far fewer cells."""
    wt = simulate(A549, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02,
                  k_cyc=K_A549, forced_D=0.60)
    mut_line = get_line("HT-29")
    mut = simulate(mut_line, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02,
                   k_cyc=calibrate_cycle_scale(mut_line), forced_D=0.60)
    assert mut.deaths[-1, 0] < 0.5 * wt.deaths[-1, 0], (
        f"p53-mutant deaths {mut.deaths[-1, 0]:.1f} not lower than wild-type "
        f"{wt.deaths[-1, 0]:.1f}")


def test_dose_response_is_monotone_and_yields_a_finite_ic50():
    """A cytotoxic dose-response must fall with dose and cross 50 %."""
    conc = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
    viab, _ = dose_response(A549, get_drug("cisplatin"), conc,
                            n_cells_per_conc=24, k_cyc=K_A549)
    assert abs(viab[0] - 1.0) < 1e-9, "control viability must be 1 by construction"
    treated = viab[1:]
    assert np.all(np.diff(treated) <= 1e-9), f"viability not monotone: {treated}"
    ic50 = ic50_from_curve(conc, treated)
    assert np.isfinite(ic50) and ic50 > 0, f"IC50 not finite: {ic50}"


def test_simulation_is_deterministic_for_a_fixed_seed():
    a = simulate(A549, None, 0.0, t_end_h=24, n_cells=32, dt_h=0.02, k_cyc=K_A549, seed=7)
    b = simulate(A549, None, 0.0, t_end_h=24, n_cells=32, dt_h=0.02, k_cyc=K_A549, seed=7)
    assert np.allclose(a.alive_weight, b.alive_weight)


def test_state_stays_bounded_and_finite():
    """No NaN or runaway: every variable stays inside its clamp."""
    res = simulate(A549, get_drug("paclitaxel"), 0.1, t_end_h=72, n_cells=32,
                   dt_h=0.02, k_cyc=K_A549)
    Y = res.final.Y
    assert np.all(np.isfinite(Y)), "non-finite state"
    assert np.all(Y >= -1e-9), "negative concentration"


def test_cycle_scale_tracks_the_doubling_time():
    """A slower line must get a smaller cycle rate scale."""
    fast, slow = get_line("HCT116"), get_line("MCF7")   # 18 h vs 29 h
    assert calibrate_cycle_scale(fast) > calibrate_cycle_scale(slow)


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
    print(f"{len(fns) - failed}/{len(fns)} engine gates pass")
    raise SystemExit(1 if failed else 0)
