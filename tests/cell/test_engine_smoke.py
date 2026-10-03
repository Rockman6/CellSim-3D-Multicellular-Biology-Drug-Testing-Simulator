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


def test_mutant_p53_resists_damage_induced_death_when_only_p53_can_kill():
    """With the p53-independent route OFF, TP53 loss of function removes
    transactivation, so the PUMA-driven intrinsic route is lost and the
    same damage kills far fewer cells. This pins the p53 axis wiring."""
    p53_only = Params(p73_gain=0.0)
    wt = simulate(A549, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02,
                  k_cyc=K_A549, forced_D=0.60, p=p53_only)
    mut_line = get_line("HT-29")
    mut = simulate(mut_line, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02,
                   k_cyc=calibrate_cycle_scale(mut_line, p53_only), forced_D=0.60, p=p53_only)
    assert mut.deaths[-1, 0] < 0.5 * wt.deaths[-1, 0], (
        f"p53-mutant deaths {mut.deaths[-1, 0]:.1f} not lower than wild-type "
        f"{wt.deaths[-1, 0]:.1f}")


def test_p53_mutant_keeps_partial_resistance_under_sustained_damage():
    """With the defaults (p73 route on, anti-apoptotic buffer on), a
    TP53-mutant line should be substantially sensitised relative to the
    p53-only model but NOT fully: real TP53-mutant cells are more
    resistant to DNA damage, not immune to it.

    Before the Bcl-2 reserve was modelled as a stoichiometric buffer the
    route removed p53 protection entirely (mutant survival 0.000), which
    is stronger than the literature supports. Pinned here so the balance
    cannot drift unnoticed in either direction.
    """
    mut = get_line("HT-29")
    on, off = Params(), Params(p73_gain=0.0)
    s_on = _surviving_fraction(mut, on, 0.30)
    s_off = _surviving_fraction(mut, off, 0.30)
    wt_on = _surviving_fraction(A549, on, 0.30)
    assert s_on < 0.5 * s_off, (
        f"p73 route barely sensitised the mutant: {s_off:.3f} -> {s_on:.3f}")
    assert s_on > wt_on, (
        f"mutant ({s_on:.3f}) should keep some advantage over wild-type "
        f"({wt_on:.3f}); total loss of p53 protection is unphysiological")


def _surviving_fraction(line, params, damage):
    """Alive weight after 72 h of sustained damage, relative to this
    line's own untreated control: what a viability assay reads."""
    k = calibrate_cycle_scale(line, params)
    ctl = simulate(line, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02, k_cyc=k, p=params)
    dmg = simulate(line, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02, k_cyc=k,
                   p=params, forced_D=damage)
    return float(dmg.alive_weight[-1, 0] / max(ctl.alive_weight[-1, 0], 1e-9))


def test_dose_response_is_monotone_and_yields_a_finite_ic50():
    """A cytotoxic dose-response must fall with dose and cross 50 %."""
    conc = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
    viab, _ = dose_response(A549, get_drug("cisplatin"), conc,
                            n_cells_per_conc=24, k_cyc=K_A549)
    assert abs(viab[0] - 1.0) < 1e-9, "control viability must be 1 by construction"
    treated = viab[1:]
    # 24 representative cells per dose give a few percent of sampling noise
    # below the kill threshold; monotone up to that noise, and a real fall.
    assert np.all(np.diff(treated) <= 0.15), f"viability rises with dose: {treated}"
    assert treated[0] > 0.8 and treated[-1] < 0.2, f"no real fall across the range: {treated}"
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


def test_zero_variability_gives_identical_cells():
    from cellsim.cell.engine import make_population
    pop = make_population(64, 1.0, np.random.default_rng(3), het_sigma=0.0)
    assert np.all(pop.het_bcl2 == 1.0) and np.all(pop.het_uptake == 1.0)
    pop = make_population(4096, 1.0, np.random.default_rng(3), het_sigma=0.3)
    assert abs(pop.het_bcl2.mean() - 1.0) < 0.03, "multipliers must average 1"
    assert 0.2 < np.log(pop.het_bcl2).std() < 0.4


def _crossing(c, v, level):
    for i in range(1, len(c)):
        if v[i] <= level < v[i - 1]:
            f = (v[i - 1] - level) / max(v[i - 1] - v[i], 1e-12)
            return float(np.exp(np.log(c[i - 1]) + f * (np.log(c[i]) - np.log(c[i - 1]))))
    return float("nan")


def test_variability_widens_the_kill_transition():
    """Fractional killing (Spencer 2009): identical cells flip from alive
    to dead within a ~1.2x dose window; protein-level variability widens
    it (measured ~2x at sigma 0.4). Still far steeper than a Hill slope
    of 1 (~32x) — see docs/VALIDATION.md, curve shape is an open item."""
    conc = np.logspace(-2, 0, 21)                     # 10 points per decade
    import dataclasses
    drug = dataclasses.replace(get_drug("doxorubicin"), k_damage_per_uM_h=0.0087)

    def width(sigma):
        viab, _ = dose_response(A549, drug, conc, n_cells_per_conc=48, k_cyc=K_A549,
                                p=Params(het_sigma=sigma))
        return _crossing(conc, viab[1:], 0.15) / _crossing(conc, viab[1:], 0.85)

    w0, w4 = width(0.0), width(0.4)
    assert w4 > 1.4 * w0, f"85%->15% width {w0:.2f}x -> {w4:.2f}x; variability did not widen it"


def test_p73_route_kills_p53_mutant_cells_but_spares_untreated_ones():
    """ATM/c-Abl -> p73 -> PUMA (Gong 1999, Agami 1999) is p53-independent:
    switching it on kills a TP53-mutant line that the p53-only model
    leaves alive, and untreated cells are untouched because the route
    only engages once ATM is genuinely activated."""
    ht29 = get_line("HT-29")
    off, on = Params(p73_gain=0.0), Params(p73_gain=1.0)
    k = calibrate_cycle_scale(ht29, off)
    off_dmg = simulate(ht29, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02, k_cyc=k,
                       forced_D=0.6, p=off)
    on_dmg = simulate(ht29, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02, k_cyc=k,
                      forced_D=0.6, p=on)
    on_ctl = simulate(ht29, None, 0.0, t_end_h=72, n_cells=32, dt_h=0.02, k_cyc=k, p=on)
    assert off_dmg.deaths[-1, 0] == 0, "p53-only model should leave HT-29 alive"
    assert on_dmg.deaths[-1, 0] > 0, "p73 route did not kill damaged p53-mutant cells"
    assert on_ctl.deaths[-1, 0] == 0, "p73 route killed untreated cells"


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
