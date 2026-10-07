#!/usr/bin/env python3
"""The notebook API returns tidy tables, and the right numbers in them.

`cellsim.api` is what a biologist actually touches, so it is gated on
two things: that every function returns one row per observation with the
columns it promises (a DataFrame when pandas is installed, a list of
dicts otherwise, with the same keys either way), and that the numbers
agree with the engine underneath rather than being independently
re-derived here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim import api  # noqa: E402
from cellsim.cell.engine import Params, calibrate_cycle_scale, ic50  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402


def _rows(table) -> list[dict]:
    """Normalise a DataFrame or a list of dicts to a list of dicts."""
    return table.to_dict("records") if hasattr(table, "to_dict") else table


def test_catalogues_list_every_line_and_drug():
    from cellsim.cell.library import CELL_LINES, DRUGS
    lines, drugs = _rows(api.lines()), _rows(api.drugs())
    assert len(lines) == len(CELL_LINES) and len(drugs) == len(DRUGS)
    assert {r["line"] for r in lines} == set(CELL_LINES)
    for r in drugs:
        assert r["mechanism"] and r["source"], r
        if r["fitted_constant"]:
            assert r["fitted_value"] is not None and r["fitted_value"] > 0, r


def test_curve_is_tidy_monotone_and_agrees_with_the_engine():
    conc = [1.0, 3.0, 10.0, 30.0, 100.0]
    rows = _rows(api.curve("A549", "cisplatin", concentrations=conc, hours=48, n_cells=24))
    assert len(rows) == len(conc), "one row per concentration"
    assert [r["conc_uM"] for r in rows] == conc
    for r in rows:
        assert r["line"] == "A549" and r["drug"] == "cisplatin"
        assert 0.0 <= r["viability"] <= 1.2, r
    viab = [r["viability"] for r in rows]
    assert viab[0] > viab[-1] + 0.3, f"curve should fall with dose: {viab}"
    # the IC50 column must be the curve's own crossing, not a separate run
    ic = rows[0]["ic50_uM"]
    assert all(abs(r["ic50_uM"] - ic) < 1e-9 for r in rows)
    below = [r["conc_uM"] for r in rows if r["viability"] < 0.5]
    if below:
        assert ic <= min(below) + 1e-9, (ic, below)


def test_a_three_hour_exposure_needs_an_unusable_concentration():
    """A spindle poison can only act on cells that attempt mitosis while
    it is present, so a 3 h window is nearly useless: Georgiadis 1997 put
    its half-kill above 32 µM, the top of their tested range, against
    0.027 µM for a continuous exposure.

    This used to assert the 3 h value was INFINITE. It is finite once
    drug retained after wash-out is modelled — the engine says ~150 µM —
    and that is the better answer: the measurement says "above 32", not
    "unreachable". What the gate checks now is that 3 h needs a
    concentration far beyond anything usable, and orders of magnitude
    above the continuous one."""
    rows = _rows(api.exposure("A549", "paclitaxel", hours=[3, 72], n_cells=24))
    short, long = rows[0], rows[1]
    assert long["reaches_half_kill"] is True and long["c50_uM"] > 0
    if short["reaches_half_kill"]:
        assert short["c50_uM"] > 32.0, (
            f"3 h half-kill at {short['c50_uM']:.3g} µM; Georgiadis saw none "
            f"below 32 µM")
        assert short["c50_uM"] > 100 * long["c50_uM"], (
            f"3 h should need orders of magnitude more than continuous: "
            f"{short['c50_uM']:.3g} vs {long['c50_uM']:.3g} µM")


def test_gr_curve_puts_a_simulation_on_the_same_axis_as_a_plate():
    """GR = 1 is untreated growth, 0 is complete cytostasis, below 0 is net
    cell loss. A measured plate and a simulation both reduced to GR can be
    compared; viability cannot be compared that way between lines, which
    is why GR exists."""
    rows = _rows(api.gr_curve("A549", "cisplatin", concentrations=[1, 3, 10, 30, 100],
                              n_cells=24))
    assert len(rows) == 5
    gr = [r["gr"] for r in rows]
    assert gr[0] > 0.8, f"a low dose should barely touch growth: {gr}"
    assert gr[-1] < 0.1, f"a high dose should at least arrest growth: {gr}"
    assert gr[0] >= gr[-1] - 0.05, f"GR should fall with dose: {gr}"
    # the control's own growth is what makes GR comparable; A549 doubles
    # every 22 h, so a 72 h assay is a little over three doublings
    doublings = rows[0]["control_doublings"]
    assert 2.5 < doublings < 4.0, doublings
    assert all(r["gr"] >= -1.0 for r in rows), "GR is floored at total loss"


def test_ic50_spread_exposes_the_engines_own_stochasticity():
    """The number a user should see before quoting an IC50. The engine
    draws a finite sample of lineages, so its IC50 moves 10-20 % between
    seeds; with heterogeneity off it is near-deterministic, which is what
    identifies the cause as the draw rather than cell count."""
    rows = _rows(api.ic50_spread("A549", "cisplatin", seeds=(1, 2, 3), n_cells=24))
    assert len(rows) == 3 and {r["seed"] for r in rows} == {1, 2, 3}
    vals = [r["ic50_uM"] for r in rows]
    assert all(np.isfinite(v) and v > 0 for v in vals), vals
    spread = rows[0]["spread_fold"]
    assert abs(spread - max(vals) / min(vals)) < 1e-9
    assert 1.0 < spread < 2.0, f"expected a real but bounded spread, got {spread}"
    # with identical cells the same measurement is near-deterministic
    flat = _rows(api.ic50_spread("A549", "cisplatin", seeds=(1, 2, 3), n_cells=24,
                                 params=Params(het_sigma=0.0)))
    assert flat[0]["spread_fold"] < 1.05, (
        f"heterogeneity should be the source of the spread: {flat[0]['spread_fold']}")
    assert flat[0]["spread_fold"] < spread


def test_washout_shows_the_drug_leaving_and_the_colony_recovering():
    rows = _rows(api.washout("A549", "cisplatin", conc_uM=20.0, exposure_h=12.0,
                             hours=96.0, n_cells=48, record_every_h=12.0))
    assert len(rows) == 9, len(rows)
    # Times must come back clean, not as the integrator's accumulated
    # clock (11.999999999999833), or grouping and plotting by time break.
    assert [r["t_h"] for r in rows] == [0, 12, 24, 36, 48, 60, 72, 84, 96], \
        [r["t_h"] for r in rows]
    assert rows[0]["drug_present"] and not rows[-1]["drug_present"]
    assert all(r["control_fold_change"] >= 1.0 for r in rows), "control must not shrink"
    # the colony is held back while the drug is on, then grows again
    trough = min(r["fold_change"] for r in rows)
    assert rows[-1]["fold_change"] > trough, "no recovery after wash-out"
    assert rows[-1]["relative_to_control"] < 1.0, "treatment should leave fewer cells"


def test_combination_reports_each_arm_against_independent_action():
    cis, pac = get_drug("cisplatin"), get_drug("paclitaxel")
    line = get_line("A549")
    k = calibrate_cycle_scale(line, Params())
    ic_c = ic50(line, cis, guess_uM=10.0, k_cyc=k, n_cells_per_conc=16, n_seeds=1)
    ic_p = ic50(line, pac, guess_uM=0.03, k_cyc=k, n_cells_per_conc=16, n_seeds=1)
    # 72 h readout, which is the condition docs/VALIDATION.md measured the
    # antagonism under; it weakens as the readout lengthens and survivors
    # resume growing.
    rows = _rows(api.combination("A549", "cisplatin", "paclitaxel", ic_c, ic_p,
                                 readout_h=72.0, n_cells=64, seeds=(1, 2)))
    assert len(rows) == 5, [r["arm"] for r in rows]
    alone = [r for r in rows if r["arm"].endswith("alone")]
    together = [r for r in rows if not r["arm"].endswith("alone")]
    assert len(alone) == 2 and len(together) == 3
    for r in alone:
        # None in the dict form; pandas renders the same absence as NaN
        v = r["ratio_to_independent"]
        assert v is None or (isinstance(v, float) and np.isnan(v)), v
    for r in together:
        assert r["ratio_to_independent"] > 0 and r["n_seeds"] == 2
        assert r["surviving_fraction"] < min(a["surviving_fraction"] for a in alone) + 0.05, (
            "the pair should not leave more cells than either drug alone")
    # The documented phenotype: a cytostatic antagonises a mitosis-specific
    # partner. Gated on the mean of the three arrangements, because this
    # test runs far fewer cells and seeds than the experiment that
    # quantifies it (scripts/experiment_combination.py).
    ratios = [r["ratio_to_independent"] for r in together]
    assert float(np.mean(ratios)) > 1.05, ratios


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
    print(f"{len(fns) - failed}/{len(fns)} api gates pass")
    raise SystemExit(1 if failed else 0)
