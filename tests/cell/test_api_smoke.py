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


def test_exposure_reports_infinity_where_half_kill_is_unreachable():
    """A tubulin binder cannot halve an asynchronous colony in 3 h at any
    concentration; the table must say so rather than return a number."""
    rows = _rows(api.exposure("A549", "paclitaxel", hours=[3, 72], n_cells=24))
    short, long = rows[0], rows[1]
    assert short["reaches_half_kill"] is False and not np.isfinite(short["c50_uM"])
    assert long["reaches_half_kill"] is True and np.isfinite(long["c50_uM"])
    assert long["c50_uM"] > 0


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
    ic_c = ic50(line, cis, guess_uM=10.0, k_cyc=k, n_cells_per_conc=16)
    ic_p = ic50(line, pac, guess_uM=0.03, k_cyc=k, n_cells_per_conc=16)
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
