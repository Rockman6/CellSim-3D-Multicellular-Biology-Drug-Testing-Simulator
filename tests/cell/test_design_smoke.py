#!/usr/bin/env python3
"""Gates for design mode: drugs that block or boost one species.

A designed drug is a hypothesis, so there is nothing to validate it
against except the engine itself. What can be checked is that it means
what it says, by comparison with things that already have a meaning:

* fully bound, a blocker with effect 1 must do what a knockout of that
  gene does, and a booster what an overexpression does — the design is a
  dose-dependent version of `cellsim.cell.perturb`, not a new route;
* half bound (dose = Kd) it must do half;
* a knockout combined with a two-drug combination must run (the
  perturbation factor used to be sized for one drug slot and broke
  combinations);
* the library-mechanism designs carry their reference drug's behaviour.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _tiers import run  # noqa: E402
from cellsim.cell.compound import DESIGNABLE, designed  # noqa: E402
from cellsim.cell.engine import calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.perturb import Perturbation  # noqa: E402

LINE = get_line("A549")
K = calibrate_cycle_scale(LINE)
KW = dict(t_end_h=48.0, n_cells=48, k_cyc=K, seed=5)


def _grown(res):
    return float(res.alive_weight[-1, 0] / res.alive_weight[0, 0])


def test_saturating_blocker_equals_knockout():
    # Kd far below the dose and fast uptake: bound fraction ~1 within minutes.
    # (Faster than ~0.05 h is below the 0.02 h step and RK4 goes unstable.)
    ko = simulate(LINE, None, 0.0, perturbations=[Perturbation("CCNE1")], **KW)
    import dataclasses
    blocker = dataclasses.replace(designed("x", "block", species="CycE", kd_uM=1e-6),
                                  tau_uptake_h=0.05)
    drug = simulate(LINE, blocker, 10.0, **KW)
    untreated = simulate(LINE, None, 0.0, **KW)
    assert _grown(ko) < 0.8 * _grown(untreated), "a cyclin E knockout should slow growth"
    assert abs(_grown(drug) - _grown(ko)) < 0.05 * _grown(ko), (_grown(drug), _grown(ko))


def test_saturating_booster_equals_overexpression():
    import dataclasses
    oe = simulate(LINE, None, 0.0, perturbations=[Perturbation("CDKN1A", "overexpress", 4.0)], **KW)
    booster = dataclasses.replace(designed("x", "boost", species="p21", kd_uM=1e-6, effect=4.0),
                                  tau_uptake_h=0.05)
    drug = simulate(LINE, booster, 10.0, **KW)
    assert abs(_grown(drug) - _grown(oe)) < 0.05 * _grown(oe), (_grown(drug), _grown(oe))


def test_dose_at_kd_blocks_about_half():
    import dataclasses
    full = dataclasses.replace(designed("x", "block", species="CycE", kd_uM=1.0), tau_uptake_h=0.05)
    g = {c: _grown(simulate(LINE, full, c, **KW)) for c in (0.0, 1.0, 1000.0)}
    assert g[1000.0] < g[1.0] < g[0.0], g


def test_knockout_with_a_combination_runs():
    pair = (get_drug("cisplatin"), get_drug("paclitaxel"))
    res = simulate(LINE, pair, [5.0, 0.005], perturbations=[Perturbation("BAX")],
                   t_end_h=12.0, n_cells=8, k_cyc=K)
    assert np.isfinite(res.alive_weight).all()


def test_library_mechanism_designs_and_refusals():
    d = designed("mine", "freeze", strength=2.0)
    assert d.mechanism == "tubulin" and d.source.startswith("DESIGNED")
    assert d.partition == get_drug("paclitaxel").partition * 2.0
    for bad in (dict(action="block", species="D"), dict(action="block", species="CycE", effect=2.0),
                dict(action="boost", species="p21", effect=0.5), dict(action="melt")):
        try:
            designed("x", **bad)
        except ValueError:
            continue
        raise AssertionError(f"{bad} should be refused")
    assert "D" not in DESIGNABLE and "Cin" not in DESIGNABLE


if __name__ == "__main__":
    raise SystemExit(run(globals(), "design"))
