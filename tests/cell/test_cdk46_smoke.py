#!/usr/bin/env python3
"""Gates for the CDK4/6 inhibitor, where Rb status IS the mechanism.

Palbociclib blocks cyclin D-CDK4/6, so Rb stays unphosphorylated and the
cell stops in G1 — and where there is no working Rb there is nothing to
stop, and the drug does nothing (Fry et al. 2004 Mol Cancer Ther 3:1427).
That makes it the second drug, after nutlin, whose line-to-line
differences come from mechanism rather than a fitted marker. Two library
lines lack Rb for different reasons: MDA-MB-468 has lost RB1, and HeLa's
HPV18 E7 holds Rb off. Both were resistant in GDSC; neither was used in
any fit.

Checked here:
* the Rb-less lines are untouched at any dose, and an RB1 knockout makes
  a sensitive line behave the same way (two routes, one answer);
* in a line with Rb, cells pile up in G1 — an arrest, not a kill;
* it is cytostatic: almost nobody dies, and p53 does not move.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _tiers import run  # noqa: E402
from cellsim.cell.engine import calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.perturb import Perturbation  # noqa: E402
from cellsim.cell.trace import trace_cells  # noqa: E402

PALBO = get_drug("palbociclib")


def _relative(line_name, conc, perturbations=None, hours=72.0, seed=2):
    line = get_line(line_name)
    k = calibrate_cycle_scale(line)
    kw = dict(t_end_h=hours, n_cells=64, k_cyc=k, seed=seed, perturbations=perturbations)
    ctl = simulate(line, None, 0.0, **kw)
    trt = simulate(line, PALBO, conc, **kw)
    return float(trt.alive_weight[-1, 0] / ctl.alive_weight[-1, 0]), trt


def test_rb_less_lines_are_untouched():
    for name in ("HeLa", "MDA-MB-468"):
        assert not get_line(name).rb_functional
        for c in (1.0, 30.0):
            rel, _ = _relative(name, c)
            assert rel > 0.95, f"{name} has no working Rb and should ignore palbociclib: {rel:.2f}"


def test_rb1_knockout_gives_the_same_resistance():
    with_rb, _ = _relative("A549", 10.0)
    without, _ = _relative("A549", 10.0, perturbations=[Perturbation("RB1")])
    assert with_rb < 0.6, f"A549 with Rb should respond: {with_rb:.2f}"
    assert without > 0.9, f"knocking out RB1 should remove the response: {without:.2f}"


def test_cells_stop_after_one_division():
    """Cells already past the restriction point finish their cycle — some
    of those divisions land as late as ~20 h — and then stop. So the test
    is on SECOND divisions: in 40 h an untreated A549 culture has many,
    a palbociclib-treated one almost none. (A G1 fraction is a poor
    measure: the engine's phase rule already calls ~77 % of untreated A549
    cells G1.)"""
    def second_divisions(drug, conc):
        tr = trace_cells("A549", drug, conc, hours=40, n_cells=40, every_h=1.0, seed=3)
        return sum(1 for c in tr.cells for e in c["events"]
                   if e["kind"] == "divided" and "generation 1)" not in e["text"])
    drug, ctl = second_divisions(PALBO, 10.0), second_divisions(None, 0.0)
    assert ctl >= 10 and drug <= 0.1 * ctl, f"second divisions in 40 h: drug {drug}, control {ctl}"


def test_cytostatic_not_cytotoxic():
    _, trt = _relative("A549", 10.0)
    dead = float(trt.deaths[-1, 0] / (trt.deaths[-1, 0] + trt.alive_weight[-1, 0]))
    assert dead < 0.05, f"palbociclib should arrest, not kill: {dead:.1%} dead"
    tr = trace_cells("A549", PALBO, 10.0, hours=24, n_cells=12, every_h=1.0, seed=4)
    p53 = np.array([c["series"]["p53"] for c in tr.cells])
    assert p53.max() < 1.6 * p53[:, 0].mean(), "a CDK4/6 inhibitor does not damage DNA"


if __name__ == "__main__":
    raise SystemExit(run(globals(), "CDK4/6"))
