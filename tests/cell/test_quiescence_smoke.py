#!/usr/bin/env python3
"""Phase-5 gates: the proliferation / quiescence split at mitotic exit.

A daughter does not always start the next cycle. CDK2 activity bifurcates
as the mother divides and the low branch sits in G0 until it is recruited
back (Spencer et al. 2013 Cell 155:369). The engine had no such state:
every cell cycled, which is the documented miss against the Cell Tracking
Challenge HeLa movie, where 40–60 % of tracked cells never divide inside
46 h and the Kaplan-Meier curve plateaus while the engine's keeps falling.

The gates below are about the mechanism being what it claims to be —
off by default, a genuine hold rather than a slow-down, reversible when
recruitment is switched on, and correctly halving population growth when
half of each generation stops. What it buys against the real movie is
measured in `scripts/validate_dish.py`, not asserted here.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, simulate  # noqa: E402
from cellsim.cell.library import get_line  # noqa: E402

LINE = get_line("A549")
HOURS = 72.0


def _growth(p: Params, seed: int = 1, n: int = 64) -> float:
    res = simulate(LINE, None, 0.0, t_end_h=HOURS, n_cells=n, p=p, seed=seed,
                   record_every_h=HOURS)
    return float(res.alive_weight[-1, 0] / res.alive_weight[0, 0])


def test_quiescence_is_off_by_default():
    """The population must carry no G0 state at all unless asked, so that
    every result predating this mechanism is untouched."""
    res = simulate(LINE, None, 0.0, t_end_h=12.0, n_cells=8, record_every_h=12.0)
    assert res.final.quiescent is None
    withq = simulate(LINE, None, 0.0, t_end_h=12.0, n_cells=8, record_every_h=12.0,
                     p=Params(quiescent_fraction=0.5))
    assert withq.final.quiescent is not None


def test_a_quiescent_population_grows_more_slowly():
    """Half of each generation dropping out must show up as slower
    growth, and all of it dropping out must stop growth after the cells
    already in cycle have finished."""
    full = _growth(Params())
    half = _growth(Params(quiescent_fraction=0.5))
    none_ = _growth(Params(quiescent_fraction=1.0))
    assert half < full * 0.75, f"half quiescent: x{half:.1f} vs x{full:.1f}"
    assert none_ < half, f"fully quiescent: x{none_:.1f} vs half x{half:.1f}"
    # with nothing re-entering, growth stops once the cells that were
    # already cycling have divided once
    assert none_ < 4.0, f"a fully quiescent population should stall, got x{none_:.1f}"


def test_recruitment_out_of_g0_restores_growth():
    """Quiescence must be a reversible hold, not a death sentence: with a
    recruitment rate, the same population grows again."""
    stuck = _growth(Params(quiescent_fraction=0.8))
    recruited = _growth(Params(quiescent_fraction=0.8, quiescence_exit_per_h=0.1))
    assert recruited > stuck * 1.3, (
        f"recruitment should restore growth: x{recruited:.1f} vs stuck x{stuck:.1f}")


def test_g0_holds_the_cycle_rather_than_slowing_it():
    """A quiescent cell keeps its cycle state; it is not merely running
    slowly. With everything in G0 and no recruitment, the cyclins of the
    cells that have arrived there must stop moving."""
    p = Params(quiescent_fraction=1.0)
    a = simulate(LINE, None, 0.0, t_end_h=48.0, n_cells=48, p=p, seed=1,
                 record_every_h=48.0)
    pop = a.final
    q = pop.quiescent & pop.alive
    assert q.sum() >= 10, f"expected a G0 pool, got {int(q.sum())}"
    # those cells sit in early G1, where fresh_cycle_state leaves them
    from cellsim.cell.engine import IX
    cyc_b = pop.Y[q, IX["CycB"]]
    assert float(cyc_b.max()) < 0.25, (
        f"a held cell must not accumulate cyclin B toward mitosis: {cyc_b.max():.3f}")


def test_quiescence_does_not_kill():
    """Specificity: cells in G0 are not dying, they are waiting."""
    res = simulate(LINE, None, 0.0, t_end_h=72.0, n_cells=48,
                   p=Params(quiescent_fraction=1.0), record_every_h=72.0)
    pop = res.final
    assert pop.alive.all(), f"{int((~pop.alive).sum())} cells died in G0 with no drug"
    assert float(res.deaths[-1, 0]) == 0.0


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
        except Exception as e:                      # noqa: BLE001
            failed += 1
            print(f"  ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} quiescence gates pass")
    raise SystemExit(1 if failed else 0)
