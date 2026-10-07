#!/usr/bin/env python3
"""Phase-3 gates: resistance that EVOLVES, not only resistance selected.

Everything before this could select from the spread a population already
had. Heterogeneity was drawn once, so every lineage that survived had
been resistant from the start. `DishParams.mutation_rate` lets a daughter
differ heritably from her mother, which is a different claim and needs
its own control.

The control that makes the claim testable is a HOMOGENEOUS start
(`Params.het_sigma = 0`): identical cells, so selection has nothing to
act on. Any resistance that appears there arose during treatment.

Kept small (short course, small lattice) so it runs in CI; the full
sweep is `scripts/experiment_evolution.py`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import DishParams, run_dish  # noqa: E402
from cellsim.cell.engine import Params, calibrate_cycle_scale, ic50  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

LINE = get_line("A549")
DRUG = get_drug("cisplatin")
K = calibrate_cycle_scale(LINE, Params())
DOSE_UM = ic50(LINE, DRUG, guess_uM=10.0, n_cells_per_conc=24, k_cyc=K, n_seeds=1)
DAYS = 5.0
MUT_RATE, MUT_SD = 0.05, 0.3

_CACHE: dict = {}


def _run(mutation: bool, treated: bool, het: float = 0.0) -> dict:
    key = (mutation, treated, het)
    if key not in _CACHE:
        p = Params(het_sigma=het, cycle_cv=0.25)
        dp = DishParams(geometry="monolayer",
                        mutation_rate=MUT_RATE if mutation else 0.0,
                        mutation_effect_sd=MUT_SD if mutation else 0.0)
        dose = (lambda t: DOSE_UM) if treated else (lambda t: 0.0)
        res = run_dish(LINE, DRUG, dose, t_end_h=DAYS * 24.0, n_seed=120, grid_sites=48,
                       record_every_h=24.0, seed=1, p=p, dp=dp)
        r, d = res.record, res.dish
        _CACHE[key] = {"live": r.n_live[-1], "uptake": r.mean_uptake[-1],
                       "reserve": r.mean_reserve[-1], "mutations": r.n_mutations[-1],
                       "bcl2": d.cells.het_bcl2.copy(), "het_uptake": d.cells.het_uptake.copy()}
    return _CACHE[key]


def test_identical_cells_without_mutation_cannot_become_resistant():
    """The control. With no spread to select from and no heritable change,
    every surviving cell must still be exactly the cell we started with."""
    out = _run(mutation=False, treated=True)
    assert out["mutations"] == 0, out["mutations"]
    if out["live"] == 0:
        return                                   # eradicated: also no resistance
    assert np.allclose(out["bcl2"], 1.0) and np.allclose(out["het_uptake"], 1.0), (
        f"traits drifted without mutation: reserve {out['reserve']:.3f}, "
        f"uptake {out['uptake']:.3f}")


def test_mutation_is_neutral_without_the_drug():
    """Specificity: mutation alone must not look like resistance. Untreated,
    the mean traits stay at 1 even though changes are accumulating."""
    out = _run(mutation=True, treated=False)
    assert out["mutations"] > 0, "no mutations occurred; the rate is too low to test"
    assert abs(out["uptake"] - 1.0) < 0.12, out["uptake"]
    assert abs(out["reserve"] - 1.0) < 0.12, out["reserve"]


def test_resistance_arises_from_an_identical_population_under_treatment():
    """The claim. From identical cells, treatment plus heritable change
    leaves survivors that accumulate markedly less drug — resistance that
    was not there to be selected."""
    treated = _run(mutation=True, treated=True)
    control = _run(mutation=True, treated=False)
    assert treated["live"] > 0, "nothing survived; nothing to measure"
    assert treated["uptake"] < 0.85 * control["uptake"], (
        f"survivors not enriched for low accumulation: {treated['uptake']:.3f} "
        f"vs untreated {control['uptake']:.3f}")


def test_evolved_resistance_stays_in_the_range_derived_lines_reach():
    """A model that evolved resistance without limit would be useless.
    Derived lines reach two- to eight-fold (McDermott 2014); the engine
    must not exceed that from a few days of selection."""
    treated = _run(mutation=True, treated=True)
    control = _run(mutation=True, treated=False)
    if treated["live"] == 0:
        return
    fold = control["uptake"] / max(treated["uptake"], 1e-9)
    assert fold < 8.0, f"implausible accumulation shift of {fold:.1f}x in {DAYS:g} days"


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
    print(f"{len(fns) - failed}/{len(fns)} evolution gates pass")
    raise SystemExit(1 if failed else 0)
