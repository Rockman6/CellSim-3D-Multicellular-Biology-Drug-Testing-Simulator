#!/usr/bin/env python3
"""Phase-2 gates: the engine must get schedule dependence right.

These are the questions mechanism answers and marker correlations do
not, and they are asked entirely *within* one cell line, so none of
them runs into the per-line prediction ceiling recorded in
docs/VALIDATION.md.

The phenotype under test is the documented contrast between a tubulin
binder and a DNA-damage agent:

* paclitaxel kills only cells that attempt mitosis while it is present,
  and only above a threshold tubulin occupancy, so its effect is
  strongly limited by exposure TIME (Lopes 1993 Cancer Chemother
  Pharmacol 32:235; Gascoigne & Taylor 2008 Cancer Cell 14:111). Its
  clinical PK/PD driver is time above a threshold concentration rather
  than AUC or peak (Gianni 1995 J Clin Oncol 13:180; Huizing 1993
  J Clin Oncol 11:2127).
* cisplatin forms adducts as soon as it is inside the cell, so it is
  comparatively concentration-driven.

Kept small (few windows, few cells) so it runs in CI; the full sweep is
`scripts/experiment_schedule.py`.
"""
from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import calibrate_cycle_scale, ic50  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.stream import Schedule, run_stream  # noqa: E402

LINE = get_line("A549")
K = calibrate_cycle_scale(LINE)
N_CELLS = 48

# These phenotypes belong to the mechanism, not to the fitted potency, so
# each drug is held at the potency its doses were designed for (the fit
# before the October 2026 refit). Cisplatin and paclitaxel enter the
# engine only as potency x concentration, so this is exactly a rescaling
# of the doses, and refitting the library cannot carry a dose across a
# threshold and break a test that is about shape.
DESIGN_POTENCY = {"cisplatin": {"k_damage_per_uM_h": 0.000951},
                  "doxorubicin": {"k_damage_per_uM_h": 0.00794},
                  "paclitaxel": {"partition": 0.153}}


def _drug(name: str):
    return dataclasses.replace(get_drug(name), **DESIGN_POTENCY[name])


def _surviving(drug, conc_uM: float, hours: float, seed: int = 1) -> float:
    treated = run_stream(LINE, drug, Schedule.parse(f"0:{conc_uM:g},{hours:g}:0"),
                         t_end_h=72.0, n_cells=N_CELLS, record_every_h=72.0, seed=seed)
    control = run_stream(LINE, None, Schedule.constant(0.0),
                         t_end_h=72.0, n_cells=N_CELLS, record_every_h=72.0, seed=seed)
    return (treated[-1]["population"]["alive_weight"]
            / max(control[-1]["population"]["alive_weight"], 1e-12))


def test_cisplatin_short_exposures_follow_concentration_times_time():
    """Ozawa et al. (Cancer Res 1989 49:3823) measured cisplatin's kill as
    governed by C x T: the concentration needed for an iso-effect halves
    when exposure doubles (slope -1 on log-log). Over 1-4 h, before
    repair has time to matter, the engine must reproduce that law."""
    cis = _drug("cisplatin")
    c1 = ic50(LINE, cis, guess_uM=400.0, exposure_h=1.0, n_cells_per_conc=16, k_cyc=K)
    c4 = ic50(LINE, cis, guess_uM=100.0, exposure_h=4.0, n_cells_per_conc=16, k_cyc=K)
    slope = math.log(c4 / c1) / math.log(4.0)
    assert -1.25 < slope < -0.75, (
        f"C50 {c1:.3g} uM at 1 h vs {c4:.3g} uM at 4 h: slope {slope:.2f}, expected ~ -1")


def test_a_short_paclitaxel_pulse_cannot_reach_half_kill_at_any_dose():
    """Time dependence in its starkest form: a tubulin binder kills only
    cells that attempt mitosis while it is present, so a 6 h pulse cannot
    halve an asynchronous colony however high the concentration (here
    ~70x the 72 h IC50), while a DNA-damage agent can do it in 1 h
    (previous gate). Georgiadis et al. 1997 (Clin Cancer Res 3:449)
    likewise found no half kill below 32 uM for 3 h paclitaxel."""
    surv = _surviving(_drug("paclitaxel"), 3.0, 6.0)
    assert surv > 0.5, f"6 h paclitaxel at 3 uM left only {surv:.3f}"


def test_spreading_a_fixed_dose_too_thin_switches_paclitaxel_off():
    """Threshold behaviour. The same AUC delivered over 72 h instead of
    24 h drops the concentration below the tubulin occupancy that
    triggers arrest, and the drug stops working — which is why AUC is
    the wrong exposure metric for this class."""
    pac = _drug("paclitaxel")
    auc = 1.8                                  # µM·h
    mid = _surviving(pac, auc / 24.0, 24.0)    # 0.075 µM for 24 h
    thin = _surviving(pac, auc / 72.0, 72.0)   # 0.025 µM for 72 h, at threshold
    assert mid < 0.7, f"paclitaxel ineffective even at the mid split ({mid:.3f})"
    assert thin > 0.9, (
        f"spreading the same AUC over 72 h should fall below threshold and do "
        f"almost nothing; got {thin:.3f}")
    assert thin > mid, "the thin split must be worse, not better"


def test_a_fixed_dose_has_an_interior_best_schedule_for_paclitaxel():
    """Both extremes fail for different reasons — too brief and few cells
    reach mitosis, too dilute and the drug never engages — so the best
    split of a fixed AUC lies in between. This is a prediction of the
    mechanism, not a fitted result."""
    pac = _drug("paclitaxel")
    auc = 1.8
    windows = (3.0, 24.0, 72.0)
    surv = [_surviving(pac, auc / w, w) for w in windows]
    assert surv[1] < surv[0] and surv[1] < surv[2], (
        f"expected an interior optimum across {windows}; got {surv}")
    assert surv[0] > 0.9 and surv[2] > 0.9, (
        f"both extremes should be nearly ineffective; got {surv[0]:.3f}, {surv[2]:.3f}")


def test_continuous_exposure_beats_an_early_washout_at_equal_concentration():
    """Sanity in the other direction: at the SAME concentration, leaving
    the drug on longer must not help less than washing it out early."""
    for name, conc in (("cisplatin", 20.0), ("doxorubicin", 0.2), ("paclitaxel", 0.1)):
        drug = _drug(name)
        early = _surviving(drug, conc, 6.0)
        full = _surviving(drug, conc, 72.0)
        assert full <= early + 1e-9, (
            f"{name}: 72 h exposure ({full:.3f}) should not leave more cells "
            f"alive than a 6 h pulse ({early:.3f})")


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
    print(f"{len(fns) - failed}/{len(fns)} schedule gates pass")
    raise SystemExit(1 if failed else 0)
