#!/usr/bin/env python3
"""Phase-2 gates: treatment selects, and the schedule decides on what.

Within-line, so the per-line prediction ceiling in docs/VALIDATION.md
does not apply. Each representative cell in the engine is a lineage
carrying its own anti-apoptotic reserve and drug accumulation, and a
lineage that survives keeps those values as it doubles — so selection
emerges without being modelled explicitly.

Selection is measured on a trait grid (cellsim.cell.selection), not by
averaging the few random lineages a harsh schedule leaves: an earlier
version of these gates did that, rested on one or two surviving
lineages, and asserted the opposite of what the exact measurement
shows. The trade-off under test, at matched total exposure:

* continuous low-dose exposure sets a damage steady state, and whether
  it stays below threshold depends on how much drug a cell takes up, so
  it selects mainly on accumulation;
* intermittent high-dose pulses push everyone over threshold briefly,
  and what decides survival is whether the apoptotic reserve can absorb
  the transient, so selection shifts toward the reserve.

Scope limits, stated so a passing test is not over-read: heterogeneity
is drawn once, so this is selection from a pre-existing tail rather than
evolution of new resistance, and lineages do not compete for space.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.selection import pulsed, summarize, survival_map  # noqa: E402

LINE = get_line("A549")
P = Params()
K = calibrate_cycle_scale(LINE, P)
HOURS = 7 * 24.0
RESERVE = np.geomspace(0.5, 3.0, 7)
UPTAKE = np.geomspace(0.15, 1.5, 7)

# These phenotypes belong to the mechanism, not to the fitted potency, so
# each drug is held at the potency its doses were designed for (the fit
# before the October 2026 refit). Cisplatin and paclitaxel enter the
# engine only as potency x concentration, so this is exactly a rescaling
# of the doses, and refitting the library cannot carry a dose across a
# threshold and break a test that is about shape.
DESIGN_POTENCY = {"cisplatin": {"k_damage_per_uM_h": 0.000951},
                  "paclitaxel": {"partition": 0.153}}
DOSE_UM = {"cisplatin": 15.0, "paclitaxel": 0.04}   # ~1x IC50 at design potency


def _drug(name: str):
    return dataclasses.replace(get_drug(name), **DESIGN_POTENCY[name])


_CACHE: dict = {}


def _arm(name: str, arm: str) -> dict:
    """Selection summary for one drug and arm (cached across gates)."""
    key = (name, arm)
    if key not in _CACHE:
        c = DOSE_UM[name]
        fn = {"continuous": lambda t_h: c, "pulsed": pulsed(3.0 * c),
              "untreated": lambda t_h: 0.0}[arm]
        S = survival_map(LINE, _drug(name), fn, hours=HOURS, reps=2, k_cyc=K, p=P,
                         reserve=RESERVE, uptake=UPTAKE)
        _CACHE[key] = summarize(S, reserve=RESERVE, uptake=UPTAKE, sigma=P.het_sigma)
    return _CACHE[key]


def test_untreated_population_is_not_selected():
    """Specificity: with no drug every grid point survives and the
    survivors are the starting population."""
    s = _arm("cisplatin", "untreated")
    assert abs(s["population_survival"] - 1.0) < 1e-9, s["population_survival"]
    assert abs(s["survivor_mean_reserve"] - s["baseline_mean_reserve"]) < 1e-9
    assert abs(s["survivor_mean_uptake"] - s["baseline_mean_uptake"]) < 1e-9


def test_treatment_enriches_the_resistant_tail():
    """Survivors of continuous exposure must accumulate less drug than
    the population that started."""
    s = _arm("cisplatin", "continuous")
    assert 0.0 < s["population_survival"] < 1.0, s["population_survival"]
    assert s["survivor_mean_uptake"] < 0.9 * s["baseline_mean_uptake"], (
        f"survivors not enriched for low drug accumulation ({s['survivor_mean_uptake']:.3f})")


def test_pulsing_kills_more_at_matched_exposure():
    """Dose intensity: the same total exposure as 3x the concentration
    for a third of the time kills more of a concentration-driven drug."""
    cont, puls = _arm("cisplatin", "continuous"), _arm("cisplatin", "pulsed")
    assert puls["population_survival"] < 0.5 * cont["population_survival"], (
        f"pulsed {puls['population_survival']:.3g} vs continuous "
        f"{cont['population_survival']:.3g}")


def test_pulsing_shifts_selection_toward_the_apoptotic_reserve():
    """The corrected direction. Continuous exposure selects mainly on
    accumulation; brief high-dose pulses select far harder on the
    reserve, because a reserve can absorb a transient but not a
    sustained insult."""
    cont, puls = _arm("cisplatin", "continuous"), _arm("cisplatin", "pulsed")
    assert puls["survivor_mean_reserve"] > cont["survivor_mean_reserve"] + 0.3, (
        f"pulsed survivors' reserve {puls['survivor_mean_reserve']:.3f} vs continuous "
        f"{cont['survivor_mean_reserve']:.3f}")
    assert puls["uptake_over_reserve"] < cont["uptake_over_reserve"], (
        f"boundary should tilt toward the reserve under pulsing: |b_u/b_r| "
        f"{puls['uptake_over_reserve']:.2f} vs continuous {cont['uptake_over_reserve']:.2f}")


def test_an_exposure_limited_drug_gains_least_from_pulsing():
    """Consistency with test_schedule_smoke.py: paclitaxel needs cells to
    transit mitosis while the drug is present, so compressing the same
    exposure into a third of the time wastes part of it. Its gain from
    pulsing must be far smaller than a concentration-driven drug's."""
    def gain(name):
        return (_arm(name, "continuous")["population_survival"]
                / max(_arm(name, "pulsed")["population_survival"], 1e-9))
    pac, cis = gain("paclitaxel"), gain("cisplatin")
    assert pac < 0.5 * cis, f"paclitaxel gains {pac:.2f}x from pulsing vs cisplatin {cis:.2f}x"


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
    print(f"{len(fns) - failed}/{len(fns)} resistance-selection gates pass")
    raise SystemExit(1 if failed else 0)
