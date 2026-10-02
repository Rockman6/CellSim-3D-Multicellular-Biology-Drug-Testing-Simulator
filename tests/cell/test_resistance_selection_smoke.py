#!/usr/bin/env python3
"""Phase-2 gates: treatment selects, and the schedule decides how hard.

Within-line, so the per-line prediction ceiling in docs/VALIDATION.md
does not apply. Each representative cell in the engine is a lineage
carrying its own anti-apoptotic reserve and drug accumulation, and a
lineage that survives keeps those values as it doubles — so selection
should emerge without being modelled explicitly.

The trade-off under test is the dose-intensity one: concentrating the
same total exposure into brief high-dose pulses kills more, and leaves a
more resistant remnant behind. How strongly depends on whether the drug
is concentration-driven or exposure-time-driven, which
`test_schedule_smoke.py` establishes separately.

Scope limits, stated so a passing test is not over-read: heterogeneity
is drawn once, so this is selection from a pre-existing tail rather than
evolution of new resistance, and lineages do not compete for space.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

LINE = get_line("A549")
P = Params()
K = calibrate_cycle_scale(LINE, P)
N_CELLS = 96
HOURS = 7 * 24.0


def _pulsed(conc, on_h, off_h, total_h):
    segs, t = [], 0.0
    while t < total_h:
        segs.append((t, conc))
        segs.append((min(t + on_h, total_h), 0.0))
        t += on_h + off_h
    return lambda th: next(c for s, c in reversed(segs) if th >= s - 1e-12)


def _run(drug, dose_fn):
    res = simulate(LINE, drug, 0.0, t_end_h=HOURS, n_cells=N_CELLS, seed=1, p=P,
                   k_cyc=K, dose_fn=dose_fn, record_every_h=24.0)
    pop = res.final
    w = pop.weight * pop.alive
    if w.sum() <= 0:
        return {"alive": 0.0, "reserve": float("nan"), "uptake": float("nan"),
                "resistance": float("nan")}
    reserve = float(np.average(pop.het_bcl2, weights=w))
    uptake = float(np.average(pop.het_uptake, weights=w))
    # A lineage can be resistant two ways: hold a bigger anti-apoptotic
    # reserve, or accumulate less drug. The combined score is the ratio,
    # so either route raises it and a comparison is not blind to one.
    return {"alive": float(w.sum()), "reserve": reserve, "uptake": uptake,
            "resistance": reserve / uptake if uptake > 0 else float("inf")}


def test_treatment_enriches_the_resistant_tail():
    """Survivors must be a biased sample: higher anti-apoptotic reserve
    and lower drug accumulation than the population that started, whose
    means are 1.0 by construction."""
    out = _run(get_drug("cisplatin"), lambda th: 15.0)
    assert out["alive"] > 0, "nothing survived; cannot measure selection"
    assert out["reserve"] > 1.02, (
        f"survivors not enriched for anti-apoptotic reserve ({out['reserve']:.3f})")
    assert out["uptake"] < 0.98, (
        f"survivors not enriched for low drug accumulation ({out['uptake']:.3f})")


def test_untreated_population_is_not_selected():
    """Specificity: with no drug the means must stay at 1.0."""
    out = _run(get_drug("cisplatin"), lambda th: 0.0)
    assert abs(out["reserve"] - 1.0) < 0.05, out["reserve"]
    assert abs(out["uptake"] - 1.0) < 0.05, out["uptake"]


def test_pulsing_kills_more_but_selects_harder():
    """Dose intensity trade-off, at MATCHED total exposure: 3x the
    concentration for a third of the time kills more of a
    concentration-driven drug, and leaves a more resistant remnant."""
    drug = get_drug("cisplatin")
    cont = _run(drug, lambda th: 15.0)
    puls = _run(drug, _pulsed(45.0, 24.0, 48.0, HOURS))
    assert puls["alive"] < cont["alive"], (
        f"pulsed should kill more at matched exposure: {puls['alive']:.3g} "
        f"vs continuous {cont['alive']:.3g}")
    if puls["alive"] > 0:
        assert puls["resistance"] > cont["resistance"], (
            f"the harder-killing arm should leave the more resistant remnant: "
            f"score {puls['resistance']:.3f} vs {cont['resistance']:.3f}")


def test_dose_intensity_decides_which_kind_of_resistance_is_selected():
    """A mechanistic prediction worth pinning. Against a brief 3x dose the
    anti-apoptotic reserve cannot save a cell, so survival depends almost
    entirely on accumulating less drug; a low continuous dose lets the
    reserve matter too. So the high-intensity arm should select far more
    sharply on uptake than on reserve."""
    drug = get_drug("cisplatin")
    cont = _run(drug, lambda th: 15.0)
    puls = _run(drug, _pulsed(45.0, 24.0, 48.0, HOURS))
    assert puls["alive"] > 0, "pulsed arm left nothing to measure"
    assert puls["uptake"] < cont["uptake"], (
        f"high-intensity dosing should select harder on drug accumulation: "
        f"{puls['uptake']:.3f} vs {cont['uptake']:.3f}")
    # and it does so without needing the reserve channel at all
    assert puls["uptake"] < 0.5, (
        f"expected strong selection on accumulation, got {puls['uptake']:.3f}")


def test_an_exposure_limited_drug_gains_least_from_pulsing():
    """Consistency with test_schedule_smoke.py: paclitaxel needs cells to
    transit mitosis while the drug is present, so compressing the same
    exposure into a third of the time wastes part of it. Its advantage
    from pulsing must be smaller than a concentration-driven drug's."""
    def gain(drug, conc):
        cont = _run(drug, lambda th, c=conc: c)
        puls = _run(drug, _pulsed(conc * 3.0, 24.0, 48.0, HOURS))
        return cont["alive"] / max(puls["alive"], 1e-9)
    pac_gain = gain(get_drug("paclitaxel"), 0.04)
    cis_gain = gain(get_drug("cisplatin"), 15.0)
    assert pac_gain < cis_gain, (
        f"paclitaxel should gain less from pulsing than cisplatin: "
        f"{pac_gain:.1f}x vs {cis_gain:.1f}x")


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
