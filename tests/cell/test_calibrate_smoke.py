#!/usr/bin/env python3
"""Phase-4 gates: tuning the engine to someone else's cells.

The strongest test available here is a round trip, because it needs no
external truth: take a drug's library constant, ask the engine what IC50
it produces, then hand that IC50 back to the calibrator and check it
recovers the constant it was never told. Anything less than that is
checking the code against itself.

The rest of the gates are about honesty rather than accuracy: an
interval that is propagated instead of invented, a refusal when the
measurement is outside what the mechanism can produce, and a flag when
the assay was too short to determine a doubling time.

Kept small (few cells, few bisection steps) so it runs in CI; a full
calibration uses more of both.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.calibrate import (calibrate_potency, doubling_time_from_controls,  # noqa: E402
                               predict_band)
from cellsim.cell.engine import Params, calibrate_cycle_scale, ic50  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

LINE = get_line("A549")
P = Params()
K = calibrate_cycle_scale(LINE, P)
N_CELLS, ITERS = 16, 8


def test_doubling_time_comes_back_out_of_the_control_wells():
    """Four-fold growth over 72 h is two doublings, so 36 h."""
    rng = np.random.default_rng(1)
    t0 = 1000 * rng.normal(1.0, 0.03, 12)
    final = 4000 * rng.normal(1.0, 0.03, 12)
    d = doubling_time_from_controls(72.0, t0, final)
    assert abs(d.hours - 36.0) < 2.0, d.summary()
    assert d.ci[0] < d.hours < d.ci[1]
    assert abs(d.doublings_observed - 2.0) < 0.1
    assert not d.note, d.note


def test_an_assay_too_short_to_determine_growth_says_so():
    """Under one doubling, the log ratio is small and its relative error
    large. The number is still returned, with the caveat attached."""
    rng = np.random.default_rng(2)
    t0 = 1000 * rng.normal(1.0, 0.05, 12)
    final = 1400 * rng.normal(1.0, 0.05, 12)
    d = doubling_time_from_controls(24.0, t0, final)
    assert d.doublings_observed < 1.0
    assert "doublings" in d.note, d.note
    assert (d.ci[1] - d.ci[0]) > 5.0, "a half-doubling assay should give a wide interval"


def test_a_control_that_did_not_grow_is_refused():
    try:
        doubling_time_from_controls(72.0, [1000.0] * 6, [980.0] * 6)
    except ValueError as e:
        assert "did not grow" in str(e), e
    else:
        raise AssertionError("a flat control must be refused")


def test_calibration_recovers_a_constant_it_was_not_told():
    """The round trip. Perturb cisplatin's potency, ask the engine for the
    IC50 that produces, then calibrate from that IC50 alone and check the
    constant comes back."""
    import dataclasses
    drug = get_drug("cisplatin")
    true_value = drug.k_damage_per_uM_h * 1.6
    perturbed = dataclasses.replace(drug, k_damage_per_uM_h=true_value)
    target = ic50(LINE, perturbed, guess_uM=10.0, n_cells_per_conc=N_CELLS, k_cyc=K, p=P)
    assert np.isfinite(target), "the perturbed drug should still have an IC50"
    cal = calibrate_potency("A549", "cisplatin", target, n_cells=N_CELLS, iters=ITERS)
    assert cal.fit_target == "k_damage_per_uM_h"
    ratio = cal.value / true_value
    assert 0.8 < ratio < 1.25, (
        f"calibration recovered {cal.value:.4g} for a true {true_value:.4g} "
        f"({ratio:.2f}x)")


def test_a_measured_interval_becomes_a_constant_interval():
    """Uncertainty must be propagated, not invented: a wider measurement
    gives a wider constant, and the point estimate sits inside."""
    narrow = calibrate_potency("A549", "cisplatin", 8.0, (7.0, 9.0),
                               n_cells=N_CELLS, iters=ITERS)
    wide = calibrate_potency("A549", "cisplatin", 8.0, (4.0, 16.0),
                             n_cells=N_CELLS, iters=ITERS)
    for cal in (narrow, wide):
        assert cal.ci[0] <= cal.value <= cal.ci[1], cal.summary()
        assert not cal.note, cal.note
    assert (wide.ci[1] / wide.ci[0]) > (narrow.ci[1] / narrow.ci[0]) * 1.5, (
        f"a four-fold measurement should give a wider constant than a 1.3-fold one: "
        f"{wide.ci[1] / wide.ci[0]:.2f} vs {narrow.ci[1] / narrow.ci[0]:.2f}")


def test_no_interval_in_means_no_interval_out():
    cal = calibrate_potency("A549", "cisplatin", 8.0, None, n_cells=N_CELLS, iters=ITERS)
    assert cal.ci == (cal.value, cal.value)
    assert "no interval" in cal.note, cal.note


def test_an_unreachable_measurement_is_refused_with_a_reason():
    """A measured IC50 the mechanism cannot produce at any constant means
    something else about the line is wrong, and saying so beats returning
    the nearest value."""
    try:
        calibrate_potency("A549", "paclitaxel", 1e-9, n_cells=N_CELLS, iters=ITERS)
    except ValueError as e:
        assert "closest the engine gets" in str(e), e
    else:
        raise AssertionError("an unreachable IC50 must be refused")


def test_predictions_come_back_as_a_band():
    """The point of calibrating with an interval: a later prediction is a
    range, and a difference smaller than that range is not a result."""
    cal = calibrate_potency("A549", "cisplatin", 8.0, (5.0, 13.0),
                            n_cells=N_CELLS, iters=ITERS)
    band = predict_band(cal, lambda d: ic50(LINE, d, guess_uM=8.0,
                                            n_cells_per_conc=N_CELLS, k_cyc=K, p=P))
    assert set(band) >= {"low", "value", "high", "spread", "range"}
    assert band["range"][0] < band["value"] < band["range"][1], band
    assert band["spread"] > 1.0, f"a 2.6-fold measurement should give a real spread: {band}"


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
        except Exception as e:                      # noqa: BLE001 - report, don't crash
            failed += 1
            print(f"  ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} calibration gates pass")
    raise SystemExit(1 if failed else 0)
