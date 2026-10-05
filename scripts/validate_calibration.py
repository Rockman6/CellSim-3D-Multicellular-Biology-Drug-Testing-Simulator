#!/usr/bin/env python3
"""Does a calibrated engine predict what it was NOT calibrated on?

This is the question the whole of Phase 4 rests on. Fitting a constant so
the simulated IC50 matches a measured one is easy and proves nothing: the
constant was chosen to make that number right. The useful claim is the
next one — that the calibrated model then says something true about a
condition the calibration never saw — and that claim has to be measured.

The design is a hold-out, with ground truth we control:

1. **Truth.** Take a drug and perturb its one fitted constant to a value
   `v_true` that the calibrator is never told. This stands in for a lab's
   cells differing from the reference line.
2. **Measurement.** Simulate the 72 h continuous dose-response those
   cells would give, add plate-reader noise, and fit it with
   `cellsim.plate.fit_4pl` exactly as a real plate would be fitted. This
   yields a measured IC50 and its confidence interval — and nothing else.
3. **Calibration.** Hand only that IC50 and interval to
   `calibrate_potency`, which returns a constant and its interval.
4. **Hold-out prediction.** Ask the calibrated model about a condition
   that played no part in any of the above: survival under a SHORT PULSE,
   washed out long before the readout. Run it at both ends of the
   calibrated interval to get a band.
5. **Score.** Compute the same quantity with `v_true` and ask whether the
   band contains it.

Repeating that gives the coverage of a calibrated prediction, which is
the number a user actually needs: not "how well does it fit" but "how
often is the truth inside the band it gives me".

Two ways this can fail, and they mean different things. If the band is
too narrow the uncertainty is understated and the model is overconfident.
If the band is enormous the calibration has not constrained anything and
the prediction is not worth making — so the width is reported beside the
coverage, and both are needed to read the result.

Usage:
    python scripts/validate_calibration.py            # ~20 min
    python scripts/validate_calibration.py --quick    # ~6 min
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.calibrate import calibrate_potency  # noqa: E402
from cellsim.cell.engine import (Params, calibrate_cycle_scale, dose_response,  # noqa: E402
                                 ic50, simulate)
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.plate import fit_4pl  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "calibration_holdout.json"
LINE_NAME, DRUG_NAME = "A549", "cisplatin"
# How far the "lab's" cells differ from the reference line. Derived
# resistant lines reach two- to eight-fold (McDermott 2014), so a factor
# of two either way is a realistic difference to ask the engine to learn.
PERTURBATIONS = (0.5, 2.0)
PLATE_NOISE_SD = 0.05           # residual sd of a decent plate (validate_plate_fit.py)
PULSE_HOURS = 12.0              # the held-out condition: a pulse, washed out
READOUT_HOURS = 72.0


def simulated_plate(line, drug, centre_uM, *, n_cells, rng, k_cyc, p):
    """The dose-response those cells would give, with plate noise on top."""
    conc = centre_uM * np.logspace(-1.5, 1.5, 10)
    viab, _ = dose_response(line, drug, conc, t_end_h=READOUT_HOURS,
                            n_cells_per_conc=n_cells, k_cyc=k_cyc, p=p)
    noisy = np.clip(viab[1:] + rng.normal(0.0, PLATE_NOISE_SD, len(conc)), 0.0, 1.5)
    return conc, noisy


def pulse_survival(line, drug, conc_uM, *, n_cells, k_cyc, p, seed=1):  # noqa: D401
    """The held-out quantity: surviving fraction after a short pulse."""
    def dose(t_h):
        return conc_uM if t_h < PULSE_HOURS - 1e-9 else 0.0
    treated = simulate(line, drug, 0.0, t_end_h=READOUT_HOURS, n_cells=n_cells,
                       k_cyc=k_cyc, p=p, dose_fn=dose, record_every_h=READOUT_HOURS,
                       seed=seed)
    control = simulate(line, None, 0.0, t_end_h=READOUT_HOURS, n_cells=n_cells,
                       k_cyc=k_cyc, p=p, record_every_h=READOUT_HOURS, seed=seed)
    return float(treated.alive_weight[-1, 0] / max(control.alive_weight[-1, 0], 1e-12))


def one_trial(perturb: float, trial: int, *, n_cells: int, cal_cells: int,
              iters: int) -> dict | None:
    line, drug = get_line(LINE_NAME), get_drug(DRUG_NAME)
    p = Params()
    k = calibrate_cycle_scale(line, p)
    field = drug.fit_target
    v_true = getattr(drug, field) * perturb
    truth_drug = dataclasses.replace(drug, **{field: v_true})
    rng = np.random.default_rng(1000 * trial + int(perturb * 10))

    # 1-2. what the lab would measure, and nothing more
    centre = ic50(line, truth_drug, guess_uM=10.0, n_cells_per_conc=n_cells, k_cyc=k, p=p)
    if not np.isfinite(centre):
        return None
    conc, viab = simulated_plate(line, truth_drug, centre, n_cells=n_cells, rng=rng,
                                 k_cyc=k, p=p)
    try:
        fit = fit_4pl(conc, viab, n_boot=300, seed=trial)
    except ValueError:
        return None
    if not np.isfinite(fit.ic50_ci[0]):
        return None

    # 3. calibration sees only the fitted IC50 and its interval
    try:
        cal = calibrate_potency(LINE_NAME, DRUG_NAME, fit.ic50_uM, fit.ic50_ci,
                                n_cells=cal_cells, iters=iters, report_ic50=False)
    except ValueError:
        return None

    # 4-5. a condition none of the above involved. The band spans the
    # calibrated interval AND seeds, because the engine's own run-to-run
    # variability is as large as the calibration's (see predict_band).
    pulse_conc = 3.0 * fit.ic50_uM
    seeds = (1, 2, 3)
    from cellsim.calibrate import predict_band
    band_out = predict_band(
        cal, lambda d, sd: pulse_survival(line, d, pulse_conc, n_cells=n_cells,
                                          k_cyc=k, p=p, seed=sd), seeds=seeds)
    lo, hi = band_out["range"]
    band = [band_out["low"], band_out["value"], band_out["high"]]
    # The truth is itself a stochastic quantity; its median over the same
    # seeds is what the band is asked to contain.
    truth_runs = [pulse_survival(line, truth_drug, pulse_conc, n_cells=n_cells,
                                 k_cyc=k, p=p, seed=sd) for sd in seeds]
    truth = float(np.median(truth_runs))
    truth_spread = float(max(truth_runs) - min(truth_runs))
    return {"perturbation": perturb, "trial": trial,
            "true_constant": float(v_true), "calibrated_constant": cal.value,
            "constant_ratio": float(cal.value / v_true),
            "constant_ci": list(cal.ci),
            "measured_ic50_uM": fit.ic50_uM, "measured_ic50_ci": list(fit.ic50_ci),
            "true_ic50_uM": float(centre),
            "holdout_truth": truth, "holdout_truth_spread": truth_spread,
            "holdout_band": [lo, hi], "holdout_point": band[1],
            "covered": bool(lo <= truth <= hi),
            "band_width": float(hi - lo),
            "point_error": float(band[1] - truth)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--trials", type=int, default=None)
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    trials = a.trials if a.trials else (4 if a.quick else 10)
    n_cells, cal_cells, iters = (24, 20, 7) if a.quick else (32, 24, 8)

    print(f"{LINE_NAME} + {DRUG_NAME}. The 'lab' line differs from the reference by "
          f"{PERTURBATIONS[0]}x and {PERTURBATIONS[1]}x in {get_drug(DRUG_NAME).fit_target}.\n"
          f"Calibration sees ONLY a noisy 72 h dose-response (sd {PLATE_NOISE_SD}); the "
          f"held-out question is survival after a {PULSE_HOURS:g} h pulse read at "
          f"{READOUT_HOURS:g} h.\n")
    rows, t0 = [], time.time()
    for perturb in PERTURBATIONS:
        for trial in range(trials):
            r = one_trial(perturb, trial, n_cells=n_cells, cal_cells=cal_cells,
                          iters=iters)
            if r is None:
                continue
            rows.append(r)
            print(f"  {perturb:>4.1f}x trial {trial}: constant recovered "
                  f"{r['constant_ratio']:.2f}x; held-out truth {r['holdout_truth']:.3f}, "
                  f"band {r['holdout_band'][0]:.3f}-{r['holdout_band'][1]:.3f} "
                  f"-> {'covered' if r['covered'] else 'MISSED'}", flush=True)

    if not rows:
        print("no trial completed", file=sys.stderr)
        return 1
    cov = float(np.mean([r["covered"] for r in rows]))
    width = float(np.median([r["band_width"] for r in rows]))
    cr = np.array([r["constant_ratio"] for r in rows])
    err = np.array([abs(r["point_error"]) for r in rows])
    print(f"\n{len(rows)} trials in {time.time() - t0:.0f} s")
    print(f"  constant recovered:        median {np.median(cr):.2f}x "
          f"(range {cr.min():.2f}-{cr.max():.2f})")
    print(f"  held-out band covered:     {cov:.0%} of trials")
    print(f"  median band width:         {width:.3f} surviving fraction")
    print(f"  median point error:        {np.median(err):.3f}")
    ts = np.array([r["holdout_truth_spread"] for r in rows])
    print(f"  the truth's OWN seed spread: median {np.median(ts):.3f} — the engine is "
          f"stochastic,\n{'':29}so no band can be narrower than this and still be honest")
    verdict = cov >= 0.7 and width < 0.5
    print(f"\nA calibrated prediction of an unseen condition is "
          f"{'USABLE' if verdict else 'NOT yet usable'}: the band contains the truth "
          f"{cov:.0%} of the time and is {width:.2f} wide.")
    print("Both numbers are needed. A band that always covers because it spans every\n"
          "possible answer has told the user nothing; one that is narrow and wrong is\n"
          "worse than nothing.")
    report = {"line": LINE_NAME, "drug": DRUG_NAME, "perturbations": list(PERTURBATIONS),
              "plate_noise_sd": PLATE_NOISE_SD, "pulse_hours": PULSE_HOURS,
              "readout_hours": READOUT_HOURS, "n_trials": len(rows),
              "coverage": cov, "median_band_width": width,
              "median_constant_ratio": float(np.median(cr)),
              "median_abs_point_error": float(np.median(err)),
              "usable": bool(verdict), "trials": rows}
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
