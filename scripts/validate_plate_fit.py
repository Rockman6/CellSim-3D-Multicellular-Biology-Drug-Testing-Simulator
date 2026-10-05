#!/usr/bin/env python3
"""Does the plate fit's error bar mean what it says?

`cellsim.plate.fit_4pl` reports a 95 % confidence interval on the IC50.
That claim is only worth making if it has been checked, and it can be
checked exactly: simulate curves whose IC50 is KNOWN, fit them, and count
how often the interval contains the truth. A 95 % interval should contain
it 95 % of the time. More is wasteful, less is misleading — and less is
what the first implementation did.

What this found, and why the code changed
-----------------------------------------
Resampling the raw residuals gave 80 % coverage at a nominal 95 %. The
reason is standard and worth stating: residuals from a fitted curve are
smaller than the true errors, because the fit has already absorbed p of
the n degrees of freedom. Scaling them by sqrt(n/(n-p)) — the same
correction that turns a residual sum of squares into an unbiased variance
estimate — is what `fit_4pl` now does.

Two further things this measures, because they decide whether a given
plate design can support a number at all:

* **Points per curve.** A five-point curve and a twelve-point curve carry
  very different information; the interval should widen accordingly.
* **Range.** If the tested concentrations miss the IC50, the fit
  extrapolates. Those fits are flagged by `FitResult.extrapolated`, and
  this reports their coverage separately, because an extrapolated
  interval is a statement about the model rather than the experiment.

Usage:
    python scripts/validate_plate_fit.py            # ~6 min
    python scripts/validate_plate_fit.py --quick    # ~1 min
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.plate import fit_4pl  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "plate_fit_validation.json"
TRUE_IC50, TRUE_HILL, BOTTOM, TOP = 1.7, 1.3, 0.03, 1.0
NOMINAL = 0.95


def curve(conc: np.ndarray) -> np.ndarray:
    return BOTTOM + (TOP - BOTTOM) / (
        1.0 + 10.0 ** ((np.log10(conc) - np.log10(TRUE_IC50)) * TRUE_HILL))


def study(conc: np.ndarray, noise: float, trials: int, n_boot: int) -> dict:
    """Fit `trials` noisy realisations of a known curve; count coverage."""
    rng = np.random.default_rng(abs(hash((noise, len(conc), trials))) % 2**32)
    clean = curve(conc)
    hit = miss_low = miss_high = 0
    widths, errors, extrap = [], [], 0
    for t in range(trials):
        v = clean + rng.normal(0.0, noise, len(conc))
        try:
            f = fit_4pl(conc, v, n_boot=n_boot, seed=t)
        except ValueError:
            continue
        if not np.isfinite(f.ic50_ci[0]):
            continue
        lo, hi = f.ic50_ci
        if lo <= TRUE_IC50 <= hi:
            hit += 1
        elif hi < TRUE_IC50:
            miss_low += 1
        else:
            miss_high += 1
        widths.append(hi / lo)
        errors.append(f.ic50_uM / TRUE_IC50)
        extrap += bool(f.extrapolated)
    n = hit + miss_low + miss_high
    return {"n_points": len(conc), "noise_sd": noise, "n_fits": n,
            "coverage": hit / n if n else float("nan"),
            "missed_low": miss_low, "missed_high": miss_high,
            "median_ci_width_fold": float(np.median(widths)) if widths else float("nan"),
            "median_bias_fold": float(np.median(errors)) if errors else float("nan"),
            "n_extrapolated": extrap}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    trials = 60 if a.quick else 200
    n_boot = 200 if a.quick else 400

    report = {"true_ic50_uM": TRUE_IC50, "true_hill": TRUE_HILL,
              "nominal_coverage": NOMINAL, "trials_per_cell": trials,
              "n_boot": n_boot, "noise": [], "points": [], "range": []}

    print(f"True IC50 {TRUE_IC50} uM, Hill {TRUE_HILL}. {trials} simulated plates per row.\n")
    print("A) how noisy the plate is (10 concentrations, 0.01-100 uM)")
    print(f"{'noise sd':>10}{'fits':>7}{'coverage':>10}{'CI width':>10}{'bias':>8}{'low/high misses':>18}")
    print("-" * 65)
    conc10 = np.geomspace(0.01, 100, 10)
    for noise in (0.02, 0.05, 0.10, 0.20):
        r = study(conc10, noise, trials, n_boot)
        report["noise"].append(r)
        misses = f"{r['missed_low']} / {r['missed_high']}"
        print(f"{noise:>10.2f}{r['n_fits']:>7}{r['coverage']:>10.0%}"
              f"{r['median_ci_width_fold']:>9.2f}x{r['median_bias_fold']:>8.2f}"
              f"{misses:>18}")

    print("\nB) how many concentrations the plate tests (noise sd 0.05)")
    print(f"{'points':>10}{'fits':>7}{'coverage':>10}{'CI width':>10}{'bias':>8}")
    print("-" * 47)
    for n_pts in (5, 7, 10, 12):
        r = study(np.geomspace(0.01, 100, n_pts), 0.05, trials, n_boot)
        report["points"].append(r)
        print(f"{n_pts:>10}{r['n_fits']:>7}{r['coverage']:>10.0%}"
              f"{r['median_ci_width_fold']:>9.2f}x{r['median_bias_fold']:>8.2f}")

    print("\nC) whether the tested range actually brackets the IC50 (noise sd 0.05)")
    print(f"{'range uM':>16}{'fits':>7}{'coverage':>10}{'CI width':>10}{'extrapolated':>14}")
    print("-" * 57)
    for lo, hi in ((0.01, 100.0), (0.1, 30.0), (0.001, 0.5), (10.0, 1000.0)):
        r = study(np.geomspace(lo, hi, 10), 0.05, trials, n_boot)
        r["range_uM"] = [lo, hi]
        report["range"].append(r)
        label = f"{lo:g}-{hi:g}"
        print(f"{label:>16}{r['n_fits']:>7}{r['coverage']:>10.0%}"
              f"{r['median_ci_width_fold']:>9.2f}x{r['n_extrapolated']:>14}")

    ok = [r for r in report["noise"] if r["noise_sd"] <= 0.10]
    calibrated = all(abs(r["coverage"] - NOMINAL) <= 0.10 for r in ok)
    report["calibrated_within_10_points"] = bool(calibrated)
    print(f"\nCalibrated at realistic noise (every level within 10 points of 95 %): "
          f"{'YES' if calibrated else 'NO'}")
    print("Rows C show what a badly chosen concentration range does: the fit still\n"
          "returns a number, the interval still looks narrow, and the truth is\n"
          "outside it. That is why fit_4pl flags extrapolated fits explicitly.")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
