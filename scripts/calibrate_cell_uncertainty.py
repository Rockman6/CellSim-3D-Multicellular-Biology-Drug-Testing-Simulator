#!/usr/bin/env python3
"""How wrong is a predicted IC50, and does the stated error bar hold?

A prediction without an interval is not usable, and an interval that has
never been scored is not an interval. This calibrates one from the
engine's own held-out errors and then tests it the only way that means
anything: on drugs whose errors were not used to build it.

Method
------
Every held-out prediction in `gdsc_validation_results.csv` has a signed
log10 error against the GDSC reference. The spread of those errors is
the prediction uncertainty. A k-fold multiplicative interval

    [predicted / 10^(k*s),  predicted * 10^(k*s)]

where s is the standard deviation of the held-out log10 errors, should
cover the reference at a stated rate.

**Leave-one-drug-out.** For each drug the width is computed from the
OTHER drugs only, and scored on this one. A width fitted on the same
errors it is scored against would be circular, and would read as
well-calibrated however badly it generalised.

Two reference points are reported beside it:

* the **replicate span** of GDSC itself (the same line and drug screened
  more than once disagree by a median 5.3x), which is the floor: no
  interval can honestly be tighter than the data's own agreement;
* the **constant-IC50 null's** intervals, built the same way. If the
  engine's intervals are not narrower at equal coverage, the mechanism
  is adding nothing to the uncertainty picture either.

Usage:
    python scripts/calibrate_cell_uncertainty.py
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_CSV = REPO_ROOT / "benchmarks" / "cell" / "gdsc_validation_results.csv"
REFERENCE_CSV = REPO_ROOT / "benchmarks" / "cell" / "gdsc_reference.csv"
OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "cell_uncertainty.json"
# Nominal two-sided coverage and the Gaussian k that would deliver it.
LEVELS = {"68%": 1.0, "90%": 1.645, "95%": 1.96}


def load_rows() -> list[dict]:
    """Held-out predictions that have a real reference IC50."""
    out = []
    with RESULTS_CSV.open(newline="") as f:
        for r in csv.DictReader(f):
            if r["role"] != "held-out" or r["reference_kind"] != "value":
                continue
            if r["predicted_ic50_uM"] in ("inf", ""):
                continue
            ref = float(r["reference_ic50_uM"])
            out.append({"drug": r["drug"], "line": r["cell_line"], "ref": ref,
                        "pred": float(r["predicted_ic50_uM"]),
                        "null": float(r["null_ic50_uM"]),
                        "err": math.log10(float(r["predicted_ic50_uM"]) / ref),
                        "null_err": math.log10(float(r["null_ic50_uM"]) / ref)})
    return out


def replicate_spread() -> float:
    """Median max/min IC50 among repeat screens of the same line and drug."""
    groups: dict[tuple, list[float]] = {}
    with REFERENCE_CSV.open(newline="") as f:
        for r in csv.DictReader(f):
            if r["ic50_in_tested_range"] == "True":
                groups.setdefault((r["drug"], r["cell_line"]), []).append(float(r["ic50_uM"]))
    spreads = [max(v) / min(v) for v in groups.values() if len(v) > 1]
    return float(np.median(spreads)) if spreads else float("nan")


def coverage(rows: list[dict], widths: dict[str, float], key: str) -> dict:
    """Share of rows whose reference falls inside the interval."""
    out = {}
    for level, s in widths.items():
        k = LEVELS[level]
        hit = [abs(r[key]) <= k * s for r in rows]
        out[level] = {"covered": int(sum(hit)), "n": len(hit),
                      "coverage": float(np.mean(hit)) if hit else float("nan"),
                      "interval_fold": float(10 ** (k * s))}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    rows = load_rows()
    if len(rows) < 8:
        print(f"only {len(rows)} held-out predictions with a reference IC50; "
              f"run scripts/validate_gdsc.py first", file=sys.stderr)
        return 2
    drugs = sorted({r["drug"] for r in rows})
    spread = replicate_spread()
    print(f"{len(rows)} held-out predictions across {len(drugs)} drugs. "
          f"GDSC's own replicate spread: {spread:.1f}x (the floor on any interval).\n")

    # Leave-one-drug-out: the width for a drug comes from the other drugs.
    scored, scored_null = [], []
    per_drug = {}
    for d in drugs:
        inside = [r for r in rows if r["drug"] == d]
        outside = [r for r in rows if r["drug"] != d]
        if not inside or len(outside) < 4:
            continue
        s = float(np.std([r["err"] for r in outside], ddof=1))
        s_null = float(np.std([r["null_err"] for r in outside], ddof=1))
        per_drug[d] = {
            "n": len(inside), "width_log10_from_other_drugs": s,
            "own_rmse_log10": float(np.sqrt(np.mean([r["err"] ** 2 for r in inside]))),
            "coverage": coverage(inside, {k: s for k in LEVELS}, "err"),
            "null_coverage": coverage(inside, {k: s_null for k in LEVELS}, "null_err")}
        for r in inside:
            scored.append({**r, "s": s})
            scored_null.append({**r, "s": s_null})

    print(f"{'drug':<16}{'n':>3}{'  width (other drugs)':>22}{'  own RMSE':>11}"
          f"{'   68%':>8}{'   90%':>8}{'   95%':>8}")
    print("-" * 76)
    for d, v in per_drug.items():
        cov = v["coverage"]
        print(f"{d:<16}{v['n']:>3}{10 ** v['width_log10_from_other_drugs']:>20.2f}x"
              f"{10 ** v['own_rmse_log10']:>10.2f}x"
              f"{cov['68%']['coverage']:>8.2f}{cov['90%']['coverage']:>8.2f}"
              f"{cov['95%']['coverage']:>8.2f}")

    summary = {}
    for level, k in LEVELS.items():
        hit = [abs(r["err"]) <= k * r["s"] for r in scored]
        hit_null = [abs(r["null_err"]) <= k * r["s"] for r in scored_null]
        width = float(np.mean([10 ** (k * r["s"]) for r in scored]))
        width_null = float(np.mean([10 ** (k * r["s"]) for r in scored_null]))
        summary[level] = {"nominal": float(level.rstrip("%")) / 100,
                          "observed": float(np.mean(hit)), "n": len(hit),
                          "mean_interval_fold": width,
                          "null_observed": float(np.mean(hit_null)),
                          "null_mean_interval_fold": width_null}
        print(f"\n{level} interval, pooled over leave-one-drug-out folds:")
        print(f"   engine: covers {np.mean(hit):.0%} of {len(hit)} held-out lines, "
              f"mean interval +-{width:.1f}x")
        print(f"   constant-IC50 null: covers {np.mean(hit_null):.0%}, "
              f"mean interval +-{width_null:.1f}x")

    report = {"n_heldout": len(rows), "drugs": drugs,
              "gdsc_replicate_spread_fold": spread,
              "levels": LEVELS, "per_drug": per_drug, "pooled": summary}
    verdict = all(abs(v["observed"] - v["nominal"]) <= 0.15 for v in summary.values())
    report["calibrated_within_15_points"] = bool(verdict)
    print(f"\nCalibrated (every level within 15 points of nominal): {'YES' if verdict else 'NO'}")
    print("An interval this wide is not a precision claim: it is the honest width given\n"
          "that most line-to-line variation is not in the markers the engine can see\n"
          "(docs/VALIDATION.md). It is reported so a user knows what the number is worth.")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
