#!/usr/bin/env python
"""Layer 1.6 — are CellSim's error bars honest? (exit criterion 4)

Criterion 4 requires calibration error <= 10 % on a HELD-OUT set. The
existing conformal test (tests/uq/test_conformal_smoke.py) proves the
split-conformal MATH is implemented correctly — but it does so on
SYNTHETIC numbers. It has never been checked whether CellSim's *real*
docking predictions carry truthful intervals. That is what this does.

Method (honest split-conformal):
  1. Dock every compound in the real calibration bundles (streptavidin /
     trypsin / EGFR), giving (predicted dG, experimental dG) pairs.
  2. Split the pooled pairs into a calibration half and a held-out test
     half (seeded, deterministic).
  3. Fit the conformal quantile q on the calibration half only.
  4. Measure empirical coverage on the held-out half: what fraction of
     true values actually land inside the claimed interval?
  5. calibration error = |empirical coverage - nominal (1 - alpha)|.

"Honest" means: if we claim a 95 % interval, the truth should land
inside ~95 % of the time. Much lower = the bars are too narrow and the
tool is overconfident (the dangerous direction). Much higher = the bars
are too wide and the tool is useless-but-safe.

Small-n caveat is reported, not hidden: with ~10 held-out points the
binomial error on a coverage estimate is large, so a Wilson 95 %
interval on the coverage itself is printed alongside.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

DEFAULT_BUNDLES = [
    "benchmarks/dock/streptavidin_calibration.yaml",
    "benchmarks/dock/trypsin_calibration.yaml",
    "benchmarks/dock/egfr_calibration.yaml",
]


def _wilson(k: int, n: int, z: float = 1.96) -> tuple:
    """Wilson score interval for a binomial proportion (small-n honest)."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - half) / d, (c + half) / d)


def collect_pairs(bundles: list, exhaustiveness: int, cpu: int) -> list:
    """Dock each bundle; return [(name, pred, expt, bundle)] for ok points."""
    from src.uq import run_calibration

    pairs = []
    for b in bundles:
        p = REPO_ROOT / b
        if not p.exists():
            print(f"  ! missing bundle {b}; skipped", flush=True)
            continue
        print(f"  docking {Path(b).stem} ...", flush=True)
        try:
            r = run_calibration(p, exhaustiveness=exhaustiveness, cpu=cpu)
        except TypeError:
            r = run_calibration(p)
        if not r.ok:
            print(f"  ! {b}: {r.reason}", flush=True)
            continue
        n_ok = 0
        for pt in r.points:
            if pt.dG_pred_kcalmol is None or pt.dG_expt_kcalmol is None:
                continue
            pairs.append((pt.name, float(pt.dG_pred_kcalmol),
                          float(pt.dG_expt_kcalmol), Path(b).stem))
            n_ok += 1
        print(f"    -> {n_ok} usable pairs", flush=True)
    return pairs


def coverage_report(pairs: list, alpha: float, seed: int) -> dict:
    """Split-conformal fit on half, measure coverage on the held-out half."""
    from src.uq import ConformalBounds

    rng = random.Random(seed)
    shuffled = pairs[:]
    rng.shuffle(shuffled)
    n = len(shuffled)
    n_cal = n // 2
    cal, test = shuffled[:n_cal], shuffled[n_cal:]
    if len(cal) < 3 or len(test) < 3:
        return {"error": f"too few pairs (n={n}) to split-conformal"}

    cb = ConformalBounds(alpha=alpha)
    info = cb.calibrate([p for _, p, _, _ in cal],
                        [t for _, _, t, _ in cal])
    cov = cb.check_coverage([p for _, p, _, _ in test],
                            [t for _, _, t, _ in test])
    nominal = 1.0 - alpha
    k = int(round(cov * len(test)))
    lo, hi = _wilson(k, len(test))
    return {
        "n_total": n,
        "n_cal": len(cal),
        "n_test": len(test),
        "alpha": alpha,
        "nominal_coverage": nominal,
        "empirical_coverage": cov,
        "coverage_wilson95": [lo, hi],
        "calibration_error": abs(cov - nominal),
        "calibration_error_gate": 0.10,
        "q_kcalmol": info.get("q_kcalmol"),
        "interval_width_kcalmol": 2 * (info.get("q_kcalmol") or 0.0),
        "test_names": [nm for nm, _, _, _ in test],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundles", nargs="*", default=DEFAULT_BUNDLES)
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="miscoverage rate (0.05 -> 95%% interval)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--out-json", default=None)
    args = ap.parse_args(argv)

    print("[uq-coverage] docking real calibration bundles")
    pairs = collect_pairs(args.bundles, args.exhaustiveness, args.cpu)
    print(f"[uq-coverage] pooled usable pairs: {len(pairs)}")
    if len(pairs) < 6:
        print("  ! not enough real pairs to measure coverage")
        return 1

    rep = coverage_report(pairs, args.alpha, args.seed)
    if "error" in rep:
        print(f"  ! {rep['error']}")
        return 1

    print()
    print(f"  claimed interval : {rep['nominal_coverage']:.0%} "
          f"(±{rep['q_kcalmol']:.2f} kcal/mol, width "
          f"{rep['interval_width_kcalmol']:.2f})")
    print(f"  calibrated on    : n={rep['n_cal']}")
    print(f"  held-out test    : n={rep['n_test']}")
    print(f"  ACTUAL coverage  : {rep['empirical_coverage']:.0%}  "
          f"(Wilson 95% CI {rep['coverage_wilson95'][0]:.0%}"
          f"-{rep['coverage_wilson95'][1]:.0%})")
    err = rep["calibration_error"]
    print(f"  calibration error: {err:.0%}  [gate <= 10%] "
          f"{'PASS' if err <= 0.10 else 'FAIL'}")
    if rep["empirical_coverage"] < rep["nominal_coverage"]:
        print("  NB under-coverage = error bars too NARROW "
              "(overconfident — the dangerous direction).")
    print(f"  n_test={rep['n_test']} is small; the Wilson interval above "
          "is the honest uncertainty on this number.")

    if args.out_json:
        Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_json).write_text(json.dumps(rep, indent=2))
        print(f"  wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
