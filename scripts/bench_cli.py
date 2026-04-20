#!/usr/bin/env python3
"""`cellsim bench` — run all bundled validation benchmarks.

Runs the mini-pose-recovery bench + every
benchmarks/dock/*_calibration.yaml and prints a single aggregate
summary table. This is CellSim's "tell me how good you currently
are" command.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def _run_mini_bench():
    from tests.dock.test_mini_bench import _load_entries, _run_one
    print("=" * 70)
    print("MINI-BENCH: pose recovery on bundled cocrystals")
    print("=" * 70)
    rows = []
    for entry in _load_entries():
        print(f"\n[bench] {entry['name']}")
        rec = _run_one(entry, exhaustiveness=64, num_modes=9, seed=1)
        rows.append(rec)
        if rec["ok"]:
            tag = "✓ PASS" if rec.get("pass_gate") else "✗ FAIL"
            print(f"  {tag}  top1={rec.get('top1_rmsd_A')} Å  "
                  f"top3b={rec.get('top3_best_rmsd_A')} Å")
        else:
            print(f"  ERROR: {rec['reason'][:80]}")
    passed = sum(1 for r in rows if r.get("pass_gate"))
    n = len(rows)
    print(f"\n[bench] mini-bench: {passed}/{n} = {passed/n*100:.0f} %")
    return dict(passed=passed, n=n, rows=rows)


def _run_calibrations(exh: int = 32, cache_path: str | None = None):
    from src.cache import Cache
    from src.uq import run_calibration

    print()
    print("=" * 70)
    print("CALIBRATIONS: CellSim ΔG vs published K_d")
    print("=" * 70)

    cal_yamls = sorted((REPO / "benchmarks" / "dock").glob("*_calibration.yaml"))
    summaries = []
    cache = Cache(cache_path) if cache_path else None
    for y in cal_yamls:
        print(f"\n[bench] {y.name}")
        r = run_calibration(
            y, exhaustiveness=exh, num_modes=9,
            seed=1, cpu=2, cache=cache)
        print(f"  {r.summary()}")
        summaries.append(dict(
            yaml=y.name,
            n_ok=r.n_ok, n_points=r.n_points,
            pearson_r=r.pearson_r,
            spearman_rho=r.spearman_rho,
            mae_kcalmol=r.mae_kcalmol,
            rmse_kcalmol=r.rmse_kcalmol,
            q95_kcalmol=r.conformal_q95_kcalmol,
        ))
    if cache:
        cache.close()
    return summaries


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exh", type=int, default=32,
                    help="Vina exhaustiveness for calibrations "
                         "(default 32; use 64 for paper-quality)")
    ap.add_argument("--cache", default=None,
                    help="SQLite cache (speeds up repeat runs)")
    ap.add_argument("--json", action="store_true",
                    help="dump full results to stdout as JSON")
    ap.add_argument("--skip-mini-bench", action="store_true")
    ap.add_argument("--skip-calibrations", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.time()
    out = {"cellsim_bench": {"started_at": time.time()}}

    if not args.skip_mini_bench:
        out["mini_bench"] = _run_mini_bench()

    if not args.skip_calibrations:
        out["calibrations"] = _run_calibrations(
            exh=args.exh, cache_path=args.cache)

    out["cellsim_bench"]["wall_seconds"] = time.time() - t0

    print()
    print("=" * 70)
    print("AGGREGATE")
    print("=" * 70)
    if "mini_bench" in out:
        m = out["mini_bench"]
        print(f"  mini-bench pose recovery: "
              f"{m['passed']}/{m['n']}  ({m['passed']/m['n']*100:.0f} %)")
    if "calibrations" in out:
        print(f"  calibrations:")
        for s in out["calibrations"]:
            if s["spearman_rho"] is not None:
                print(f"    {s['yaml']:<38s} n={s['n_ok']}/{s['n_points']}  "
                      f"Pearson r = {s['pearson_r']:+.2f}  "
                      f"Spearman ρ = {s['spearman_rho']:+.2f}  "
                      f"MAE = {s['mae_kcalmol']:.2f}")
    print(f"  wall = {out['cellsim_bench']['wall_seconds']:.1f} s")
    print()
    print("  See BENCHMARKS.md for full context + caveats.")

    if args.json:
        print("\n" + json.dumps(out, indent=2, default=str))

    return 0


if __name__ == "__main__":
    sys.exit(main())
