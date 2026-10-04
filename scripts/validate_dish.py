#!/usr/bin/env python3
"""Validate the spatial dish against measured cells, in two settings.

1. **A sparse monolayer: the Cell Tracking Challenge HeLa movie.** Each
   of its two 46 h sequences is re-run as a closed field of the same size
   seeded with the same number of cells (scripts/ctc_reference.py derives
   the targets from the curated tracks). Compared: the count trajectory,
   and the cell-cycle durations a time-lapse would record through the
   same 46 h window (mean and CV of complete cycles, Kaplan-Meier share
   still undivided at 20 and 30 h), with the engine's cells on one clock
   (cycle_cv 0) and with measured cycle-time variability (0.25).

2. **A spheroid: DLD-1 under Grimes et al. 2014.** Their oxygen numbers
   are the dish's inputs, so the diffusion limit is a consistency check,
   not a test. What is tested is behaviour those inputs do not fix: when
   the core first turns necrotic relative to the 233 um limit, and how the
   necrotic radius tracks the outer radius afterwards against Grimes's
   own anoxic-core relation. The radial growth speed after onset is the
   CALIBRATION target for the dish's mechanical reach (DishParams.
   push_sites), reported with every value tried so the choice is visible.

Usage:
    python scripts/validate_dish.py            # ~40 min
    python scripts/validate_dish.py --quick    # ~8 min, fewer seeds and days
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import DishParams, run_dish  # noqa: E402
from cellsim.cell.engine import Params  # noqa: E402
from cellsim.cell.library import get_line  # noqa: E402

CTC_JSON = REPO_ROOT / "benchmarks" / "cell" / "ctc_hela_reference.json"
OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "dish_validation.json"
HELA_SPACING_UM = 20.0          # HeLa footprint in sparse culture; sets confluence only


def km_undivided(births_and_divisions, censored_ages, ages=(20.0, 30.0)) -> dict:
    """Kaplan-Meier share of cells still undivided at each age."""
    durations = [d - b for b, d in births_and_divisions]
    times = np.array(durations + list(censored_ages))
    events = np.array([1] * len(durations) + [0] * len(censored_ages))
    out = {}
    for a in ages:
        S = 1.0
        for t in np.unique(times[(events == 1) & (times <= a)]):
            S *= 1.0 - np.sum((times == t) & (events == 1)) / np.sum(times >= t)
        out[f"{a:g}"] = float(S)
    return out


def ctc_comparison(seeds, cycle_cvs) -> dict:
    ref = json.loads(CTC_JSON.read_text())
    hela = get_line("HeLa")
    out = {}
    for seq, r in ref["sequences"].items():
        nx = int(round(r["field_um"][0] / HELA_SPACING_UM))
        ny = int(round(r["field_um"][1] / HELA_SPACING_UM))
        t_end = (r["frames"] - 1) * r["frame_interval_h"]
        counts_ref = np.array(r["counts"], float)
        out[seq] = {"reference": {"counts_first_last": [r["counts"][0], r["counts"][-1]],
                                  "complete_cycle_mean_h": r["complete_cycle_mean_h"],
                                  "complete_cycle_cv": r["complete_cycle_cv"],
                                  "km_undivided": {k: r["km_undivided_at_age_h"][k] for k in ("20", "30")}},
                    "field_sites": [nx, ny], "runs": {}}
        for cv in cycle_cvs:
            rows = []
            for seed in seeds:
                res = run_dish(hela, t_end_h=t_end, n_seed=r["counts"][0], grid_sites=(nx, ny),
                               record_every_h=r["frame_interval_h"], seed=seed,
                               p=Params(cycle_cv=cv),
                               dp=DishParams(geometry="monolayer", spacing_um=HELA_SPACING_UM,
                                             field_dt_h=0.5, dt_h=0.02))
                n = np.array(res.record.n_live, float)
                m = min(len(n), len(counts_ref))
                rel = np.abs(n[:m] - counts_ref[:m]) / counts_ref[:m]
                d = res.dish
                cycles = d.cycles
                durations = np.array([b - a for a, b in cycles]) if cycles else np.array([np.nan])
                censored = [d.t_h - b for b in d.birth_t if np.isfinite(b)]
                rows.append({"final_count": float(n[-1]), "mean_abs_rel_error": float(rel.mean()),
                             "complete_cycle_mean_h": float(np.nanmean(durations)),
                             "complete_cycle_cv": float(np.nanstd(durations, ddof=1) / np.nanmean(durations))
                             if len(durations) > 2 else float("nan"),
                             "km_undivided": km_undivided(cycles, censored)})
            agg = {k: float(np.mean([row[k] for row in rows]))
                   for k in ("final_count", "mean_abs_rel_error", "complete_cycle_mean_h",
                             "complete_cycle_cv")}
            agg["km_undivided"] = {k: float(np.mean([row["km_undivided"][k] for row in rows]))
                                   for k in ("20", "30")}
            out[seq]["runs"][f"cycle_cv={cv:g}"] = {"seeds": list(seeds), "mean": agg, "per_seed": rows}
    return out


def spheroid_run(push: float, *, n_seed: int, days: float, grid: int, seed: int = 2) -> dict:
    dp = DishParams(geometry="spheroid", push_sites=push)
    t0 = time.time()
    res = run_dish(get_line("DLD-1"), t_end_h=days * 24.0, n_seed=n_seed, grid_sites=grid,
                   record_every_h=24.0, seed=seed, p=Params(cycle_cv=0.25), dp=dp)
    rec = res.record
    R = np.array(rec.radius_um)
    Rn = np.array(rec.necrotic_radius_um)
    t = np.array(rec.t_h) / 24.0
    onset = next((i for i, v in enumerate(rec.n_necrotic) if v > 0), None)
    wall = (grid / 2.0) * dp.spacing_um
    usable = R < 0.85 * wall                    # stop before the lattice walls interfere
    after = [i for i in range(len(t)) if onset is not None and i >= onset and usable[i]]
    speed = (float(np.polyfit(t[after], R[after], 1)[0]) if len(after) >= 3 else float("nan"))
    # Grimes's anoxic-core relation: r_l^2 = R^2 - r_n^2 - 2 r_n^3 (1/r_n - 1/R)
    rl = [math.sqrt(max(R[i] ** 2 - Rn[i] ** 2 - 2 * Rn[i] ** 3 * (1 / Rn[i] - 1 / R[i]), 0.0))
          for i in after if Rn[i] > 0]
    return {"push_sites": push, "n_seed": n_seed, "days": days, "wall_s": round(time.time() - t0, 1),
            "day": t.tolist(), "radius_um": R.tolist(), "necrotic_radius_um": Rn.tolist(),
            "quiescent_frac": rec.quiescent_frac, "hypoxic_frac": rec.hypoxic_frac,
            "onset_day": float(t[onset]) if onset is not None else None,
            "radius_at_onset_um": float(R[onset]) if onset is not None else None,
            "radial_speed_after_onset_um_per_day": speed,
            "implied_diffusion_limit_um": float(np.median(rl)) if rl else None}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    seeds = (1, 2) if a.quick else (1, 2, 3, 4, 5)
    report = {"ctc_hela": ctc_comparison(seeds, (0.0, 0.25))}
    for seq, r in report["ctc_hela"].items():
        ref = r["reference"]
        print(f"CTC HeLa sequence {seq} ({r['field_sites'][0]}x{r['field_sites'][1]} sites): "
              f"cells {ref['counts_first_last'][0]} -> {ref['counts_first_last'][1]}, complete cycles "
              f"{ref['complete_cycle_mean_h']:.1f} h (CV {ref['complete_cycle_cv']:.2f}), "
              f"undivided at 20/30 h {ref['km_undivided']['20']:.2f}/{ref['km_undivided']['30']:.2f}")
        for label, run in r["runs"].items():
            m = run["mean"]
            print(f"   engine {label:<13} -> {m['final_count']:.0f} cells (mean |rel err| "
                  f"{m['mean_abs_rel_error']:.2f}); complete cycles {m['complete_cycle_mean_h']:.1f} h "
                  f"(CV {m['complete_cycle_cv']:.2f}); undivided at 20/30 h "
                  f"{m['km_undivided']['20']:.2f}/{m['km_undivided']['30']:.2f}")

    days = 6 if a.quick else 12
    pushes = (1.0, 2.0) if a.quick else (1.0, 2.0, 3.0, float("inf"))
    report["spheroid_calibration"] = []
    print("\nDLD-1 spheroid from 3000 cells; Grimes 2014: ~15 um/day after the core turns anoxic")
    for push in pushes:
        r = spheroid_run(push, n_seed=3000, days=days if math.isfinite(push) else min(days, 6),
                         grid=64)
        report["spheroid_calibration"].append(r)
        print(f"   reach {push:>4g} sites: necrosis from day {r['onset_day']} (R {r['radius_at_onset_um']:.0f} um), "
              f"then {r['radial_speed_after_onset_um_per_day']:.1f} um/day; implied diffusion limit "
              f"{r['implied_diffusion_limit_um']} um ({r['wall_s']:.0f} s)")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
