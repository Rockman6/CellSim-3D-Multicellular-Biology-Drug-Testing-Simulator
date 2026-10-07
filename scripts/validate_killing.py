#!/usr/bin/env python3
"""Does simulated immune killing behave the way live imaging says it does?

Weigelin et al. 2021 (Nat Commun 12:5217) imaged cytotoxic T cells
attacking tumour cells and found a single contact rarely kills. Each
delivers a SUBLETHAL hit the target repairs, and death needs about three
hits within roughly 50 minutes of one another, usually from several CTLs
transiting between targets. That replaced the older one-lethal-hit
picture, which is what this engine implemented first.

The checks below follow from the MECHANISM rather than from any rate, so
they cannot be tuned into passing:

1. COOPERATIVITY. If hits must land close together, then at low effector
   density most hits are wasted on cells that recover, and the number of
   hits spent per kill must fall steeply as density rises.
2. THE REPAIR WINDOW MATTERS. Shortening the window must make killing
   harder and lengthening it easier, with everything else fixed.
3. ONE-HIT KILLING IS THE OLD MODEL. Setting hits_to_kill = 1 must
   restore efficient killing at low density.
4. NO FREE KILLS. Kills must never exceed hits / hits_to_kill.

What this does NOT check is the absolute lysis at a given E:T against a
measured chromium-release curve. No such curve has been digitised here,
so the rates remain literature values and the magnitudes unvalidated.

Usage:
    python scripts/validate_killing.py [--quick]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import Dish, DishParams  # noqa: E402
from cellsim.cell.engine import Params  # noqa: E402
from cellsim.cell.library import get_line  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "killing_validation.json"


def assay(et, *, hits_to_kill=3, decay_h=50.0 / 60.0, hours=24.0, seeds=(1, 2),
          n_tumour=400, grid=40):
    live, kills, hits = [], [], []
    for sd in seeds:
        dp = DishParams(geometry="monolayer", spacing_um=20.0, field_dt_h=0.5,
                        dt_h=0.02, effector_hits_to_kill=hits_to_kill,
                        effector_hit_decay_h=decay_h)
        d = Dish(get_line("A549"), n_seed=n_tumour, grid_sites=(grid, grid),
                 seed=int(sd), p=Params(), dp=dp)
        d.add_effectors(int(round(et * n_tumour)))
        for _ in range(int(round(hours / dp.field_dt_h))):
            d.advance(None)
        live.append(int(d.cells.alive.sum()))
        kills.append(d.n_killed_by_effectors)
        hits.append(d.n_hits_delivered)
    return (float(np.mean(live)), float(np.mean(kills)), float(np.mean(hits)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    seeds = (1,) if a.quick else (1, 2, 3)
    ratios = (0.0, 0.1, 0.3, 1.0, 3.0)

    print("A549 monolayer, 24 h, three-hit additive killing (Weigelin et al. 2021).\n")
    base = assay(0.0, seeds=seeds)[0]
    print(f"{'E:T':>6}{'live':>8}{'lysis':>8}{'kills':>8}{'hits':>8}{'hits/kill':>11}")
    rows = []
    for et in ratios:
        live, kills, hits = assay(et, seeds=seeds)
        hk = hits / kills if kills else float("nan")
        rows.append({"et": et, "live": live, "lysis": 1 - live / base,
                     "kills": kills, "hits": hits, "hits_per_kill": hk})
        print(f"{et:>6}{live:>8.0f}{1 - live / base:>8.0%}{kills:>8.0f}"
              f"{hits:>8.0f}{hk:>11.2f}")

    print("\nEffect of the repair window at E:T 1.0 (50 min is the measured value):")
    windows = []
    for mins in (15, 30, 50, 120):
        live, kills, hits = assay(1.0, decay_h=mins / 60.0, seeds=seeds)
        windows.append({"window_min": mins, "live": live, "kills": kills,
                        "lysis": 1 - live / base})
        print(f"  {mins:>4} min: lysis {1 - live / base:>5.0%}, kills {kills:>5.0f}")

    print("\nOne lethal hit per contact (the older picture this replaced):")
    single = []
    for et in (0.1, 0.3, 1.0):
        live, kills, hits = assay(et, hits_to_kill=1, seeds=seeds)
        single.append({"et": et, "live": live, "lysis": 1 - live / base, "kills": kills})
        print(f"  E:T {et}: lysis {1 - live / base:>5.0%}  (three-hit model: "
              f"{next(r['lysis'] for r in rows if r['et'] == et):.0%})")

    low = next(r for r in rows if r["et"] == 0.3)
    high = next(r for r in rows if r["et"] == 3.0)
    c1 = (np.isfinite(low["hits_per_kill"]) and np.isfinite(high["hits_per_kill"])
          and low["hits_per_kill"] > high["hits_per_kill"] * 3)
    w = {d["window_min"]: d for d in windows}
    c2 = w[15]["lysis"] <= w[50]["lysis"] <= w[120]["lysis"] + 1e-9
    c3 = all(s["lysis"] > next(r["lysis"] for r in rows if r["et"] == s["et"])
             for s in single if s["et"] <= 0.3)
    c4 = all(r["kills"] <= r["hits"] / 3.0 + 1e-6 for r in rows)

    print("\nMechanism checks:")
    for ok, text in ((c1, "hits spent per kill falls steeply with density (cooperativity)"),
                     (c2, "a shorter repair window makes killing harder, a longer one easier"),
                     (c3, "one-hit killing is more efficient at low density, as the old model was"),
                     (c4, "no kill without its three hits")):
        print(f"  [{'PASS' if ok else 'FAIL'}] {text}")
    allok = bool(c1 and c2 and c3 and c4)
    print(f"\n{'Mechanism reproduced.' if allok else 'A MECHANISM CHECK FAILED.'} "
          f"The absolute lysis at a given E:T is NOT validated:\nno measured "
          f"chromium-release curve has been digitised, so the rates stay literature\n"
          f"values and the magnitudes unchecked.")
    report = {"rows": rows, "repair_window": windows, "single_hit": single,
              "mechanism_checks_pass": allok}
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
