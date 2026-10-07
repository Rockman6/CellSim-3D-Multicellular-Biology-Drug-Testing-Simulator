#!/usr/bin/env python3
"""A held-out test of HeLa's cycle time, on a movie never used to fit it.

HeLa's 19 h cycle and 50 % quiescent fraction were calibrated on the
Cell Tracking Challenge's Fluo-N2DL-HeLa sequences (scripts/
validate_dish.py), so agreement with those is not evidence. This uses
the CTC's other HeLa dataset, DIC-C2DH-HeLa: different cultures,
different microscopy, 84 frames of 10 min (14 h).

Fourteen hours is too short to see a daughter complete a cycle, so it
cannot test the quiescent fraction. It CAN test the cycle time, through
the fraction of cells present at the start that divide within the
window: a 19 h cycle and the old 31 h one disagree about it.

Cells present at the start are never quiescent in the engine (that is
decided at mitotic exit), so this checks the cycle time ALONE.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import DishParams, run_dish  # noqa: E402
from cellsim.cell.engine import Params  # noqa: E402
from cellsim.cell.library import get_line  # noqa: E402

ZIP = REPO_ROOT / "data" / "ctc" / "DIC-C2DH-HeLa.zip"
URL = "https://data.celltrackingchallenge.net/training-datasets/DIC-C2DH-HeLa.zip"
OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "hela_heldout.json"
WINDOW_H = 14.0


def measured() -> tuple[int, int]:
    """Founders that divided, and founders observed for the whole movie.
    One that leaves the field before dividing is censored, not counted
    as a non-divider."""
    if not ZIP.exists():
        ZIP.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {URL} (42 MB) ...", flush=True)
        urllib.request.urlretrieve(URL, ZIP)
    z = zipfile.ZipFile(ZIP)
    divided = observed = 0
    for seq in ("01", "02"):
        rows = [list(map(int, l.split())) for l in
                z.read(f"DIC-C2DH-HeLa/{seq}_GT/TRA/man_track.txt").decode().strip().split("\n")]
        founders = {r[0] for r in rows if r[1] == 0 and r[3] == 0}
        parents = {r[3] for r in rows if r[3] != 0}
        last = max(r[2] for r in rows)
        censored = {r[0] for r in rows
                    if r[0] in founders and r[0] not in parents and r[2] < last}
        obs = founders - censored
        divided += len(obs & parents)
        observed += len(obs)
    return divided, observed


def simulated(line, runs: int) -> np.ndarray:
    out = []
    for seed in range(1, runs + 1):
        r = run_dish(line, t_end_h=WINDOW_H, n_seed=10, grid_sites=(40, 40),
                     record_every_h=WINDOW_H, seed=seed, p=Params(cycle_cv=0.25),
                     dp=DishParams(geometry="monolayer", spacing_um=20.0,
                                   field_dt_h=0.5, dt_h=0.02))
        # with no deaths and at most one division each in 14 h, the
        # number of extra cells is the number of founders that divided
        out.append(min(max((r.record.n_live[-1] - 10) / 10.0, 0.0), 1.0))
    return np.array(out)


def p_at_least(k: int, n: int, p: float) -> float:
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=int, default=20)
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    k, n = measured()
    print(f"DIC-C2DH-HeLa: {k} of {n} founders divided within {WINDOW_H:g} h "
          f"({k / n:.2f}). Never used for calibration.\n")
    hela = get_line("HeLa")
    cases = {"current (19 h cycle)": hela,
             "old (31 h cycle)": dataclasses.replace(hela, doubling_time_h=31.0,
                                                     quiescent_fraction=0.0)}
    report = {"measured_divided": k, "measured_observed": n, "window_h": WINDOW_H,
              "models": {}}
    for label, line in cases.items():
        f = simulated(line, a.runs)
        p = p_at_least(k, n, float(f.mean()))
        verdict = "rejected" if p < 0.05 else "consistent"
        report["models"][label] = {"mean": float(f.mean()), "sd": float(f.std()),
                                   "p_at_least_measured": p, "verdict": verdict}
        print(f"  {label:<22} {f.mean():.2f} +- {f.std():.2f}   "
              f"P(>= {k}/{n}) = {p:.3f}  -> {verdict}")
    cur, old = report["models"]["current (19 h cycle)"], report["models"]["old (31 h cycle)"]
    ok = cur["verdict"] == "consistent" and old["verdict"] == "rejected"
    report["current_supported_old_rejected"] = ok
    print(f"\n{'Held-out data supports the 19 h cycle and rejects 31 h.' if ok else 'NOT the expected split.'}")
    print("This checks the cycle time only; the 50 % quiescent fraction is untested here.")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1) + "\n")
        print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
