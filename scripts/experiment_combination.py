#!/usr/bin/env python3
"""Does the ORDER of two drugs change the outcome — and is the answer real?

Phase 2, within-line, and the question that most needs a mechanism: a
correlation over cell lines cannot answer it at all, because both arms
use the same two drugs at the same doses in the same line. Only the
sequence differs.

Two things are tested, and they came out differently.

**Antagonism (robust).** Cisplatin damages DNA, driving p53 -> p21 and
arresting the cycle; paclitaxel kills only cells that attempt mitosis
while it is present. A cytostatic therefore removes the very cells a
phase-specific partner needs, and the pair should kill less than
independent action predicts. It does, by 1.3-1.5x, in every arrangement.

**Sequence (small, dose-dependent, and mostly a readout effect).** At a
72 h readout, giving the platinum first can look better. Measured over
five seeds it is small and depends on dose: with each drug at its own
72 h IC50 there is no ordering effect at any readout (0.99-1.01); with
cisplatin at 1.25x its IC50 the platinum-first arm leaves ~12 % fewer
cells at 72 h (0.88 +- 0.06), and the gap shrinks to ~3 % by 120-144 h. The cause is kinetic, not
sequential — cisplatin acts slowly (adducts, then p53, then
commitment), so placing it second truncates its effect inside a fixed
window. A single 128-cell run carries noise of +-0.02-0.03 on this
ratio, which is the size of the effect itself, so single-seed readings
of it (including this script's first version) are not evidence.

Stated plainly, the engine does NOT reproduce a biological
platinum/taxane sequence dependence. It reproduces the antagonism and it
exposes the timing confound.

Nothing here is fitted: the arrest comes from the p53 -> p21 axis ported
from the C++ prototype, the mitosis requirement from the tubulin
occupancy threshold.

Usage:
    python scripts/experiment_combination.py           # ~10 min
    python scripts/experiment_combination.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale, ic50, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "combination_experiment.json"
TOTAL_H = 72.0
HALF = TOTAL_H / 2
# Doses as multiples of each drug's own 72 h IC50 on the line under test
# (computed by the engine at run time, so a refit of the library cannot
# change what the design means): active alone over 36 h, far from
# saturating. The second design raises cisplatin, which is where an
# ordering effect appears at all.
DESIGNS = {"both at 1x IC50": (1.0, 1.0), "cisplatin 1.25x, paclitaxel 1x": (1.25, 1.0)}
SEEDS = (1, 2, 3, 4, 5)
READOUTS_H = (72.0, 96.0, 120.0, 144.0)


def _trajectory(line, drugs, dose_fn, *, n_cells, k, p, seed):
    """Alive weight at each readout in READOUTS_H (one run to the last)."""
    res = simulate(line, drugs, 0.0, t_end_h=READOUTS_H[-1], n_cells=n_cells, k_cyc=k, p=p,
                   dose_fn=dose_fn, record_every_h=24.0, seed=seed)
    idx = [int(round(t / 24.0)) for t in READOUTS_H]
    return res.alive_weight[idx, 0]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", default="A549")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    n_cells = 96 if a.quick else 256
    line = get_line(a.line)
    p = Params()
    k = calibrate_cycle_scale(line, p)
    cis, pac = get_drug("cisplatin"), get_drug("paclitaxel")
    pair = [cis, pac]
    ic_cis = ic50(line, cis, guess_uM=10.0, k_cyc=k, p=p)
    ic_pac = ic50(line, pac, guess_uM=0.03, k_cyc=k, p=p)
    print(f"{line.name}: engine 72 h IC50 cisplatin {ic_cis:.4g} uM, paclitaxel {ic_pac:.4g} uM.")
    print(f"Each drug is given for 36 h at the same concentration in every arm; only\n"
          f"WHEN differs. {len(SEEDS)} seeds x {n_cells} cells; mean +- sd over seeds.\n")
    report = {"line": line.name, "n_cells": n_cells, "seeds": list(SEEDS),
              "readouts_h": list(READOUTS_H), "ic50_uM": {"cisplatin": ic_cis, "paclitaxel": ic_pac},
              "designs": {}}
    for design, (xc, xp) in DESIGNS.items():
        CIS, PAC = xc * ic_cis, xp * ic_pac
        arms = {
            "cisplatin alone":         ([cis], lambda t: (CIS if t < HALF else 0.0,)),
            "paclitaxel alone":        ([pac], lambda t: (PAC if t < HALF else 0.0,)),
            "simultaneous":            (pair, lambda t: (CIS, PAC) if t < HALF else (0.0, 0.0)),
            "cisplatin -> paclitaxel": (pair, lambda t: (CIS, 0.0) if t < HALF else
                                        ((0.0, PAC) if t < TOTAL_H else (0.0, 0.0))),
            "paclitaxel -> cisplatin": (pair, lambda t: (0.0, PAC) if t < HALF else
                                        ((CIS, 0.0) if t < TOTAL_H else (0.0, 0.0))),
        }
        surv = {label: [] for label in arms}
        for seed in SEEDS:
            ctrl = _trajectory(line, None, None, n_cells=n_cells, k=k, p=p, seed=seed)
            for label, (drugs, fn) in arms.items():
                tr = _trajectory(line, drugs, fn, n_cells=n_cells, k=k, p=p, seed=seed)
                surv[label].append((tr / np.maximum(ctrl, 1e-12)).tolist())
        S = {label: np.array(v) for label, v in surv.items()}      # (seeds, readouts)
        indep = S["cisplatin alone"] * S["paclitaxel alone"]
        antag = {label: S[label] / indep for label in
                 ("simultaneous", "cisplatin -> paclitaxel", "paclitaxel -> cisplatin")}
        order = S["cisplatin -> paclitaxel"] / S["paclitaxel -> cisplatin"]
        print(f"{design}: cisplatin {CIS:.4g} uM, paclitaxel {PAC:.4g} uM")
        print("  combination / independent expectation at 72 h (>1 = antagonism):")
        for label, r in antag.items():
            print(f"    {label:<26}{r[:, 0].mean():>6.2f} +- {r[:, 0].std(ddof=1):.2f}")
        print("  ordering ratio, cisplatin-first / paclitaxel-first:")
        print("    readout h " + "".join(f"{t:>14.0f}" for t in READOUTS_H))
        print("    ratio     " + "".join(f"{m:>8.3f} +-{sd:.3f}" for m, sd in
                                         zip(order.mean(0), order.std(0, ddof=1))))
        print()
        report["designs"][design] = {
            "dose_x_ic50": {"cisplatin": xc, "paclitaxel": xp},
            "conc_uM": {"cisplatin": CIS, "paclitaxel": PAC},
            "surviving_fraction": {label: v for label, v in surv.items()},
            "antagonism_72h_mean": {label: float(r[:, 0].mean()) for label, r in antag.items()},
            "antagonism_72h_sd": {label: float(r[:, 0].std(ddof=1)) for label, r in antag.items()},
            "ordering_ratio_mean": order.mean(0).tolist(),
            "ordering_ratio_sd": order.std(0, ddof=1).tolist()}
    print("Antagonism is far outside seed noise in every arrangement, and grows with\n"
          "the cisplatin dose. The ordering effect is absent with both drugs at their\n"
          "IC50, appears only when cisplatin is above it, and fades as the readout\n"
          "lengthens: mostly a timing effect of a slow drug placed second.")

    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
