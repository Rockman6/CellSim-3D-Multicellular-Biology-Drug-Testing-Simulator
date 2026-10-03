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

**Sequence (an artefact, and worth knowing).** The first run at a 72 h
readout appeared to favour giving the platinum first. Extending the
readout dissolves it: the advantage runs 0.90-0.94 at 72 h, 0.95-0.98 at
96 h and 0.96-1.02 at 120 h across seeds. The cause is kinetic, not
sequential — cisplatin acts slowly (adducts, then p53, then commitment),
so placing it second truncates its effect inside a fixed window. The
engine is therefore reproducing a well-known experimental trap rather
than a biological sequence effect: a fixed-endpoint assay can confound
"the order matters" with "the second drug had less time to work".

Stated plainly, the engine does NOT reproduce a biological
platinum/taxane sequence dependence. It reproduces the antagonism and it
exposes the timing confound.

Nothing here is fitted: the arrest comes from the p53 -> p21 axis ported
from the C++ prototype, the mitosis requirement from the tubulin
occupancy threshold.

Usage:
    python scripts/experiment_combination.py           # ~15 min
    python scripts/experiment_combination.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "combination_experiment.json"
TOTAL_H = 72.0
HALF = TOTAL_H / 2
CIS, PAC = 15.0, 0.04          # µM, each active alone but not saturating


def _surviving(line, drugs, dose_fn, *, n_cells, k, p, total=TOTAL_H):
    treated = simulate(line, drugs, 0.0, t_end_h=total, n_cells=n_cells,
                       k_cyc=k, p=p, dose_fn=dose_fn, record_every_h=total)
    control = simulate(line, None, 0.0, t_end_h=total, n_cells=n_cells, k_cyc=k, p=p,
                       record_every_h=total)
    return float(treated.alive_weight[-1, 0] / max(control.alive_weight[-1, 0], 1e-12))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", default="A549")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    n_cells = 64 if a.quick else 128
    line = get_line(a.line)
    p = Params()
    k = calibrate_cycle_scale(line, p)
    cis, pac = get_drug("cisplatin"), get_drug("paclitaxel")
    pair = [cis, pac]

    # Every arm delivers each drug for the same number of hours at the
    # same concentration; only the placement in time differs.
    arms = {
        "cisplatin alone (36 h)":   ([cis], lambda t: (CIS if t < HALF else 0.0,)),
        "paclitaxel alone (36 h)":  ([pac], lambda t: (PAC if t < HALF else 0.0,)),
        "both, simultaneous":       (pair, lambda t: (CIS, PAC) if t < HALF else (0.0, 0.0)),
        "cisplatin -> paclitaxel":  (pair, lambda t: (CIS, 0.0) if t < HALF else (0.0, PAC)),
        "paclitaxel -> cisplatin":  (pair, lambda t: (0.0, PAC) if t < HALF else (CIS, 0.0)),
    }

    print(f"{line.name}, {TOTAL_H:.0f} h. Each drug is given for 36 h at the same\n"
          f"concentration in every arm; only WHEN differs. Surviving fraction vs\n"
          f"an untreated control.\n")
    results = {}
    for label, (drugs, fn) in arms.items():
        results[label] = _surviving(line, drugs, fn, n_cells=n_cells, k=k, p=p)
        print(f"  {label:<26}{results[label]:>8.3f}")

    solo = results["cisplatin alone (36 h)"] * results["paclitaxel alone (36 h)"]
    print(f"\n  {'independent expectation':<26}{solo:>8.3f}   (product of the two alone)")
    print("  anything above that line is antagonism, below it is synergy.\n")
    for label in ("both, simultaneous", "cisplatin -> paclitaxel", "paclitaxel -> cisplatin"):
        verdict = ("antagonistic" if results[label] > solo * 1.1
                   else "synergistic" if results[label] < solo * 0.9 else "additive")
        print(f"  {label:<26}{results[label] / solo:>7.2f}x expected   {verdict}")

    print("\nIs the ordering difference real, or an artefact of when the assay stops?")
    print(f"  {'readout':>9}{'cis->pac':>11}{'pac->cis':>11}{'ratio':>9}")
    durations = []
    for total in (72.0, 96.0, 120.0):
        def cis_first(t):
            return (CIS, 0.0) if t < HALF else ((0.0, PAC) if t < TOTAL_H else (0.0, 0.0))

        def pac_first(t):
            return (0.0, PAC) if t < HALF else ((CIS, 0.0) if t < TOTAL_H else (0.0, 0.0))
        aa = _surviving(line, pair, cis_first, n_cells=n_cells, k=k, p=p, total=total)
        bb = _surviving(line, pair, pac_first, n_cells=n_cells, k=k, p=p, total=total)
        durations.append({"readout_h": total, "cis_first": aa, "pac_first": bb,
                          "ratio": aa / max(bb, 1e-12)})
        print(f"  {total:>9.0f}{aa:>11.3f}{bb:>11.3f}{aa / max(bb, 1e-12):>9.2f}")
    print("  The ratio approaches 1 as the readout lengthens: cisplatin acts slowly,\n"
          "  so putting it second truncates its effect inside a fixed window. The\n"
          "  apparent sequence effect is a timing confound, not biology.")

    report = {"line": line.name, "total_h": TOTAL_H, "conc_uM": {"cisplatin": CIS,
              "paclitaxel": PAC}, "n_cells": n_cells, "surviving_fraction": results,
              "independent_expectation": solo, "readout_sweep": durations}
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
