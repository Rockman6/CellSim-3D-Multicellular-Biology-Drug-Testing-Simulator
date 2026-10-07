#!/usr/bin/env python3
"""Does lowering medium glucose do to a simulated spheroid what it does
to a real one?

WHAT CAN AND CANNOT BE CHECKED HERE. Mueller-Klieser, Freyer and
Sutherland grew EMT6/Ro spheroids in 0.8-16.5 mM glucose at 5 % and 20 %
oxygen and measured the viable rim (Br J Cancer 1986; Cancer Res 46:3504
and 46:3513). Their quantitative tables are images in the scanned
originals and could not be extracted, so the rim thickness in micrometres
at each glucose concentration is NOT available to compare against. What
the abstracts state unambiguously, and what this script therefore tests,
is the set of DIRECTIONAL claims:

  1. lowering medium glucose substantially reduces the viable rim;
  2. restricting oxygen reduces it too, at a given glucose;
  3. necrosis can develop from a lack of either.

A direction is a weaker claim than a number, and a model that gets the
direction right can still have the magnitude badly wrong. This script is
therefore a FALSIFICATION test, not a calibration: it can show the
glucose coupling is backwards or absent, and it cannot show it is
quantitatively right.

THE SPHEROID HAS TO BE BIG ENOUGH, which the first version of this
script got wrong and is worth recording. Glucose penetrates to a depth
L = sqrt(2 D C / q): about 1170 um at 16.5 mM, 680 um at 5.5 mM and
260 um at 0.8 mM. A spheroid much smaller than L has no gradient to
speak of however right the constants are, so the first run — on 160 um
spheroids — found glucose made no difference anywhere in the published
range and "passed" only because it also tested 0.1 mM, far below what
anyone measured. At a 310 um radius the same code depletes the centre by
41 % at 0.8 mM and 2 % at 16.5 mM, which is the selectivity the papers
describe. The runs below are therefore seeded large.

A KNOWN MISSING MECHANISM, which the same papers report and this engine
does not have: as glucose falls, cellular RESPIRATION RISES — the cells
compensate by burning more oxygen. Here oxygen consumption is fixed, so
a simulated spheroid starved of glucose does not consume oxygen faster,
and its oxygen profile is unchanged. That is a real gap and the reason
the two fields here are coupled only through growth and death rather
than through each other's consumption.

Usage:
    python scripts/validate_glucose.py            # ~10 min
    python scripts/validate_glucose.py --quick    # ~3 min
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import DishParams, run_dish  # noqa: E402
from cellsim.cell.engine import Params  # noqa: E402
from cellsim.cell.library import get_line  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "glucose_validation.json"
# The concentrations Mueller-Klieser et al. used, plus the two lower ones
# where this engine's thresholds bite.
GLUCOSE_mM = (16.5, 5.5, 1.65, 0.8)
# 20 % and 5 % oxygen, as partial pressures in the medium.
O2_mmHg = {"20 %": 152.0, "5 %": 38.0}


def run(line, glucose, o2_mmHg, *, days, n_seed, grid, seeds):
    out = []
    for sd in seeds:
        dp = DishParams(geometry="spheroid", glucose_medium_mM=glucose,
                        o2_medium_mmHg=o2_mmHg)
        r = run_dish(line, t_end_h=days * 24.0, n_seed=n_seed, grid_sites=grid,
                     record_every_h=24.0, seed=int(sd), p=Params(cycle_cv=0.25), dp=dp)
        rec = r.record
        viable = float(rec.radius_um[-1] - rec.necrotic_radius_um[-1])
        out.append({"live": int(rec.n_live[-1]), "necrotic": int(rec.n_necrotic[-1]),
                    "radius_um": float(rec.radius_um[-1]),
                    "necrotic_um": float(rec.necrotic_radius_um[-1]),
                    "viable_rim_um": viable,
                    "quiescent": float(rec.quiescent_frac[-1])})
    return {k: float(np.mean([o[k] for o in out])) for k in out[0]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    # Seeded large on purpose: see the note above about penetration depth.
    days, n_seed, grid = (4.0, 30000, 90) if a.quick else (6.0, 40000, 100)
    seeds = (1,) if a.quick else (1, 2)
    line = get_line("DLD-1")

    print(f"DLD-1 spheroid, {days:g} days from {n_seed} cells, mean of {len(seeds)} seeds.")
    print("Testing the DIRECTION of three published claims; the magnitudes in the\n"
          "1986 papers are image-only tables and are not available to compare against.\n")
    report = {"days": days, "n_seed": n_seed, "seeds": list(seeds), "rows": []}
    for o2_label, o2 in O2_mmHg.items():
        print(f"  oxygen {o2_label} ({o2:.0f} mmHg)")
        print(f"{'glucose':>10}{'live':>8}{'radius':>9}{'viable rim':>12}"
              f"{'necrotic core':>15}{'quiescent':>11}")
        for g in GLUCOSE_mM:
            r = run(line, g, o2, days=days, n_seed=n_seed, grid=grid, seeds=seeds)
            r.update({"glucose_mM": g, "o2_mmHg": o2, "o2_label": o2_label})
            report["rows"].append(r)
            print(f"{g:>9.2f} {r['live']:>8.0f}{r['radius_um']:>9.0f}"
                  f"{r['viable_rim_um']:>12.0f}{r['necrotic_um']:>15.0f}"
                  f"{r['quiescent']:>11.2f}")
        print()

    rows = report["rows"]

    def at(g, lab):
        return next(r for r in rows if r["glucose_mM"] == g and r["o2_label"] == lab)

    # Claim 1, tested WITHIN the published range (16.5 -> 0.8 mM), not by
    # reaching for a concentration below anything anyone measured.
    hi, lo = at(16.5, "20 %"), at(0.8, "20 %")
    c1 = lo["viable_rim_um"] < hi["viable_rim_um"] * 0.95
    # Claim 2: restricting oxygen reduces it too, at a given glucose.
    c2 = at(5.5, "5 %")["viable_rim_um"] < at(5.5, "20 %")["viable_rim_um"] * 0.98
    # Claim 3: necrosis develops from a lack of either.
    c3_glc = at(0.8, "20 %")["necrotic"] > at(16.5, "20 %")["necrotic"]
    c3_o2 = at(16.5, "5 %")["necrotic"] > at(16.5, "20 %")["necrotic"]
    # Monotonicity: the response should not reverse as glucose falls.
    rim20 = [at(g, "20 %")["viable_rim_um"] for g in GLUCOSE_mM]
    c4 = all(rim20[i] >= rim20[i + 1] * 0.98 for i in range(len(rim20) - 1))

    print("Directional claims (Mueller-Klieser, Freyer & Sutherland 1986):")
    for ok, text in ((c1, "lowering glucose reduces growth / the viable rim"),
                     (c2, "restricting oxygen reduces it too, at a given glucose"),
                     (c3_glc, "lack of GLUCOSE alone can produce necrosis"),
                     (c3_o2, "lack of OXYGEN alone can produce necrosis"),
                     (c4, "the response is monotone in glucose (no reversal)")):
        print(f"  [{'PASS' if ok else 'FAIL'}] {text}")
    allok = all((c1, c2, c3_glc, c3_o2, c4))
    report["directional_claims_pass"] = bool(allok)
    print(f"\n{'All directions reproduced.' if allok else 'A DIRECTION IS WRONG.'} This is "
          f"falsification, not calibration:\nthe magnitudes remain unchecked, and the "
          f"missing respiration-compensation\nmechanism (cells burn more oxygen as "
          f"glucose falls) means the oxygen profile\nof a starved spheroid is wrong here "
          f"whatever the glucose numbers do.")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
