#!/usr/bin/env python3
"""Can any per-line input move the predicted IC50 enough to matter?

The Phase-1 exit gate fails because the engine ties a constant-IC50
null at telling cell lines apart. This asks why, by sweeping each
per-line field across the full range the curated lines actually span
and measuring how far the predicted IC50 moves, against how far the
real IC50s differ between lines.

A field whose sweep moves the IC50 by 1.0x cannot explain a 10x spread
no matter how accurately it is measured, so this says where per-line
information has to enter and where it would be wasted.

Usage:
    python scripts/experiment_leverage.py            # ~25 min
    python scripts/experiment_leverage.py --quick    # ~8 min
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import calibrate_cycle_scale, dose_response, ic50_from_curve  # noqa: E402
from cellsim.cell.library import CELL_LINES, get_drug  # noqa: E402

REFERENCE_CSV = REPO_ROOT / "benchmarks" / "cell" / "gdsc_reference.csv"
OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "leverage_experiment.json"

# Concentration grids wide enough that the IC50 stays bracketed as the
# knob is swept, so a moved IC50 is measured rather than clipped.
GRIDS = {"cisplatin": (-1, 3), "doxorubicin": (-3, 2), "paclitaxel": (-4, 1)}

# (field, values, how the range was chosen)
KNOBS = [
    ("doubling_time_h", (18.0, 28.0, 38.0),
     "full range of the curated lines, 18-38 h"),
    ("bcl2_level", (0.5, 1.0, 2.0, 4.0),
     "8-fold span of anti-apoptotic reserve"),
    ("repair_rate_per_h", (math.log(2) / 32, math.log(2) / 8, math.log(2) / 2),
     "repair half-life 32 h to 2 h"),
]


def observed_spread() -> dict[str, float]:
    """Max/min in-range GDSC IC50 per drug — what a model must explain."""
    vals: dict[str, list[float]] = {}
    with REFERENCE_CSV.open(newline="") as f:
        for r in csv.DictReader(f):
            if r["ic50_in_tested_range"] == "True":
                vals.setdefault(r["drug"], []).append(float(r["ic50_uM"]))
    return {d: max(v) / min(v) for d, v in vals.items() if len(v) > 1}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    n_cells, n_conc = (16, 13) if a.quick else (24, 17)

    base = CELL_LINES["A549"]
    obs = observed_spread()
    report: dict = {"n_cells_per_dose": n_cells, "observed_spread": obs, "knobs": {}}

    print(f"{'per-line field':<20}{'drug':<13}{'IC50 span over the sweep':>26}"
          f"{'needed':>9}   predicted IC50 (uM) at each value")
    for knob, values, why in KNOBS:
        report["knobs"][knob] = {"values": list(values), "range_rationale": why, "drugs": {}}
        for drug_name, (lo, hi) in GRIDS.items():
            conc = np.logspace(lo, hi, n_conc)
            ics = []
            for v in values:
                line = dataclasses.replace(base, name="probe", **{knob: v})
                viab, _ = dose_response(line, get_drug(drug_name), conc,
                                        n_cells_per_conc=n_cells,
                                        k_cyc=calibrate_cycle_scale(line))
                ic = ic50_from_curve(conc, viab[1:])
                ics.append(float(ic) if np.isfinite(ic) else None)
            good = [x for x in ics if x]
            span = (max(good) / min(good)) if len(good) > 1 else float("nan")
            report["knobs"][knob]["drugs"][drug_name] = {
                "ic50_uM": ics, "span": span, "observed_span": obs.get(drug_name)}
            shown = " ".join(f"{x:.3g}" if x else "n/a" for x in ics)
            print(f"{knob:<20}{drug_name:<13}{span:>24.1f}x"
                  f"{obs.get(drug_name, float('nan')):>8.1f}x   {shown}")

    print("\nReading: a field whose sweep moves the IC50 far less than lines actually\n"
          "differ cannot explain line-to-line variation, however well it is measured.")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
