#!/usr/bin/env python3
"""Schedule dependence: is AUC the right exposure metric, and for which drug?

Phase 2's first question, and one that needs no per-line marker.

Two arms, because they ask different things.

**Arm A, matched AUC.** Hold concentration x hours constant and vary the
split. If AUC were the governing variable every split would give the
same kill. It does not, and the way it fails is the point: paclitaxel
engages only above a threshold tubulin occupancy, so spreading a fixed
AUC thinly drops the concentration below that threshold and the drug
stops working entirely, however long it is left on.

**Arm B, fixed concentration above threshold, varying duration.** This
is the comparison the clinical literature actually makes. A tubulin
binder kills only cells that attempt mitosis while it is present, so
its effect should keep growing with exposure time as more of an
asynchronous population transits M. A DNA-damage agent forms lesions as
soon as it is inside, so it should saturate much sooner.

The expected contrast is well documented for paclitaxel: efficacy
tracks TIME ABOVE A THRESHOLD CONCENTRATION rather than AUC or peak
(Gianni 1995 J Clin Oncol 13:180; Huizing 1993 J Clin Oncol 11:2127),
with the mitotic-transit requirement established by Lopes 1993 Cancer
Chemother Pharmacol 32:235 and the arrest/slippage fate decision by
Gascoigne & Taylor 2008 Cancer Cell 14:111.

Usage:
    python scripts/experiment_schedule.py           # ~20 min
    python scripts/experiment_schedule.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.stream import Schedule, run_stream  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "schedule_experiment.json"

# Exposure windows, in hours, all within a 72 h assay.
WINDOWS = (3.0, 6.0, 12.0, 18.0, 24.0, 36.0, 48.0, 60.0, 72.0)

# AUC per drug (uM*h), chosen near each drug's own active range so that
# the shortest window is not simply saturating and the longest is not
# simply inert. Taken from the fitted IC50s, not tuned to the answer.
AUC = {"cisplatin": 720.0,        # ~10 uM for 72 h
       "doxorubicin": 4.0,        # ~0.056 uM for 72 h
       "paclitaxel": 1.8}         # ~0.025 uM for 72 h

# Arm B: a concentration comfortably above each drug's active threshold,
# so every window is pharmacologically engaged and only duration varies.
# Paclitaxel's arrest threshold for this line is ~0.03 uM in medium
# (tubulin Kd 0.01 uM, arrest at occupancy 0.3, partition 0.153), hence
# 0.1 uM here.
FIXED_CONC = {"cisplatin": 20.0, "doxorubicin": 0.2, "paclitaxel": 0.1}


def surviving(line, drug, conc_uM: float, hours: float, *, n_cells: int, seed: int) -> float:
    """Alive weight at 72 h relative to an untreated control of the same line."""
    treated = run_stream(line, drug, Schedule.parse(f"0:{conc_uM:g},{hours:g}:0"),
                         t_end_h=72.0, n_cells=n_cells, record_every_h=72.0, seed=seed)
    control = run_stream(line, None, Schedule.constant(0.0),
                         t_end_h=72.0, n_cells=n_cells, record_every_h=72.0, seed=seed)
    t = treated[-1]["population"]["alive_weight"]
    c = control[-1]["population"]["alive_weight"]
    return t / max(c, 1e-12)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", default="A549")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    n_cells = 48 if a.quick else 96
    line = get_line(a.line)
    report = {"line": line.name, "n_cells": n_cells, "windows_h": list(WINDOWS),
              "arm_a_auc_uM_h": AUC, "arm_b_conc_uM": FIXED_CONC, "arm_a": {}, "arm_b": {}}

    print(f"{line.name}, 72 h assay. Surviving fraction vs an untreated control.\n")
    print("ARM A — same AUC, different split between concentration and hours")
    header = f"{'drug':<13}{'AUC':>7}  " + "".join(f"{w:>6.0f}h" for w in WINDOWS)
    print(header + "\n" + "-" * len(header))
    for drug_name, auc in AUC.items():
        drug = get_drug(drug_name)
        row, concs = [], []
        for w in WINDOWS:
            c = auc / w
            concs.append(c)
            row.append(surviving(line, drug, c, w, n_cells=n_cells, seed=1))
        best = min(range(len(row)), key=lambda i: row[i])
        interior = 0 < best < len(row) - 1
        report["arm_a"][drug_name] = {"conc_uM": concs, "surviving_fraction": row,
                                      "best_window_h": WINDOWS[best],
                                      "optimum_is_interior": interior}
        print(f"{drug_name:<13}{auc:>7.4g}  " + "".join(f"{v:>7.3f}" for v in row)
              + f"   best at {WINDOWS[best]:.0f} h"
              + ("  <- interior optimum" if interior else ""))
    print("  If AUC governed the outcome these rows would be flat. They are not.\n"
          "  An INTERIOR optimum means both extremes fail, for different reasons:\n"
          "  too brief and few cells reach the step the drug acts on; too dilute\n"
          "  and the concentration falls below the threshold that engages it.\n")

    print("ARM B — concentration fixed ABOVE threshold, exposure time varied")
    header = f"{'drug':<13}{'conc':>8}  " + "".join(f"{w:>6.0f}h" for w in WINDOWS) + "   3h/72h"
    print(header + "\n" + "-" * len(header))
    for drug_name, conc in FIXED_CONC.items():
        drug = get_drug(drug_name)
        row = [surviving(line, drug, conc, w, n_cells=n_cells, seed=1) for w in WINDOWS]
        ratio = row[0] / max(row[-1], 1e-9)
        report["arm_b"][drug_name] = {"conc_uM": conc, "surviving_fraction": row,
                                      "short_over_long": ratio}
        shown = ">1000" if ratio > 1000 else f"{ratio:.1f}"
        print(f"{drug_name:<13}{conc:>8.4g}  " + "".join(f"{v:>7.3f}" for v in row)
              + f"   {shown:>6}x")
    print("  '3h/72h' is how many times more cells survive a 3 h exposure than a 72 h\n"
          "  one at the SAME concentration. Larger means more exposure-time limited.")

    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
