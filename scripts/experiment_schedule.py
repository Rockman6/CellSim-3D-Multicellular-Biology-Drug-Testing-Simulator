#!/usr/bin/env python3
"""Schedule dependence: is AUC the right exposure metric, and for which drug?

Phase 2's first question, and one that needs no per-line marker.

The measurement is the ISO-EFFECT CURVE: for each exposure time T the
concentration C50(T) that halves the 72 h cell count when the drug is
washed out at T. Everything a schedule comparison can ask is in that
one curve, and it does not depend on a chosen dose:

* If killing were governed by AUC, C50(T) * T would be constant and
  log C50 against log T would have slope -1. Ozawa et al. measured
  exactly that for cell-cycle phase-non-specific agents, cisplatin
  included (Cancer Chemother Pharmacol 1988 21:185; Cancer Res 1989
  49:3823), and time dependence instead for phase-specific ones.
* A drug that needs cells to reach a phase while it is present cannot
  reach half kill at all below some exposure time, however high the
  concentration (C50 = inf), and gains steeply with time above it.
  That is the documented taxane behaviour: efficacy follows time above
  a threshold concentration (Gianni 1995 J Clin Oncol 13:180), 24 h
  kill saturates above ~50 nM, and 24 -> 72 h exposure raises
  cytotoxicity 5-200-fold (Liebmann 1993 Br J Cancer 68:1104).
* The best way to split a fixed AUC is the T that minimises C50(T) * T.
  An interior minimum means both extremes fail, for different reasons.

Two assay-matched comparisons are reported beside it, with their caveats:

* Liebmann's plateau: paclitaxel survival after a 24 h exposure should
  stop falling once the concentration is well above threshold.
* Georgiadis et al. 1997 (Clin Cancer Res 3:449), 14 NSCLC lines, MTT:
  median paclitaxel IC50 >32 uM, 9.4 uM and 0.027 uM for 3, 24 and
  120 h exposures. Here A549 (an NSCLC line) is read at 120 h after the
  same three exposures. MTT measures metabolic activity, which
  mitotically arrested and slipped cells keep, while the engine counts
  live cells, so the comparison is directional, not quantitative.

Usage:
    python scripts/experiment_schedule.py           # ~10 min
    python scripts/experiment_schedule.py --quick
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale, dose_response, ic50  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "schedule_experiment.json"

# Exposure times, all washed out inside a 72 h assay.
EXPOSURES_H = (1.0, 2.0, 3.0, 6.0, 12.0, 24.0, 48.0, 72.0)
IC50_GUESS_UM = {"cisplatin": 10.0, "doxorubicin": 0.05, "paclitaxel": 0.03}
# Short-exposure window over which a C x T law is tested (Ozawa's range).
CXT_WINDOW_H = (1.0, 6.0)


def iso_effect(line, drug, ic72: float, *, n_cells: int, k: float, p: Params,
               t_end_h: float = 72.0, exposures=EXPOSURES_H) -> list[float]:
    """C50 for each exposure time (inf where half kill is unreachable)."""
    out = []
    for T in exposures:
        guess = ic72 * t_end_h / T if T < t_end_h else ic72
        out.append(ic50(line, drug, guess_uM=guess, exposure_h=None if T >= t_end_h else T,
                        t_end_h=t_end_h, n_cells_per_conc=n_cells, k_cyc=k, p=p))
    return out


def loglog_slope(ts, cs) -> float:
    """Least-squares slope of log C50 on log T over finite points."""
    pts = [(math.log(t), math.log(c)) for t, c in zip(ts, cs) if np.isfinite(c)]
    if len(pts) < 2:
        return float("nan")
    x, y = np.array(pts).T
    return float(np.polyfit(x, y, 1)[0])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", default="A549")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    n_cells = 24 if a.quick else 48
    line = get_line(a.line)
    p = Params()
    k = calibrate_cycle_scale(line, p)
    report = {"line": line.name, "n_cells_per_conc": n_cells, "exposures_h": list(EXPOSURES_H),
              "cxt_window_h": list(CXT_WINDOW_H), "drugs": {}}

    print(f"{line.name}. C50(T): concentration that halves the 72 h cell count when the drug\n"
          f"is washed out at T. 'inf' = no concentration reaches half kill.\n")
    for name, guess in IC50_GUESS_UM.items():
        drug = get_drug(name)
        ic72 = ic50(line, drug, guess_uM=guess, n_cells_per_conc=n_cells, k_cyc=k, p=p)
        c50 = iso_effect(line, drug, ic72, n_cells=n_cells, k=k, p=p)
        cxt = [c * t if np.isfinite(c) else float("inf") for c, t in zip(c50, EXPOSURES_H)]
        lo, hi = CXT_WINDOW_H
        win = [(t, c) for t, c in zip(EXPOSURES_H, c50) if lo <= t <= hi]
        slope = loglog_slope([t for t, _ in win], [c for _, c in win])
        finite = [t for t, c in zip(EXPOSURES_H, c50) if np.isfinite(c)]
        t_min = min(finite) if finite else float("inf")
        best = int(np.argmin(cxt))
        report["drugs"][name] = {
            "ic50_72h_uM": ic72, "c50_uM": c50, "c50_x_T_uM_h": cxt,
            "short_exposure_loglog_slope": slope, "shortest_exposure_reaching_half_kill_h": t_min,
            "best_split_h": EXPOSURES_H[best],
            "best_split_is_interior": 0 < best < len(EXPOSURES_H) - 1}
        print(f"{name}  (72 h continuous IC50 {ic72:.4g} uM)")
        print("   T (h)   " + "".join(f"{t:>9.0f}" for t in EXPOSURES_H))
        print("   C50 uM  " + "".join(f"{c:>9.3g}" if np.isfinite(c) else f"{'inf':>9}" for c in c50))
        print("   C50*T   " + "".join(f"{v:>9.3g}" if np.isfinite(v) else f"{'inf':>9}" for v in cxt))
        s_txt = f"{slope:+.2f}" if np.isfinite(slope) else "n/a (no half kill)"
        print(f"   slope of log C50 on log T over {lo:g}-{hi:g} h: {s_txt}   (-1 = AUC-governed)")
        print(f"   shortest exposure that can reach half kill: {t_min:g} h;"
              f" cheapest split of a fixed AUC: {EXPOSURES_H[best]:g} h"
              + ("  <- interior" if 0 < best < len(EXPOSURES_H) - 1 else "") + "\n")

    # Liebmann 1993: 24 h paclitaxel kill saturates above ~50 nM.
    pac = get_drug("paclitaxel")
    plateau_c = np.array([0.05, 0.2, 1.0, 5.0])
    viab, _ = dose_response(line, pac, plateau_c, t_end_h=72.0, n_cells_per_conc=2 * n_cells,
                            k_cyc=k, p=p, exposure_h=24.0)
    plateau = viab[1:].tolist()
    report["liebmann_24h_plateau"] = {"conc_uM": plateau_c.tolist(), "surviving_fraction": plateau,
                                      "spread": float(max(plateau) - min(plateau))}
    print("Liebmann 1993 plateau: paclitaxel 24 h exposure, read at 72 h")
    print("   conc uM " + "".join(f"{c:>8g}" for c in plateau_c))
    print("   surv    " + "".join(f"{v:>8.3f}" for v in plateau)
          + f"   spread {max(plateau) - min(plateau):.3f} (published: no extra kill above ~50 nM)\n")

    # Georgiadis 1997 protocol: 3 / 24 / 120 h exposure, read at 120 h.
    ic120 = ic50(line, pac, guess_uM=0.03, t_end_h=120.0, n_cells_per_conc=n_cells, k_cyc=k, p=p)
    geo = iso_effect(line, pac, ic120, n_cells=n_cells, k=k, p=p, t_end_h=120.0,
                     exposures=(3.0, 24.0, 120.0))
    report["georgiadis_protocol"] = {"exposures_h": [3, 24, 120], "read_at_h": 120,
                                     "engine_c50_uM": geo,
                                     "published_nsclc_median_uM": [">32", 9.4, 0.027]}
    print("Georgiadis 1997 protocol: paclitaxel, exposure 3 / 24 / 120 h, read at 120 h")
    print("   engine (A549)       " + "  ".join(f"{c:.3g}" if np.isfinite(c) else "inf" for c in geo))
    print("   NSCLC median (MTT)  >32  9.4  0.027   (MTT also counts arrested and slipped cells)")

    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
