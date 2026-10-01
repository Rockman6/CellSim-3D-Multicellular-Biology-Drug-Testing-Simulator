#!/usr/bin/env python3
"""Phase-1 validation: fit one potency constant per drug, predict every line.

Protocol (docs/PLAN.md, Phase 1 item 4):

1. Reference per (drug, line) from ``benchmarks/cell/gdsc_reference.csv``:
   the geometric mean of the GDSC screens whose IC50 lies inside their
   tested range; the acceptance span is the min-max of those screens,
   widened to at least 3x either side of the geometric mean. A line with
   no in-range screen only has a lower bound (the screen never reached
   50 % kill at its top dose); there the gate is "predicted IC50 is not
   more than 3x below that top dose".
2. Fit: for each drug, find by bisection the value of its single
   ``fit_target`` constant that reproduces the reference IC50 of each
   *fit line*, and take the geometric mean across fit lines. Fit lines
   are the clean TP53 wild-type lines (A549, MCF7) that have an in-range
   reference. Nothing else is tuned.
3. Predict every line with that constant and grade it. Lines not used in
   the fit are genuine out-of-sample tests, including both TP53-mutant
   lines and HeLa (wild-type TP53 degraded by HPV E6).
4. Ordering: Spearman rank correlation between predicted and reference
   IC50 across lines that have an in-range reference.

GDSC viability is relative cell number at 72 h, so it counts cytostasis
as well as death; the engine's readout (alive weight relative to the
untreated control at 72 h) is the same quantity.

Usage:
    python scripts/validate_gdsc.py               # ~5-10 min on a laptop
    python scripts/validate_gdsc.py --quick       # fewer cells, coarser fit
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import (  # noqa: E402
    calibrate_cycle_scale, dose_response, ic50_from_curve,
)
from cellsim.cell.library import CELL_LINES, DRUGS  # noqa: E402

REFERENCE_CSV = REPO_ROOT / "benchmarks" / "cell" / "gdsc_reference.csv"
RESULTS_CSV = REPO_ROOT / "benchmarks" / "cell" / "gdsc_validation_results.csv"
SUMMARY_JSON = REPO_ROOT / "benchmarks" / "cell" / "gdsc_validation_summary.json"
FIT_LINES = ("A549", "MCF7")   # clean TP53 wild-type
MIN_SPAN = 3.0                 # acceptance at least 3x either side of the geometric mean


# ── reference ─────────────────────────────────────────────────────────
@dataclasses.dataclass
class Target:
    kind: str            # "value" | "lower_bound"
    center: float        # geometric mean of in-range IC50s (value) or top dose (lower_bound)
    lo: float
    hi: float
    n_screens: int


def load_targets(path: Path) -> dict[tuple[str, str], Target]:
    raw: dict[tuple[str, str], dict] = {}
    with path.open(newline="") as f:
        for r in csv.DictReader(f):
            d = raw.setdefault((r["drug"], r["cell_line"]), {"in": [], "tops": []})
            d["tops"].append(float(r["max_conc_uM"]))
            if r["ic50_in_tested_range"] == "True":
                d["in"].append(float(r["ic50_uM"]))
    out = {}
    for key, d in raw.items():
        if d["in"]:
            gm = math.exp(sum(math.log(v) for v in d["in"]) / len(d["in"]))
            out[key] = Target("value", gm, min(min(d["in"]), gm / MIN_SPAN),
                              max(max(d["in"]), gm * MIN_SPAN), len(d["in"]))
        else:
            top = max(d["tops"])
            out[key] = Target("lower_bound", top, top / MIN_SPAN, math.inf, len(d["tops"]))
    return out


# ── engine wrappers ───────────────────────────────────────────────────
def predicted_ic50(line_name: str, drug, center_uM: float, *, n_cells: int,
                   k_cyc: dict[str, float], seed: int = 1) -> float:
    """IC50 from a 13-point log grid spanning 4 decades either side."""
    conc = center_uM * np.logspace(-4, 4, 13)
    viab, _ = dose_response(CELL_LINES[line_name], drug, conc,
                            n_cells_per_conc=n_cells, k_cyc=k_cyc[line_name], seed=seed)
    return ic50_from_curve(conc, viab[1:])


def with_gain(drug, gain: float):
    return dataclasses.replace(drug, **{drug.fit_target: gain})


def fit_gain_for_line(line_name: str, drug, target_uM: float, *, n_cells: int,
                      k_cyc: dict[str, float], iters: int) -> tuple[float, float]:
    """Bisection in log10 of the fitted constant. The direction (does the
    IC50 fall or rise as the constant grows?) is read off the bracket ends,
    so a potency gain and an affinity constant are fitted the same way."""
    g0 = getattr(drug, drug.fit_target)
    lo, hi = math.log10(g0) - 5, math.log10(g0) + 5
    ic_lo = predicted_ic50(line_name, with_gain(drug, 10 ** lo), target_uM, n_cells=n_cells, k_cyc=k_cyc)
    ic_hi = predicted_ic50(line_name, with_gain(drug, 10 ** hi), target_uM, n_cells=n_cells, k_cyc=k_cyc)
    falling = ic_lo >= ic_hi                     # IC50 falls as the constant rises
    weak_end, strong_end = (ic_lo, ic_hi) if falling else (ic_hi, ic_lo)
    if strong_end > target_uM:   # even the most potent setting cannot kill enough
        return math.nan, strong_end
    if weak_end < target_uM:     # even the least potent setting kills too much
        return math.nan, weak_end
    ic = math.nan
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ic = predicted_ic50(line_name, with_gain(drug, 10 ** mid), target_uM, n_cells=n_cells, k_cyc=k_cyc)
        too_weak = ic > target_uM
        if too_weak == falling:      # move toward more potency
            lo = mid
        else:
            hi = mid
    return 10 ** (0.5 * (lo + hi)), ic


def spearman(a: list[float], b: list[float]) -> float:
    def ranks(x):
        order = np.argsort(x)
        r = np.empty(len(x))
        r[order] = np.arange(len(x))
        return r
    if len(a) < 3:
        return math.nan
    ra, rb = ranks(np.asarray(a, float)), ranks(np.asarray(b, float))
    return float(np.corrcoef(ra, rb)[0, 1])


# ── main ──────────────────────────────────────────────────────────────
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="16 cells per dose, 8 bisection steps")
    ap.add_argument("--cells", type=int, default=None, help="cells per dose (default 32)")
    ap.add_argument("--drugs", default=",".join(DRUGS), help="comma-separated subset")
    ap.add_argument("--no-write", action="store_true", help="do not write the results files")
    a = ap.parse_args(argv)
    n_cells = a.cells or (16 if a.quick else 32)
    iters = 8 if a.quick else 12

    targets = load_targets(REFERENCE_CSV)
    lines = sorted({line for (_, line) in targets})
    k_cyc = {ln: calibrate_cycle_scale(CELL_LINES[ln]) for ln in lines}
    t_start = time.time()
    rows, summary = [], {"n_cells_per_dose": n_cells, "bisection_steps": iters,
                         "fit_lines": list(FIT_LINES), "drugs": {}}

    for drug_name in [d.strip() for d in a.drugs.split(",") if d.strip()]:
        drug = DRUGS[drug_name]
        drug_lines = [ln for ln in lines if (drug_name, ln) in targets]
        fit_on = [ln for ln in FIT_LINES
                  if (drug_name, ln) in targets and targets[(drug_name, ln)].kind == "value"]
        print(f"\n== {drug_name}: fit {drug.fit_target} on {fit_on or 'nothing'}")
        per_line_gain = {}
        for ln in fit_on:
            g, ic = fit_gain_for_line(ln, drug, targets[(drug_name, ln)].center,
                                      n_cells=n_cells, k_cyc=k_cyc, iters=iters)
            per_line_gain[ln] = g
            print(f"   {ln:<11} gain {g:.4g}  (IC50 at fit {ic:.3g} µM, target "
                  f"{targets[(drug_name, ln)].center:.3g})")
        good = [g for g in per_line_gain.values() if np.isfinite(g) and g > 0]
        if not good:
            print("   no fit possible; skipping predictions for this drug")
            summary["drugs"][drug_name] = {"fit_target": drug.fit_target, "gain": None,
                                           "per_line_gain": per_line_gain, "fit_on": fit_on}
            continue
        gain = math.exp(sum(math.log(g) for g in good) / len(good))
        fitted = with_gain(drug, gain)
        print(f"   fitted {drug.fit_target} = {gain:.4g}")

        ref_vals, pred_vals = [], []
        for ln in drug_lines:
            t = targets[(drug_name, ln)]
            ic = predicted_ic50(ln, fitted, t.center, n_cells=n_cells, k_cyc=k_cyc)
            ok = (t.lo <= ic <= t.hi) if t.kind == "value" else (ic >= t.lo)
            fold = (ic / t.center) if (t.kind == "value" and np.isfinite(ic)) else math.nan
            role = "fit" if ln in fit_on else "held-out"
            rows.append({"drug": drug_name, "cell_line": ln,
                         "p53_functional": CELL_LINES[ln].p53_functional, "role": role,
                         "reference_kind": t.kind, "reference_ic50_uM": round(t.center, 5),
                         "accept_lo_uM": round(t.lo, 5),
                         "accept_hi_uM": (round(t.hi, 5) if np.isfinite(t.hi) else "inf"),
                         "n_screens": t.n_screens,
                         "predicted_ic50_uM": (round(ic, 5) if np.isfinite(ic) else "inf"),
                         "fold_error": (round(fold, 3) if np.isfinite(fold) else ""),
                         "pass": bool(ok)})
            if t.kind == "value":
                ref_vals.append(t.center)
                pred_vals.append(ic if np.isfinite(ic) else 1e9)
            mark = "PASS" if ok else "MISS"
            span = (f"[{t.lo:.3g}, {t.hi:.3g}]" if t.kind == "value" else f">= {t.lo:.3g}")
            print(f"   {mark}  {ln:<11} {role:<8} ref {t.center:>8.3g} µM  span {span:<20} "
                  f"pred {ic:>8.3g} µM")
        rho = spearman(ref_vals, pred_vals)
        # Null model: one IC50 for every line, the geometric mean of the fit
        # lines' references. If the engine cannot beat this, it is not yet
        # telling lines apart; it is only reproducing the drug's potency scale.
        null_ic = math.exp(sum(math.log(targets[(drug_name, ln)].center) for ln in fit_on)
                           / len(fit_on))
        drug_rows = [r for r in rows if r["drug"] == drug_name]
        null_pass = 0
        for r in drug_rows:
            t = targets[(drug_name, r["cell_line"])]
            ok_null = (t.lo <= null_ic <= t.hi) if t.kind == "value" else (null_ic >= t.lo)
            r["null_pass"] = bool(ok_null)
            null_pass += ok_null
        rmse_model = (math.sqrt(sum(math.log10(p_ / r_) ** 2 for p_, r_ in zip(pred_vals, ref_vals))
                                / len(ref_vals)) if ref_vals else math.nan)
        rmse_null = (math.sqrt(sum(math.log10(null_ic / r_) ** 2 for r_ in ref_vals)
                               / len(ref_vals)) if ref_vals else math.nan)
        summary["drugs"][drug_name] = {
            "fit_target": drug.fit_target, "gain": gain, "per_line_gain": per_line_gain,
            "fit_on": fit_on, "spearman_vs_reference": rho,
            "log10_rmse_model": rmse_model, "log10_rmse_null": rmse_null,
            "null_ic50_uM": null_ic,
            "pass": sum(r["pass"] for r in drug_rows), "null_pass": null_pass,
            "n": len(drug_rows)}
        print(f"   Spearman (lines with in-range reference) = {rho:.2f}")
        print(f"   log10 RMSE: engine {rmse_model:.2f} vs constant-IC50 null {rmse_null:.2f}"
              f"   span passes: engine {summary['drugs'][drug_name]['pass']}/{len(drug_rows)}, "
              f"null {null_pass}/{len(drug_rows)}")

    held = [r for r in rows if r["role"] == "held-out"]
    print(f"\nheld-out predictions inside the GDSC span: engine "
          f"{sum(r['pass'] for r in held)}/{len(held)}, constant-IC50 null "
          f"{sum(r.get('null_pass', False) for r in held)}/{len(held)}   "
          f"(all rows {sum(r['pass'] for r in rows)}/{len(rows)}, {time.time() - t_start:.0f} s)")
    print("The span gate is at least 9x wide, so it checks the potency scale, not whether "
          "the engine tells lines apart; compare the log10 RMSE against the null for that.")
    summary["held_out_pass"] = sum(r["pass"] for r in held)
    summary["held_out_null_pass"] = sum(r.get("null_pass", False) for r in held)
    summary["held_out_n"] = len(held)
    if not a.no_write and rows:
        with RESULTS_CSV.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=float) + "\n")
        print(f"wrote {RESULTS_CSV.relative_to(REPO_ROOT)} and "
              f"{SUMMARY_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
