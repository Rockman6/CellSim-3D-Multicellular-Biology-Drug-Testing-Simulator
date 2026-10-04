#!/usr/bin/env python3
"""GDSC validation: fit one potency constant per drug, predict every line.

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
    Params, calibrate_cycle_scale, dose_response, ic50_from_curve,
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
    top: float = math.inf  # highest concentration any screen tested


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
        top = max(d["tops"])
        if d["in"]:
            gm = math.exp(sum(math.log(v) for v in d["in"]) / len(d["in"]))
            out[key] = Target("value", gm, min(min(d["in"]), gm / MIN_SPAN),
                              max(max(d["in"]), gm * MIN_SPAN), len(d["in"]), top)
        else:
            out[key] = Target("lower_bound", top, top / MIN_SPAN, math.inf, len(d["tops"]), top)
    return out


# ── engine wrappers ───────────────────────────────────────────────────
def predicted_ic50(line_name: str, drug, center_uM: float, *, n_cells: int,
                   k_cyc: dict[str, float], seed: int = 1, params: Params = Params()) -> float:
    """IC50 in two passes: a 13-point log grid over 4 decades either side
    to bracket the crossing, then 9 points across the bracket. The single
    coarse pass this replaced (4.6-fold steps) moved a fitted line's own
    prediction by up to 15 % between runs."""
    line = CELL_LINES[line_name]
    conc = center_uM * np.logspace(-4, 4, 13)
    viab, _ = dose_response(line, drug, conc, n_cells_per_conc=n_cells,
                            k_cyc=k_cyc[line_name], seed=seed, p=params)
    first = ic50_from_curve(conc, viab[1:])
    if not np.isfinite(first) or first <= conc[0]:
        return first
    i = int(np.searchsorted(conc, first))
    fine = np.geomspace(conc[max(0, i - 1)], conc[min(len(conc) - 1, i)], 9)
    viab, _ = dose_response(line, drug, fine, n_cells_per_conc=n_cells,
                            k_cyc=k_cyc[line_name], seed=seed, p=params)
    refined = ic50_from_curve(fine, viab[1:])
    return refined if np.isfinite(refined) else first


def with_gain(drug, gain: float):
    return dataclasses.replace(drug, **{drug.fit_target: gain})


def fit_gain_for_line(line_name: str, drug, target_uM: float, *, n_cells: int,
                      k_cyc: dict[str, float], iters: int,
                      params: Params = Params()) -> tuple[float, float]:
    """Bisection in log10 of the fitted constant. The direction (does the
    IC50 fall or rise as the constant grows?) is read off the bracket ends,
    so a potency gain and an affinity constant are fitted the same way."""
    g0 = getattr(drug, drug.fit_target)
    lo, hi = math.log10(g0) - 5, math.log10(g0) + 5
    ic_lo = predicted_ic50(line_name, with_gain(drug, 10 ** lo), target_uM, n_cells=n_cells, k_cyc=k_cyc,
                           params=params)
    ic_hi = predicted_ic50(line_name, with_gain(drug, 10 ** hi), target_uM, n_cells=n_cells, k_cyc=k_cyc,
                           params=params)
    falling = ic_lo >= ic_hi                     # IC50 falls as the constant rises
    weak_end, strong_end = (ic_lo, ic_hi) if falling else (ic_hi, ic_lo)
    if strong_end > target_uM:   # even the most potent setting cannot kill enough
        return math.nan, strong_end
    if weak_end < target_uM:     # even the least potent setting kills too much
        return math.nan, weak_end
    ic = math.nan
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ic = predicted_ic50(line_name, with_gain(drug, 10 ** mid), target_uM, n_cells=n_cells, k_cyc=k_cyc,
                           params=params)
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


# ── validation run ────────────────────────────────────────────────────
def sign_test_p(wins: int, losses: int) -> float:
    """Two-sided exact sign test. Comparing two RMSE numbers hides whether
    the difference is noise; this asks the paired question line by line."""
    n = wins + losses
    if n == 0:
        return float("nan")
    k = min(wins, losses)
    tail = sum(math.comb(n, i) for i in range(0, k + 1))
    return min(1.0, 2 * tail / 2 ** n)


def _rmse(pred: list[float], ref: list[float]) -> float:
    if not ref:
        return math.nan
    return math.sqrt(sum(math.log10(p_ / r_) ** 2 for p_, r_ in zip(pred, ref)) / len(ref))


def run_validation(*, params: Params = Params(), n_cells: int = 32, iters: int = 12,
                   drug_names: tuple[str, ...] = tuple(DRUGS), verbose: bool = True
                   ) -> tuple[list[dict], dict]:
    """Fit on FIT_LINES, predict every line; return (rows, summary)."""
    say = print if verbose else (lambda *a, **k: None)
    targets = load_targets(REFERENCE_CSV)
    lines = sorted({line for (_, line) in targets})
    k_cyc = {ln: calibrate_cycle_scale(CELL_LINES[ln], params) for ln in lines}
    t_start = time.time()
    rows: list[dict] = []
    summary: dict = {"n_cells_per_dose": n_cells, "bisection_steps": iters,
                     "fit_lines": list(FIT_LINES), "params": dataclasses.asdict(params),
                     "drugs": {}}

    for drug_name in drug_names:
        drug = DRUGS[drug_name]
        drug_lines = [ln for ln in lines if (drug_name, ln) in targets]
        fit_on = [ln for ln in FIT_LINES
                  if (drug_name, ln) in targets and targets[(drug_name, ln)].kind == "value"]
        say(f"\n== {drug_name}: fit {drug.fit_target} on {fit_on or 'nothing'}")
        per_line_gain = {}
        for ln in fit_on:
            g, ic = fit_gain_for_line(ln, drug, targets[(drug_name, ln)].center,
                                      n_cells=n_cells, k_cyc=k_cyc, iters=iters, params=params)
            per_line_gain[ln] = g
            say(f"   {ln:<11} gain {g:.4g}  (IC50 at fit {ic:.3g} µM, target "
                f"{targets[(drug_name, ln)].center:.3g})")
        good = [g for g in per_line_gain.values() if np.isfinite(g) and g > 0]
        if not good:
            say("   no fit possible; skipping predictions for this drug")
            summary["drugs"][drug_name] = {"fit_target": drug.fit_target, "gain": None,
                                           "per_line_gain": per_line_gain, "fit_on": fit_on}
            continue
        gain = math.exp(sum(math.log(g) for g in good) / len(good))
        fitted = with_gain(drug, gain)
        say(f"   fitted {drug.fit_target} = {gain:.4g}")
        # Null model: one IC50 for every line, the geometric mean of the fit
        # lines' references. If the engine cannot beat this, it is not yet
        # telling lines apart; it is only reproducing the drug's potency scale.
        null_ic = math.exp(sum(math.log(targets[(drug_name, ln)].center) for ln in fit_on)
                           / len(fit_on))

        drug_rows = []
        for ln in drug_lines:
            t = targets[(drug_name, ln)]
            ic = predicted_ic50(ln, fitted, t.center, n_cells=n_cells, k_cyc=k_cyc,
                                params=params)
            ok = (t.lo <= ic <= t.hi) if t.kind == "value" else (ic >= t.lo)
            ok_null = (t.lo <= null_ic <= t.hi) if t.kind == "value" else (null_ic >= t.lo)
            fold = (ic / t.center) if (t.kind == "value" and np.isfinite(ic)) else math.nan
            role = "fit" if ln in fit_on else "held-out"
            # Sensitive = half kill reached inside the tested range. This is
            # the readout that still counts lines the screen never killed.
            observed_sensitive = t.kind == "value"
            predicted_sensitive = bool(np.isfinite(ic) and ic <= t.top)
            null_sensitive = bool(null_ic <= t.top)
            row = {"drug": drug_name, "cell_line": ln,
                   "p53_functional": CELL_LINES[ln].p53_functional, "role": role,
                   "reference_kind": t.kind, "reference_ic50_uM": round(t.center, 5),
                   "accept_lo_uM": round(t.lo, 5),
                   "accept_hi_uM": (round(t.hi, 5) if np.isfinite(t.hi) else "inf"),
                   "n_screens": t.n_screens,
                   "predicted_ic50_uM": (round(ic, 5) if np.isfinite(ic) else "inf"),
                   "null_ic50_uM": round(null_ic, 5),
                   "fold_error": (round(fold, 3) if np.isfinite(fold) else ""),
                   "pass": bool(ok), "null_pass": bool(ok_null),
                   "top_dose_uM": t.top, "observed_sensitive": observed_sensitive,
                   "predicted_sensitive": predicted_sensitive,
                   "class_correct": predicted_sensitive == observed_sensitive,
                   "null_class_correct": null_sensitive == observed_sensitive}
            rows.append(row)
            drug_rows.append(row)
            mark = "PASS" if ok else "MISS"
            span = (f"[{t.lo:.3g}, {t.hi:.3g}]" if t.kind == "value" else f">= {t.lo:.3g}")
            say(f"   {mark}  {ln:<11} {role:<8} ref {t.center:>8.3g} µM  span {span:<20} "
                f"pred {ic:>8.3g} µM")

        def vals(rs):
            ref = [r["reference_ic50_uM"] for r in rs if r["reference_kind"] == "value"]
            pred = [(r["predicted_ic50_uM"] if r["predicted_ic50_uM"] != "inf" else 1e9)
                    for r in rs if r["reference_kind"] == "value"]
            return pred, ref

        pred_all, ref_all = vals(drug_rows)
        pred_ho, ref_ho = vals([r for r in drug_rows if r["role"] == "held-out"])
        rho = spearman(ref_all, pred_all)
        d = {"fit_target": drug.fit_target, "gain": gain, "per_line_gain": per_line_gain,
             "fit_on": fit_on, "spearman_vs_reference": rho, "null_ic50_uM": null_ic,
             "log10_rmse_model": _rmse(pred_all, ref_all),
             "log10_rmse_null": _rmse([null_ic] * len(ref_all), ref_all),
             "log10_rmse_model_heldout": _rmse(pred_ho, ref_ho),
             "log10_rmse_null_heldout": _rmse([null_ic] * len(ref_ho), ref_ho),
             "n_heldout_in_range": len(ref_ho),
             "pass": sum(r["pass"] for r in drug_rows),
             "null_pass": sum(r["null_pass"] for r in drug_rows), "n": len(drug_rows),
             "class_correct": sum(r["class_correct"] for r in drug_rows),
             "null_class_correct": sum(r["null_class_correct"] for r in drug_rows),
             "n_observed_resistant": sum(not r["observed_sensitive"] for r in drug_rows)}
        summary["drugs"][drug_name] = d
        say(f"   Spearman (lines with in-range reference) = {rho:.2f}")
        say(f"   log10 RMSE all lines: engine {d['log10_rmse_model']:.2f} vs null "
            f"{d['log10_rmse_null']:.2f};  held-out only: engine "
            f"{d['log10_rmse_model_heldout']:.2f} vs null {d['log10_rmse_null_heldout']:.2f} "
            f"(n={len(ref_ho)})   span passes: engine {d['pass']}/{d['n']}, "
            f"null {d['null_pass']}/{d['n']}")
        say(f"   sensitive-or-resistant calls: engine {d['class_correct']}/{d['n']}, null "
            f"{d['null_class_correct']}/{d['n']} ({d['n_observed_resistant']} lines never "
            f"reached half kill in the screen)")

    held = [r for r in rows if r["role"] == "held-out"]
    # Paired per-line comparison against the null. This is the exit gate:
    # comparing two RMSE numbers passed on a 0.01 difference once, which
    # was noise, so the gate now asks whether the engine beats the null on
    # MORE LINES than it loses on, and whether that could be chance.
    wins = losses = 0
    for r in held:
        if r["reference_kind"] != "value" or r["predicted_ic50_uM"] == "inf":
            continue
        ref = float(r["reference_ic50_uM"])
        e_model = abs(math.log10(float(r["predicted_ic50_uM"]) / ref))
        e_null = abs(math.log10(float(r["null_ic50_uM"]) / ref))
        if e_model < e_null:
            wins += 1
        elif e_null < e_model:
            losses += 1
    p_value = sign_test_p(wins, losses)
    discriminates = wins > losses and p_value < 0.05
    beats = [name for name, d in summary["drugs"].items()
             if d.get("gain") and d["n_heldout_in_range"] > 0
             and d["log10_rmse_model_heldout"] < d["log10_rmse_null_heldout"]]
    summary.update({
        "paired_vs_null_wins": wins, "paired_vs_null_losses": losses,
        "paired_vs_null_p": p_value, "discriminates_lines": discriminates,
        "held_out_pass": sum(r["pass"] for r in held),
        "held_out_null_pass": sum(r["null_pass"] for r in held),
        "held_out_n": len(held),
        "held_out_class_correct": sum(r["class_correct"] for r in held),
        "held_out_null_class_correct": sum(r["null_class_correct"] for r in held),
        "drugs_beating_null_heldout": beats,
        "exit_gate_met": discriminates,
        "wall_s": round(time.time() - t_start, 1)})
    say(f"\nheld-out sensitive/resistant calls: engine {summary['held_out_class_correct']}/"
        f"{len(held)}, constant-IC50 null {summary['held_out_null_class_correct']}/{len(held)}")
    say(f"held-out predictions inside the GDSC span: engine "
        f"{summary['held_out_pass']}/{len(held)}, constant-IC50 null "
        f"{summary['held_out_null_pass']}/{len(held)}   ({summary['wall_s']:.0f} s)")
    say(f"\nPaired per-line test vs the constant-IC50 null, held-out lines only:")
    say(f"  engine closer on {wins} lines, null closer on {losses}  ->  sign-test p = {p_value:.3f}")
    say(f"Phase-1 exit gate (engine beats the null on more held-out lines than it "
        f"loses on, p < 0.05): {'MET' if discriminates else 'NOT MET'}")
    if not discriminates:
        say("  i.e. the engine reproduces each drug's potency scale but does not yet "
            "tell cell lines apart better than a single constant.")
    say(f"  (per-drug RMSE comparison, which is not the gate: {beats or 'none'} "
        f"had the lower held-out RMSE)")
    return rows, summary


# ── main ──────────────────────────────────────────────────────────────
def _parse_param(text: str) -> tuple[str, float]:
    key, _, val = text.partition("=")
    names = {f.name for f in dataclasses.fields(Params)}
    if key not in names or not val:
        raise argparse.ArgumentTypeError(f"--param expects NAME=VALUE with NAME one of the "
                                         f"engine Params fields; got {text!r}")
    return key, float(val)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="16 cells per dose, 8 bisection steps")
    ap.add_argument("--cells", type=int, default=None, help="cells per dose (default 32)")
    ap.add_argument("--drugs", default=",".join(DRUGS), help="comma-separated subset")
    ap.add_argument("--param", type=_parse_param, action="append", default=[],
                    metavar="NAME=VALUE", help="override an engine parameter (repeatable)")
    ap.add_argument("--no-write", action="store_true", help="do not write the results files")
    a = ap.parse_args(argv)
    params = Params(**dict(a.param))
    rows, summary = run_validation(
        params=params, n_cells=a.cells or (16 if a.quick else 32),
        iters=8 if a.quick else 12,
        drug_names=tuple(d.strip() for d in a.drugs.split(",") if d.strip()))
    print("The span check is at least 9x wide, so it tests the potency scale, not "
          "line-to-line discrimination; the paired sign test above is the exit gate.")
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
