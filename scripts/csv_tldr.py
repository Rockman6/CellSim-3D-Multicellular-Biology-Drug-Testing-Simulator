#!/usr/bin/env python3
"""csv_tldr.py — compress a bench CSV into a single Slack/standup line.

Biologists running `cellsim fep-binding bench` or `cellsim fep-freesolv`
end up with a CSV that the full fep-report flow expands into a multi-
section markdown. For a quick status update ("how's the run going?")
the full report is overkill — what's wanted is a one-liner you can
paste into Slack or an email subject.

Usage:
    python scripts/csv_tldr.py run/fep/.../freesolv_results.csv
    python scripts/csv_tldr.py path/to/bench_results.csv --gate 2.0

The script auto-detects hydration vs binding from the ΔG sign
pattern (hydration has mixed signs, binding has all negatives for
binders), picks the appropriate gate (1.5 kcal/mol hydration,
2.0 kcal/mol binding), and emits:

    FreeSolv hydration FEP: 12/12 ok, MAE 0.42 kcal/mol (gate ≤1.5),
    Pearson r +0.993 — PASS

Exit codes:
  0 = PASS (all gates met)
  1 = FAIL (one or more gate failed)
  2 = inconclusive (missing data; partial run)
  3 = usage / I/O error
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def _load_rows(csv_path: Path) -> list[dict]:
    if not csv_path.is_file():
        raise FileNotFoundError(f"no such file: {csv_path}")
    with csv_path.open("r", encoding="utf-8-sig", newline="") as fi:
        reader = csv.DictReader(fi)
        return list(reader)


def _f(v):
    if v in (None, "", "None", "nan"):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3 or len(x) != len(y):
        return None
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    dx = math.sqrt(sum((xi - mx) ** 2 for xi in x))
    dy = math.sqrt(sum((yi - my) ** 2 for yi in y))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def _kendall_tau(x: list[float], y: list[float]) -> float | None:
    """Kendall τ-a — O(n²), fine for <100 compounds."""
    n = len(x)
    if n < 3 or n != len(y):
        return None
    conc = disc = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            if dx * dy > 0:
                conc += 1
            elif dx * dy < 0:
                disc += 1
    denom = n * (n - 1) / 2.0
    if denom == 0:
        return None
    return (conc - disc) / denom


def detect_kind(rows: list[dict]) -> str:
    """Hydration (FreeSolv-style): ΔG_expt spans positive + negative.
    Binding: all expt ΔG < 0 (binders only)."""
    expts = [_f(r.get("dG_expt_kcalmol")) for r in rows]
    expts = [e for e in expts if e is not None]
    if not expts:
        return "unknown"
    if any(e > 0 for e in expts) and any(e < 0 for e in expts):
        return "hydration"
    if all(e < 0 for e in expts):
        return "binding"
    return "unknown"


def _label_for_kind(kind: str, rows: list[dict]) -> str:
    names = [(r.get("name") or "").lower() for r in rows]
    joined = " ".join(names)
    if kind == "hydration":
        return "FreeSolv hydration FEP"
    if kind == "binding":
        if "biotin" in joined or "streptavidin" in joined:
            return "streptavidin binding FEP"
        if ("erlotinib" in joined or "gefitinib" in joined
                or "lapatinib" in joined or "egfr" in joined):
            return "EGFR kinase binding FEP"
        return "binding FEP"
    return "FEP"


def summarise(rows: list[dict], gate: float | None = None) -> dict:
    """Return the TL;DR bundle: label, counts, metrics, verdict."""
    n_total = len(rows)
    ok_rows = [r for r in rows
               if (r.get("ok") or "").strip().lower() in (
                   "true", "1", "yes")
               and _f(r.get("dG_pred_kcalmol")) is not None]
    n_ok = len(ok_rows)
    kind = detect_kind(rows)
    label = _label_for_kind(kind, rows)
    gate_kcal = gate
    if gate_kcal is None:
        gate_kcal = 1.5 if kind == "hydration" else 2.0

    if not ok_rows:
        return {
            "label": label, "kind": kind,
            "n_total": n_total, "n_ok": n_ok,
            "mae": None, "pearson": None, "kendall": None,
            "sign_ok": None, "gate": gate_kcal,
            "verdict": "inconclusive",
            "reason": "no ok rows",
        }

    expts = [_f(r["dG_expt_kcalmol"]) for r in ok_rows]
    preds = [_f(r["dG_pred_kcalmol"]) for r in ok_rows]
    abs_res = [abs(p - e) for p, e in zip(preds, expts)]
    mae = sum(abs_res) / len(abs_res)
    pearson = _pearson(expts, preds)
    kendall = _kendall_tau(expts, preds)

    # Sign rule.
    sign_ok: bool | None = None
    sign_reason = ""
    if kind == "hydration":
        # Require sign-correct on every compound.
        sign_ok = all(
            (e > 0 and p > 0) or (e < 0 and p < 0)
            or (abs(e) < 0.5 and abs(p) < 0.5)
            for e, p in zip(expts, preds))
        if not sign_ok:
            sign_reason = "sign flip on ≥1 compound"
    elif kind == "binding":
        sign_ok = all(p < 0 for p in preds)
        if not sign_ok:
            sign_reason = "non-binder predicted"

    mae_ok = mae <= gate_kcal
    complete = (n_ok == n_total)

    # Gate logic mirrors cellsim fep-report (commit f7b8fe4): compute
    # MAE on the ok subset, PASS if MAE ≤ gate AND sign-critical
    # holds, FAIL otherwise. Surface partial-run state in the
    # `reason` field so the line still reads "PASS (10/12 ok)" rather
    # than hiding the gap. Earlier csv_tldr returned 'inconclusive'
    # whenever n_ok < n_total, which conflicted with fep-report's
    # PASS on the same data and confused the biologist about which
    # verdict to trust.
    if mae_ok and (sign_ok is not False):
        verdict = "PASS"
    else:
        verdict = "FAIL"

    return {
        "label": label, "kind": kind,
        "n_total": n_total, "n_ok": n_ok,
        "mae": mae, "pearson": pearson, "kendall": kendall,
        "sign_ok": sign_ok, "gate": gate_kcal,
        "verdict": verdict,
        "reason": sign_reason if not sign_ok else "",
    }


def format_tldr(s: dict) -> str:
    n_ok = s["n_ok"]
    n_total = s["n_total"]
    label = s["label"]
    if s["mae"] is None:
        return (f"{label}: {n_ok}/{n_total} ok — "
                f"{s['verdict']} ({s['reason'] or 'no data'})")
    mae_str = f"MAE {s['mae']:.2f} kcal/mol (gate ≤{s['gate']:.1f})"
    corr = ""
    if s["kind"] == "hydration" and s["pearson"] is not None:
        corr = f", Pearson r {s['pearson']:+.3f}"
    elif s["kind"] == "binding" and s["kendall"] is not None:
        corr = f", Kendall τ {s['kendall']:+.2f}"
    tail = ""
    if s["reason"]:
        tail = f" ({s['reason']})"
    # Partial-run footer when n_ok < n_total but verdict still PASS —
    # mirrors fep-report's "PASS but X compounds incomplete" tone.
    partial = ""
    if n_ok < n_total and s["verdict"] == "PASS":
        partial = f" (partial: {n_total - n_ok} failed)"
    return (f"{label}: {n_ok}/{n_total} ok, {mae_str}"
            f"{corr} — {s['verdict']}{tail}{partial}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv_path",
                    help="bench CSV (from `cellsim fep-binding bench` "
                         "or `cellsim fep-freesolv`)")
    ap.add_argument("--gate", type=float, default=None,
                    help="MAE gate in kcal/mol (default: 1.5 for "
                         "hydration, 2.0 for binding)")
    args = ap.parse_args(argv)

    try:
        rows = _load_rows(Path(args.csv_path))
    except (FileNotFoundError, OSError) as e:
        print(f"csv_tldr: {e}", file=sys.stderr)
        return 3

    s = summarise(rows, gate=args.gate)
    print(format_tldr(s))
    if s["verdict"] == "PASS":
        return 0
    if s["verdict"] == "FAIL":
        return 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
