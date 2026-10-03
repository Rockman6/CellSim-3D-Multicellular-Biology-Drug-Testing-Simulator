#!/usr/bin/env python3
"""Leave-one-mutant-out test of the p53-independent (p73) death route.

Question: does adding ATM/c-Abl -> p73 -> PUMA (Params.p73_gain) make the
engine predict TP53-mutant lines better *out of sample*?

For each candidate gain the full validation is re-run: drug constants
are refitted on the wild-type lines A549 + MCF7 (the route also acts in
wild-type cells), then every line is predicted. Then, for each fold:

    train mutant -> pick the gain with the lowest log10 error on it
    test mutant  -> score that gain there, against
                    (a) the p53-only engine (gain 0), and
                    (b) the constant-IC50 null.

The test mutant never influences the gain it is scored with. Only drugs
with an in-range GDSC reference on that line count (doxorubicin and
paclitaxel; cisplatin has only lower bounds on the mutant lines).

Usage:
    python scripts/experiment_p73.py            # ~15-20 min
    python scripts/experiment_p73.py --quick    # ~6 min, coarser
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params  # noqa: E402

_spec = importlib.util.spec_from_file_location("validate_gdsc", REPO_ROOT / "scripts" / "validate_gdsc.py")
vg = importlib.util.module_from_spec(_spec)
sys.modules["validate_gdsc"] = vg
_spec.loader.exec_module(vg)

GAINS = (0.0, 0.05, 0.1, 0.25, 0.5, 1.0)
OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "p73_loo_experiment.json"


def mutant_lines_with_data() -> tuple[str, ...]:
    """TP53-mutant lines that have at least one in-range GDSC reference,
    read from the reference extract so adding a line to the library is
    enough to widen the test."""
    from cellsim.cell.library import CELL_LINES
    targets = vg.load_targets(vg.REFERENCE_CSV)
    out = sorted({ln for (_, ln), t in targets.items()
                  if t.kind == "value" and not CELL_LINES[ln].p53_functional})
    return tuple(out)


def sq_errors(rows: list[dict], line: str, column: str) -> list[float]:
    out = []
    for r in rows:
        if r["cell_line"] == line and r["reference_kind"] == "value":
            pred = r[column]
            pred = 1e9 if pred == "inf" else float(pred)
            out.append(math.log10(pred / r["reference_ic50_uM"]) ** 2)
    return out


def rmse(errs: list[float]) -> float:
    return math.sqrt(sum(errs) / len(errs)) if errs else math.nan


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    n_cells, iters = (16, 8) if a.quick else (32, 12)
    global MUTANTS
    MUTANTS = mutant_lines_with_data()
    print(f"TP53-mutant lines with an in-range GDSC reference: {', '.join(MUTANTS)}")

    runs: dict[float, list[dict]] = {}
    print(f"{'p73_gain':>8}  " + "  ".join(f"{m:>11}" for m in MUTANTS)
          + "   (log10 RMSE per TP53-mutant line, over all drugs with an in-range reference)")
    for g in GAINS:
        rows, _ = vg.run_validation(params=Params(p73_gain=g), n_cells=n_cells, iters=iters,
                                    verbose=False)
        runs[g] = rows
        cells = [rmse(sq_errors(rows, ln, "predicted_ic50_uM")) for ln in MUTANTS]
        print(f"{g:>8}  " + "  ".join(f"{c:>11.3f}" for c in cells))

    null_rows = runs[0.0]
    print(f"{'null':>8}  " + "  ".join(
        f"{rmse(sq_errors(null_rows, ln, 'null_ic50_uM')):>11.3f}" for ln in MUTANTS))

    folds = []
    for test in MUTANTS:
        # Leave-one-out: choose the gain on every OTHER mutant line pooled,
        # then score it on the one left out. The test line never
        # contributes to the choice.
        train = [m for m in MUTANTS if m != test]
        def train_rmse(g, train=train):
            errs = []
            for ln in train:
                errs += sq_errors(runs[g], ln, "predicted_ic50_uM")
            return rmse(errs)
        best = min(GAINS, key=train_rmse)
        test_best = rmse(sq_errors(runs[best], test, "predicted_ic50_uM"))
        test_p53 = rmse(sq_errors(runs[0.0], test, "predicted_ic50_uM"))
        test_null = rmse(sq_errors(null_rows, test, "null_ic50_uM"))
        folds.append({"train": train, "test": test, "chosen_gain": best,
                      "test_rmse_with_route": test_best, "test_rmse_p53_only": test_p53,
                      "test_rmse_null": test_null,
                      "route_helps": test_best < test_p53, "beats_null": test_best < test_null})
        print(f"  held out {test:<12} chosen gain {best:<5} -> with route {test_best:.3f}, "
              f"p53-only {test_p53:.3f}, null {test_null:.3f}")

    helps = all(f["route_helps"] for f in folds)
    beats = all(f["beats_null"] for f in folds)
    verdict = ("route improves out-of-sample prediction in both folds"
               + (" and beats the null" if beats else ", but does not beat the null")
               if helps else "route does not improve out-of-sample prediction in both folds")
    print(f"\nverdict: {verdict}")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps({
            "gains": list(GAINS), "n_cells_per_dose": n_cells, "bisection_steps": iters,
            "mutant_lines": list(MUTANTS),
            "per_gain_rmse": {str(g): {ln: rmse(sq_errors(runs[g], ln, "predicted_ic50_uM"))
                                       for ln in MUTANTS} for g in GAINS},
            "null_rmse": {ln: rmse(sq_errors(null_rows, ln, "null_ic50_uM"))
                          for ln in MUTANTS},
            "folds": folds, "verdict": verdict}, indent=2) + "\n")
        print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
