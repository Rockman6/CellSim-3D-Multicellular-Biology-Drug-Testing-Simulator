#!/usr/bin/env python3
"""The library must carry the constants the validation fitted.

`scripts/validate_gdsc.py` fits one potency constant per drug and writes
it to `benchmarks/cell/gdsc_validation_summary.json`; the published
scores are computed with those fitted values. Everything a user runs —
the Lab, the CLI, `cellsim.api` — reads `cellsim/cell/library.py`
instead. Nothing tied the two together, and in October 2026 seven of
the ten drugs in the library still carried their provisional
placeholders: nutlin-3a 77x too potent (an A549 IC50 of 0.06 uM against
GDSC's 4-7 uM), SN-38 and gemcitabine over a hundred-fold too weak. The
validation page graded the right numbers while the product ran the wrong
ones.

This fails the moment the two disagree by more than 2 %.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _tiers import run  # noqa: E402
from cellsim.cell.library import DRUGS  # noqa: E402

SUMMARY = REPO_ROOT / "benchmarks" / "cell" / "gdsc_validation_summary.json"
TOLERANCE = 0.02


def test_every_fitted_constant_is_in_the_library():
    fitted = json.loads(SUMMARY.read_text())["drugs"]
    wrong = []
    for name, fit in fitted.items():
        lib = getattr(DRUGS[name], fit["fit_target"])
        if abs(lib / fit["gain"] - 1.0) > TOLERANCE:
            wrong.append(f"{name}.{fit['fit_target']}: library {lib:.4g}, fitted {fit['gain']:.4g}")
    assert not wrong, "library constants differ from the GDSC fit:\n  " + "\n  ".join(wrong)


def test_every_library_drug_with_a_fit_target_was_fitted():
    fitted = json.loads(SUMMARY.read_text())["drugs"]
    missing = [n for n, d in DRUGS.items() if d.fit_target and n not in fitted]
    assert not missing, f"drugs with a fit_target but no fit in the summary: {missing}"


def test_the_fit_used_the_engine_parameters_now_in_force():
    """A fit made under different engine parameters is a fit of a
    different model."""
    import dataclasses
    from cellsim.cell.engine import Params
    used = json.loads(SUMMARY.read_text())["params"]
    now = dataclasses.asdict(Params())
    differ = {k: (used.get(k), v) for k, v in now.items() if used.get(k) != v}
    assert not differ, f"Params changed since the GDSC fit (fit, now): {differ}"


if __name__ == "__main__":
    raise SystemExit(run(globals(), "library constants"))
