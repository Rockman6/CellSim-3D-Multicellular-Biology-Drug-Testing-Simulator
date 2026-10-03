#!/usr/bin/env python3
"""The Phase-1 validation protocol's reference logic (no simulation).

`scripts/validate_gdsc.py` turns GDSC screens into acceptance spans. A
bug there would silently grade every prediction wrong, so the rules are
pinned here against the committed reference extract.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

_spec = importlib.util.spec_from_file_location(
    "validate_gdsc", REPO_ROOT / "scripts" / "validate_gdsc.py")
vg = importlib.util.module_from_spec(_spec)
sys.modules["validate_gdsc"] = vg   # dataclasses resolves annotations via sys.modules
_spec.loader.exec_module(vg)

TARGETS = vg.load_targets(vg.REFERENCE_CSV)


def test_reference_covers_the_three_drugs():
    drugs = {d for (d, _) in TARGETS}
    assert drugs == {"cisplatin", "doxorubicin", "paclitaxel"}, drugs


def test_span_is_at_least_three_fold_either_side():
    for key, t in TARGETS.items():
        if t.kind == "value":
            assert t.lo <= t.center / vg.MIN_SPAN * (1 + 1e-9), key
            assert t.hi >= t.center * vg.MIN_SPAN * (1 - 1e-9), key


def test_span_contains_every_in_range_replicate():
    """The doxorubicin replicates disagree by up to 9.3x; the span must
    hold both, so no prediction is graded tighter than the data."""
    t = TARGETS[("doxorubicin", "MCF7")]          # 0.010 and 0.093 µM screens
    assert t.kind == "value" and t.n_screens == 2
    assert t.lo <= 0.0101 and t.hi >= 0.093
    assert math.isclose(t.center, math.sqrt(0.01003 * 0.09307), rel_tol=1e-2)


def test_screens_that_never_reached_half_kill_give_only_a_lower_bound():
    """Cisplatin on MCF7: all three screens' IC50s lie above their top dose."""
    t = TARGETS[("cisplatin", "MCF7")]
    assert t.kind == "lower_bound"
    assert t.center == 10.0 and math.isinf(t.hi)
    assert math.isclose(t.lo, 10.0 / vg.MIN_SPAN)


def test_out_of_range_screens_are_ignored_when_an_in_range_one_exists():
    """Doxorubicin on MDA-MB-231: 1.26 µM sits above its 1.024 µM top and
    must not pull the reference; 0.108 µM is the in-range screen."""
    t = TARGETS[("doxorubicin", "MDA-MB-231")]
    assert t.kind == "value" and math.isclose(t.center, 0.1081, rel_tol=1e-2)


def test_fit_lines_are_clean_p53_wild_type():
    from cellsim.cell.library import CELL_LINES
    for ln in vg.FIT_LINES:
        assert CELL_LINES[ln].p53_functional, ln
    assert "HeLa" not in vg.FIT_LINES, "HeLa p53 is degraded by HPV E6; keep it held out"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} target-derivation gates pass")
    raise SystemExit(1 if failed else 0)
