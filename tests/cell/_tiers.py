"""Two tiers of cell test, and which one runs where.

A MECHANISM gate checks that the code does what it claims: a knockout of
a gene the engine does not model is refused, effectors with no kill rate
kill nothing, an untouched dish is bit-identical to before. These are
fast and catch most regressions, so they run on every pull request.

A VALIDATION gate checks the engine against a measurement: several
six-day spheroids at different glucose levels, a calibration recovered
from a simulated plate. They take minutes each and rarely change between
pull requests, so they run nightly. A regression they would catch is
still caught within a day.

Nothing is weakened by this. Every gate still runs; the split changes
WHEN the slow ones run, not WHAT they check. The tempting alternative —
shrinking a slow test until it is fast — is how two tests in this suite
once passed while measuring nothing, because they had been made too
small for the effect to appear.

Running a test file directly runs everything. CI's per-PR job sets
CELLSIM_TESTS=mechanism to skip the validation tier.
"""
from __future__ import annotations

import os

TIER = os.environ.get("CELLSIM_TESTS", "all").strip().lower()


def validation(fn):
    """Mark a test as a slow check against measurement (nightly tier)."""
    fn.tier = "validation"
    return fn


def skipped(fn) -> bool:
    """True when this run should leave `fn` for the nightly job."""
    return TIER == "mechanism" and getattr(fn, "tier", "mechanism") == "validation"


def run(namespace: dict, label: str) -> int:
    """The shared runner: every test_ function in `namespace`, honouring
    the tier, reporting skips explicitly so a quick run never looks like a
    full one."""
    fns = [v for k, v in sorted(namespace.items()) if k.startswith("test_")]
    failed = skipped_n = 0
    for fn in fns:
        if skipped(fn):
            skipped_n += 1
            print(f"  skip {fn.__name__}  (validation tier: nightly)")
            continue
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
        except Exception as e:                      # noqa: BLE001 - report, don't crash
            failed += 1
            print(f"  ERROR {fn.__name__}: {type(e).__name__}: {e}")
    ran = len(fns) - skipped_n
    note = f", {skipped_n} left for the nightly run" if skipped_n else ""
    print(f"{ran - failed}/{ran} {label} gates pass{note}")
    return 1 if failed else 0
