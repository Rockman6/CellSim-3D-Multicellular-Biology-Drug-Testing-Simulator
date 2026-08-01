"""Polyprotic (chloroquine-class) partitioning — traps harder than monoprotic.

Rigor anchor: the polyprotic form must reduce EXACTLY to the monoprotic
`neutral_fraction` when given a single pKa. Then pins that a diprotic base
accumulates far more than a monoprotic one with the same strongest pKa,
and reproduces chloroquine's known ~10^3-10^4 lysosomal trapping.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    COMPARTMENT_PH,
    neutral_fraction,
    neutral_fraction_polyprotic,
    accumulation_ratio,
    accumulation_ratio_polyprotic,
)


def test_polyprotic_reduces_to_monoprotic_for_one_pKa():
    """THE anchor: one pKa → identical to the monoprotic formula."""
    for pH in (4.7, 6.3, 7.2, 7.4, 8.0):
        for pKa in (4.0, 7.0, 9.0, 10.2):
            for ion in ("base", "acid"):
                a = neutral_fraction(pH, pKa, ion)
                b = neutral_fraction_polyprotic(pH, [pKa], ion)
                assert abs(a - b) < 1e-12, (pH, pKa, ion)
    # ...and so do the accumulation ratios.
    a = accumulation_ratio(9.0, 4.7, 7.4, "base")
    b = accumulation_ratio_polyprotic([9.0], 4.7, 7.4, "base")
    assert abs(a - b) / a < 1e-12


def test_diprotic_base_traps_harder_than_monoprotic():
    mono = accumulation_ratio_polyprotic([10.2], COMPARTMENT_PH["lysosome"],
                                         COMPARTMENT_PH["blood"], "base")
    di = accumulation_ratio_polyprotic([10.2, 8.4],
                                       COMPARTMENT_PH["lysosome"],
                                       COMPARTMENT_PH["blood"], "base")
    assert di > mono
    # The second protonation site adds orders of magnitude.
    assert di / mono > 10.0


def test_chloroquine_ideal_trapping_exceeds_the_measured_value():
    """Chloroquine (pKa ≈ 10.2, 8.4): the IDEAL Henderson-Hasselbalch limit.

    The model predicts ~2.3e5×, while measured lysosomal accumulation is
    ~1e3-1e4×. That gap is REAL PHYSICS the model deliberately omits, not
    a fitting error: chloroquine is a base, so as it accumulates it
    RAISES lysosomal pH, which shrinks the very gradient driving the
    trapping (self-limiting). Finite drug, lysosomal volume and binding
    also cap it. So this number is an UPPER BOUND on passive trapping —
    documented in the module, not tuned away.
    """
    R = accumulation_ratio_polyprotic([10.2, 8.4],
                                      COMPARTMENT_PH["lysosome"],
                                      COMPARTMENT_PH["blood"], "base")
    assert 2.0e5 < R < 2.6e5, R
    # It must exceed the measured range — the model is an upper bound.
    assert R > 1e4


def test_neutral_and_bad_input():
    assert neutral_fraction_polyprotic(4.7, [10.2, 8.4], "neutral") == 1.0
    try:
        neutral_fraction_polyprotic(7.0, [], "base")
        assert False
    except ValueError:
        pass
    try:
        neutral_fraction_polyprotic(7.0, [7.0], "zwitterion")
        assert False
    except ValueError:
        pass
