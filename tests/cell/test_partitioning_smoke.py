"""pH partitioning / ion trapping: why the concentration inside differs.

Pins the Henderson-Hasselbalch neutral fraction, the accumulation
direction (weak base → acidic compartment, weak acid → basic), the
neutral-drug null case, and a lysosomotropic-base example with the
right order of magnitude.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    COMPARTMENT_PH,
    neutral_fraction,
    accumulation_ratio,
    partition_across_membrane,
)


def test_half_ionised_at_pH_equals_pKa():
    assert abs(neutral_fraction(7.0, 7.0, "base") - 0.5) < 1e-12
    assert abs(neutral_fraction(5.0, 5.0, "acid") - 0.5) < 1e-12


def test_neutral_drug_never_partitions():
    assert neutral_fraction(4.0, 8.0, "neutral") == 1.0
    r = accumulation_ratio(pKa=8.0, pH_in=4.7, pH_out=7.4, ion_type="neutral")
    assert abs(r - 1.0) < 1e-12


def test_weak_base_accumulates_in_acidic_lysosome():
    """A weak base concentrates where it is more protonated (low pH)."""
    R = accumulation_ratio(pKa=9.0, pH_in=COMPARTMENT_PH["lysosome"],
                           pH_out=COMPARTMENT_PH["blood"], ion_type="base")
    assert R > 1.0
    # pKa 9, blood 7.4 → lysosome 4.7: base is ~fully protonated both sides
    # but far more trapped in the lysosome; ratio is large (>100×).
    assert R > 100.0


def test_weak_acid_is_excluded_from_acidic_compartment():
    """A weak acid does the opposite — accumulates in BASIC compartments."""
    R_acidic = accumulation_ratio(pKa=4.0, pH_in=COMPARTMENT_PH["lysosome"],
                                  pH_out=COMPARTMENT_PH["blood"], ion_type="acid")
    assert R_acidic < 1.0
    R_basic = accumulation_ratio(pKa=4.0, pH_in=COMPARTMENT_PH["mitochondrion"],
                                 pH_out=COMPARTMENT_PH["blood"], ion_type="acid")
    assert R_basic > 1.0


def test_partition_scales_concentration():
    res = partition_across_membrane(
        1e-7, pKa=9.0, pH_in=COMPARTMENT_PH["lysosome"], ion_type="base")
    assert abs(res.total_in_M - res.accumulation_ratio * 1e-7) < 1e-20
    assert res.total_in_M > res.total_out_M
    assert "accumulates INSIDE" in res.summary()


def test_direction_flips_with_pH_gradient_sign():
    # Same base, but target compartment more BASIC than blood → excluded.
    R = accumulation_ratio(pKa=9.0, pH_in=COMPARTMENT_PH["mitochondrion"],
                           pH_out=COMPARTMENT_PH["blood"], ion_type="base")
    assert R < 1.0


def test_rejects_bad_inputs():
    try:
        neutral_fraction(7.0, 7.0, "zwitterion")
        assert False
    except ValueError:
        pass
    try:
        partition_across_membrane(-1.0, pKa=9.0, pH_in=4.7)
        assert False
    except ValueError:
        pass
