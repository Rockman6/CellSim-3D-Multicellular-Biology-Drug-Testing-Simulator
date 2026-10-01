"""Intracellular binding sink: free vs total drug, and the buffering it does.

Pins the free/total round-trip (invertibility), the no-sink and
strong-sink limits, saturation, and that the sink lowers the free fraction
the target actually sees.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    BindingSink,
    total_from_free,
    free_from_total,
)


def test_no_sink_means_free_equals_total():
    sink = BindingSink(capacity_M=0.0, Kd_M=1e-6)
    r = free_from_total(1e-7, sink)
    assert abs(r.free_M - r.total_M) < 1e-20
    assert abs(r.free_fraction - 1.0) < 1e-12


def test_free_total_round_trip_is_consistent():
    """total_from_free and free_from_total must invert each other."""
    sink = BindingSink(capacity_M=1e-4, Kd_M=1e-6)
    free = 3e-7
    total = total_from_free(free, sink).total_M
    back = free_from_total(total, sink).free_M
    assert abs(back - free) / free < 1e-9


def test_strong_sink_buffers_free_far_below_total():
    """A large, tight sink holds most of the drug bound (fu ≪ 1)."""
    sink = BindingSink(capacity_M=1e-3, Kd_M=1e-8)   # abundant, tight
    r = free_from_total(1e-6, sink)
    assert r.free_M < r.total_M
    assert r.free_fraction < 0.01


def test_sink_saturates_and_free_fraction_recovers():
    """Once free ≫ Ks the sink is full, so fu climbs back toward 1."""
    sink = BindingSink(capacity_M=1e-7, Kd_M=1e-8)
    low = total_from_free(1e-10, sink)     # free ≪ Ks
    high = total_from_free(1e-5, sink)     # free ≫ Ks, sink swamped
    assert high.free_fraction > low.free_fraction
    assert high.sink_saturated_fraction > 0.9
    assert low.sink_saturated_fraction < high.sink_saturated_fraction


def test_bound_plus_free_equals_total():
    sink = BindingSink(capacity_M=5e-5, Kd_M=2e-7)
    r = free_from_total(1e-6, sink)
    assert abs((r.free_M + r.bound_M) - r.total_M) < 1e-18


def test_rejects_bad_inputs():
    try:
        BindingSink(capacity_M=1e-4, Kd_M=0.0)
        assert False
    except ValueError:
        pass
    try:
        free_from_total(-1.0, BindingSink(capacity_M=1e-4, Kd_M=1e-6))
        assert False
    except ValueError:
        pass
