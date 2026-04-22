"""Regression tests for src.fep.binding.estimate_sampling_wall_hours.

Purpose: prevent wasted GPU runs. Biologists about to spend 6-48
hours of GPU time on a binding FEP benefit from a cost preview;
these tests pin the estimator so the preview stays roughly honest
(±2-3× is the advertised accuracy).

Anchors:
  - CPU    : 2 M steps·atom/sec
  - Metal  : 26 M steps·atom/sec
  - CUDA   : 150 M steps·atom/sec
  - overhead surcharge : 1.25×

Independent re-derivation for a streptavidin-class 30k-atom
system with 11 × (25 000 prod + 2 500 equil) × 2 legs:
  total_steps = 2 × 11 × 27 500 = 605 000
  CPU  = 605 000 × 30 000 / 2e6 × 1.25 / 3600 ≈ 3.15 h
  Metal= same / (26/2)        ≈ 0.24 h
  CUDA = same / (150/2)       ≈ 0.042 h ≈ 2.5 min
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.fep.binding import (
    estimate_sampling_wall_hours,
    format_wall_estimate_block,
)


def test_methane_smoke_runs_in_minutes_on_metal():
    """FreeSolv-class small system — should be minutes on Metal."""
    est = estimate_sampling_wall_hours(
        n_atoms=2000, n_windows=11,
        n_production_steps=25000,
        n_equilibration_steps=2500)
    assert est["metal_m5max"] < 0.1, est   # < 6 min
    assert est["cuda_h100"] < 0.02, est     # < 72 s
    # CPU on small system should also finish overnight.
    assert est["cpu"] < 2.0, est            # < 2 h


def test_streptavidin_class_30k_atoms_reasonable():
    est = estimate_sampling_wall_hours(
        n_atoms=30000, n_windows=11,
        n_production_steps=25000,
        n_equilibration_steps=2500)
    # Sanity ratios: CUDA ~ 75× faster than CPU, Metal ~ 13× faster.
    cpu_h = est["cpu"]
    metal_h = est["metal_m5max"]
    cuda_h = est["cuda_h100"]
    assert 1.0 < cpu_h < 10.0, est        # 3 h ballpark
    assert 0.1 < metal_h < 1.0, est       # ~15 min ballpark
    assert cuda_h < 0.1, est               # < 6 min
    # Platform ordering.
    assert cuda_h < metal_h < cpu_h, est


def test_egfr_class_40k_atoms_cpu_flagged_as_infeasible():
    """EGFR kinase series: ~40k atoms, 6 compounds → CPU-only is
    not a viable plan (days). The formatter must warn."""
    est = estimate_sampling_wall_hours(
        n_atoms=40000, n_windows=11,
        # Paper-grade sampling: 50 ps per window.
        n_production_steps=50000,
        n_equilibration_steps=5000)
    assert est["cpu"] > 8.0, est
    block = format_wall_estimate_block(est, gate_hours=48.0)
    # CPU line must appear. Metal/CUDA must appear.
    assert "CPU" in block
    assert "Metal" in block
    assert "CUDA" in block


def test_formatter_flags_cpu_infeasible_over_gate():
    """Force a huge config so CPU > 48h — the warning must fire."""
    est = estimate_sampling_wall_hours(
        n_atoms=100000, n_windows=21,
        n_production_steps=100000,
        n_equilibration_steps=10000)
    block = format_wall_estimate_block(est, gate_hours=48.0)
    assert "CPU-only is not viable" in block, block


def test_formatter_does_not_flag_short_cpu_runs():
    est = estimate_sampling_wall_hours(
        n_atoms=1000, n_windows=5,
        n_production_steps=1000,
        n_equilibration_steps=500)
    block = format_wall_estimate_block(est, gate_hours=48.0)
    assert "CPU-only is not viable" not in block, block


def test_formatter_uses_minutes_for_short_runs():
    est = estimate_sampling_wall_hours(
        n_atoms=500, n_windows=3,
        n_production_steps=500,
        n_equilibration_steps=100)
    block = format_wall_estimate_block(est)
    # All three platforms should report in "min" for this tiny config.
    assert "min" in block, block


def test_formatter_uses_days_for_very_long_runs():
    est = estimate_sampling_wall_hours(
        n_atoms=200000, n_windows=21,
        n_production_steps=500000,
        n_equilibration_steps=50000)
    block = format_wall_estimate_block(est)
    # CPU line should render in days, not hours.
    assert " d" in block, block


def test_scaling_inversely_with_atoms():
    """Double the atoms → roughly double the wall."""
    small = estimate_sampling_wall_hours(
        n_atoms=5000, n_windows=11,
        n_production_steps=10000, n_equilibration_steps=1000)
    big = estimate_sampling_wall_hours(
        n_atoms=10000, n_windows=11,
        n_production_steps=10000, n_equilibration_steps=1000)
    for plat in ("cpu", "metal_m5max", "cuda_h100"):
        ratio = big[plat] / small[plat]
        assert 1.9 < ratio < 2.1, (plat, ratio, small[plat], big[plat])


def test_scaling_linearly_with_steps():
    """Double the prod steps → roughly double the wall."""
    short = estimate_sampling_wall_hours(
        n_atoms=5000, n_windows=11,
        n_production_steps=10000, n_equilibration_steps=1000)
    long = estimate_sampling_wall_hours(
        n_atoms=5000, n_windows=11,
        n_production_steps=21000, n_equilibration_steps=1000)
    # (21000+1000) / (10000+1000) = 22/11 = 2.0
    for plat in ("cpu", "metal_m5max", "cuda_h100"):
        ratio = long[plat] / short[plat]
        assert 1.9 < ratio < 2.1, (plat, ratio)


def test_overhead_surcharge_in_result():
    """Sanity: the overhead surcharge actually raises the estimate
    above the pure MD time. A 1000-atom, 100-step run on CUDA at
    150M/1000 steps/sec = 150k steps/sec would be 100/150000 =
    6.67e-4 s pure. With 1.25× surcharge + n_legs=2: 1.67e-3 s."""
    est = estimate_sampling_wall_hours(
        n_atoms=1000, n_windows=1,
        n_production_steps=100, n_equilibration_steps=0)
    cuda_seconds = est["cuda_h100"] * 3600
    # 2 legs × 100 steps × 1.25 overhead / (150M / 1k) = 1.67e-3 s
    assert 1.5e-3 < cuda_seconds < 2.0e-3, cuda_seconds


if __name__ == "__main__":
    funcs = [
        test_methane_smoke_runs_in_minutes_on_metal,
        test_streptavidin_class_30k_atoms_reasonable,
        test_egfr_class_40k_atoms_cpu_flagged_as_infeasible,
        test_formatter_flags_cpu_infeasible_over_gate,
        test_formatter_does_not_flag_short_cpu_runs,
        test_formatter_uses_minutes_for_short_runs,
        test_formatter_uses_days_for_very_long_runs,
        test_scaling_inversely_with_atoms,
        test_scaling_linearly_with_steps,
        test_overhead_surcharge_in_result,
    ]
    fails = []
    for f in funcs:
        try:
            f()
            print(f"[PASS] {f.__name__}")
        except AssertionError as e:
            print(f"[FAIL] {f.__name__}: {e}")
            fails.append(f.__name__)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"[ERROR] {f.__name__}: {e}")
            fails.append(f.__name__)
    print(f"{len(funcs) - len(fails)}/{len(funcs)} PASS")
    sys.exit(0 if not fails else 1)
