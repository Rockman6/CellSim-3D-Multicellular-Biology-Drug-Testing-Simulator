"""Layer 1.3 end-to-end hydration ΔG smoke.

Verifies the full pipeline:
  SMILES → vacuum + solvated alchemical systems →
  Langevin+MBAR on both legs → ΔG_hyd

This is a PIPELINE smoke, not an accuracy gate. The default
smoke parameters (5 windows × 4 samples) give a wildly
undersampled number — the actual Milestone A accuracy gate
(FreeSolv MAE ≤ 1.5 kcal/mol) needs ≥50 ps per window on a
GPU runner and is a separate workflow_dispatch target.

This test asserts only that:
  - the call returns ok=True
  - ΔG_hyd is finite
  - vacuum leg ≈ 0 for methane (physics sanity)
  - ΔG_hyd SIGN is correct for methane (+, hydrophobic)
  - wall time stays under 3 min

The sign check (test_methane_hydration_sign_is_positive) is load-
bearing: the M5 Max pilot-1 run returned methane = -1.85 kcal/mol
because compute_hydration_dg had a spurious minus on the
composition formula. No test pinned the sign, so the Milestone A
tag shipped with that bug intact. Without a sign test, a future
regression could re-introduce the same flip silently.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


def test_methane_hydration_pipeline_runs():
    from src.fep import compute_hydration_dg
    r = compute_hydration_dg(
        "C", n_windows=5,
        n_production_steps=200, n_equilibration_steps=50,
        sample_stride=50, seed=1)
    print(f"  {r.summary()}")
    assert r.ok, r.reason
    assert r.dG_hydration_kcalmol is not None
    assert math.isfinite(r.dG_hydration_kcalmol)
    assert math.isfinite(r.uncertainty_kcalmol)
    assert r.wall_seconds is not None and r.wall_seconds < 180, (
        f"pipeline wall = {r.wall_seconds:.1f}s; "
        "expected < 3 min at smoke params.")
    # Vacuum decoupling for neutral methane must be ≈ 0 — no
    # λ-dependent self-interaction change. This is the load-
    # bearing physics sanity the sampling gate is also pinning.
    assert abs(r.dG_vacuum_decouple_kcalmol) < 1.0, (
        f"ΔG_vac = {r.dG_vacuum_decouple_kcalmol:+.2f} kcal/mol "
        "for methane should be ≈ 0; sampling pipeline broken.")


def test_methane_hydration_sign_is_positive():
    """Methane is hydrophobic — ΔG_hyd > 0 (FreeSolv expt: +2.00).

    Pins the sign of compute_hydration_dg's composition formula.
    This is the test that *would have caught* the Milestone A
    pilot-1 sign bug (M5 Max returned -1.85 for methane;
    milestone-a-pilot-1 tag shipped without ever checking sign).

    Runs slightly more sampling than the scaffold smoke so the
    sign is robust against noise: 7 windows × 500 prod steps ≈
    1-2 min on CPU. Tolerance is generous (>0.2 kcal/mol) because
    the true +2.00 converges only at production sampling; we just
    need the SIGN to be right.
    """
    from src.fep import compute_hydration_dg
    r = compute_hydration_dg(
        "C", n_windows=7,
        n_production_steps=500, n_equilibration_steps=100,
        sample_stride=50, seed=1)
    print(f"  {r.summary()}")
    assert r.ok, r.reason
    assert r.dG_hydration_kcalmol is not None
    assert math.isfinite(r.dG_hydration_kcalmol)
    # Gate: methane must be predicted hydrophobic (positive).
    # Margin 0.2 kcal/mol to avoid flaky near-zero signs.
    assert r.dG_hydration_kcalmol > 0.2, (
        f"ΔG_hyd(methane) = {r.dG_hydration_kcalmol:+.3f} "
        "kcal/mol; methane is hydrophobic, so ΔG_hyd must be "
        ">+0.2. A negative value indicates the Milestone A sign "
        "bug has returned (see BENCHMARKS.md § 'Milestone A "
        "post-mortem'). Fix: check the composition formula in "
        "compute_hydration_dg.")
    # Tightness gate (locked in after split-schedule fix b89dd51):
    # methane should predict within 1.0 kcal/mol of FreeSolv +2.00
    # at these smoke params (measured +2.05 at seed=1, residual
    # 0.05). If this fails, the schedule change has regressed and
    # water-penetration is back. Margin ±1 kcal/mol absorbs
    # seed-to-seed Langevin noise at 7×500 = 3 500 steps/leg.
    assert abs(r.dG_hydration_kcalmol - 2.00) < 1.0, (
        f"ΔG_hyd(methane) = {r.dG_hydration_kcalmol:+.3f} vs "
        "FreeSolv expt +2.00; residual > 1.0 kcal/mol at smoke "
        "params is outside the expected seed-noise band. The "
        "split-schedule fix (BENCHMARKS.md § 'root cause') may "
        "have regressed — check src/fep/sampling.py's "
        "_split_lambda_schedule and the (lam_e, lam_s) unpacks.")


if __name__ == "__main__":
    funcs = [
        test_methane_hydration_pipeline_runs,
        test_methane_hydration_sign_is_positive,
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
