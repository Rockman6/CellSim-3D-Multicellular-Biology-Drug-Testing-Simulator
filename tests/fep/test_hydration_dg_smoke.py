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
  - wall time stays under 3 min
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


if __name__ == "__main__":
    try:
        test_methane_hydration_pipeline_runs()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
