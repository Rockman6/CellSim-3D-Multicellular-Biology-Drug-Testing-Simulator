"""Layer 1.6 conformal smoke — split-conformal CI on synthetic data.

Verifies src.uq.conformal.ConformalBounds:
  - calibrate runs on N=50 synthetic (pred, truth) pairs
  - returned quantile q is positive and finite
  - interval(pred) brackets pred with half-width q
  - empirical coverage on N=200 held-out test set is >= 0.90
    (theoretical target 0.95; allow some slack for finite-N jitter)

Run:
    conda activate cellsim
    python tests/uq/test_conformal_smoke.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.uq import ConformalBounds  # noqa: E402


def test_conformal_synthetic():
    import numpy as np
    rng = np.random.default_rng(42)
    # Calibration set (Vina-like noise around truth).
    n_cal = 50
    true_cal = rng.uniform(-10.0, -4.0, size=n_cal)
    pred_cal = true_cal + rng.normal(0.0, 0.8, size=n_cal)

    cb = ConformalBounds(alpha=0.05)
    info = cb.calibrate(pred_cal.tolist(), true_cal.tolist())
    print(f"[conformal] calibrate: {info}")

    assert cb.q is not None and cb.q > 0.0 and cb.q < 10.0, (
        f"bad q = {cb.q}")
    assert info["n_cal"] == n_cal
    assert info["coverage_target"] == 0.95

    # Interval centred on prediction with half-width q.
    lo, hi = cb.interval(-7.0)
    half = (hi - lo) / 2.0
    assert abs(half - cb.q) < 1e-9, f"half={half} vs q={cb.q}"
    assert lo == -7.0 - cb.q
    assert hi == -7.0 + cb.q

    # Empirical coverage on a fresh test set.
    n_test = 200
    true_test = rng.uniform(-10.0, -4.0, size=n_test)
    pred_test = true_test + rng.normal(0.0, 0.8, size=n_test)
    cov = cb.check_coverage(pred_test.tolist(), true_test.tolist())
    print(f"[conformal] n_test={n_test}  coverage={cov:.3f}  "
          f"target=0.95")
    # Conformal guarantees >= (1-alpha) in expectation; with finite
    # N the observed coverage fluctuates. Accept >= 0.90 as the
    # smoke gate.
    assert cov >= 0.90, f"coverage {cov:.3f} below 0.90 gate"


if __name__ == "__main__":
    try:
        test_conformal_synthetic()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
