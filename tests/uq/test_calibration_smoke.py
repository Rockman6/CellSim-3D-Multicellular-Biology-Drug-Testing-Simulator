"""Layer 1.7 calibration smoke — CellSim vs published K_d.

Runs src.uq.calibration on the bundled streptavidin set (4
published binders with K_d spanning 10^-14 to 10^-5 M).

Gates:
  - docker succeeds on >= 3/4 entries
  - Spearman ρ >= 0.8 (CellSim ranks by affinity usably)
  - MAE reported (no absolute gate — Vina is known to
    systematically underestimate tight binders)
  - Conformal q95 finite and positive (calibration produced a
    real number for downstream interval use)

Run:
    conda activate cellsim
    python tests/uq/test_calibration_smoke.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cache import Cache  # noqa: E402
from src.uq import run_calibration  # noqa: E402


CAL_YAML = REPO_ROOT / "benchmarks" / "dock" / "streptavidin_calibration.yaml"


def test_calibration_streptavidin():
    assert CAL_YAML.exists(), f"missing {CAL_YAML}"
    with tempfile.TemporaryDirectory(prefix="cellsim-calib-") as tmp:
        cache = Cache(Path(tmp) / "c.sqlite")
        r = run_calibration(
            CAL_YAML, exhaustiveness=16, num_modes=3,
            seed=1, cpu=2, cache=cache)
        print(r.summary())

        assert r.n_ok >= 3, f"only {r.n_ok}/{r.n_points} docked"
        assert r.ok, f"run_calibration not ok: {r.reason}"
        assert r.spearman_rho is not None
        assert r.spearman_rho >= 0.8, (
            f"Spearman ρ = {r.spearman_rho:+.3f} < 0.8; "
            "CellSim ranks too poorly for triage on this set")
        assert r.mae_kcalmol is not None
        assert 0.0 < r.mae_kcalmol < 30.0
        assert r.conformal_q95_kcalmol is not None
        assert r.conformal_q95_kcalmol > 0.0
        cache.close()


if __name__ == "__main__":
    try:
        test_calibration_streptavidin()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
