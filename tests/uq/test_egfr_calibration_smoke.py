"""Layer 1.7 EGFR kinase calibration smoke.

Pins the honest Vina-on-kinases finding: Vina under-scores
ATP-hinge H-bonded quinazolines and ranks tight binders worse
than weak parent 4-anilinoquinazoline.

Gates (all must hold — these are regression witnesses, not
aspirational targets):
  - calibration harness runs ok on 6/6 compounds
  - Spearman ρ is NEGATIVE (rank-order inverted)
  - MAE is ≥ 1.5 kcal/mol (under-scoring of tight binders)

If any of these gates *improves* (ρ becomes positive, MAE drops),
that is a signal that either (a) Vina/Meeko upstream silently
improved, (b) the receptor prep changed, or (c) the YAML changed
— in any case, BENCHMARKS.md should be re-run and updated.

Run:
    conda activate cellsim
    python tests/uq/test_egfr_calibration_smoke.py
"""

from __future__ import annotations

import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cache import Cache  # noqa: E402
from src.uq.calibration import run_calibration  # noqa: E402


YAML = REPO_ROOT / "benchmarks" / "dock" / "egfr_calibration.yaml"


def test_egfr_calibration_honest_failure():
    with tempfile.TemporaryDirectory(prefix="cellsim-egfr-") as tmp:
        cache = Cache(Path(tmp) / "c.sqlite")
        t0 = time.time()
        # Low exhaustiveness to keep CI wall-time manageable; the
        # finding is robust to exh=32 (verified manually).
        r = run_calibration(
            YAML, exhaustiveness=8, num_modes=9,
            seed=1, cpu=2, cache=cache)
        wall = time.time() - t0
        print(f"[egfr-calib] wall={wall:.1f} s")
        print(r.summary())
        cache.close()

        assert r.ok, f"calibration failed: {r.reason}"
        assert r.n_ok == r.n_points == 6, (
            f"expected 6/6 compounds docked; got {r.n_ok}/{r.n_points}")
        assert r.spearman_rho is not None
        # Spearman ρ should be negative — this is the honest Vina
        # kinase-failure finding. Allow up to 0.1 slack above zero in
        # case a future meeko/vina patch nudges it; anything >+0.1
        # means the known failure mode has changed and BENCHMARKS.md
        # needs re-checking.
        assert r.spearman_rho < 0.1, (
            f"Spearman ρ = {r.spearman_rho:.3f}; expected < 0.1 "
            "(known Vina-on-kinase saturation). If this now passes "
            "> 0.1 the finding changed — re-run and update "
            "BENCHMARKS.md EGFR section.")
        assert r.mae_kcalmol is not None and r.mae_kcalmol >= 1.5, (
            f"MAE = {r.mae_kcalmol:.2f} kcal/mol; expected ≥ 1.5 "
            "(Vina under-scores tight kinase binders). If this now "
            "drops below 1.5 the upstream scoring improved — "
            "re-run and update BENCHMARKS.md.")


if __name__ == "__main__":
    try:
        test_egfr_calibration_honest_failure()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
