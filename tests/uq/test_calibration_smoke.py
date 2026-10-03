"""Layer 1.7 calibration smoke — CellSim vs published K_d.

Runs cellsim.uq.calibration on the bundled streptavidin set (4
published binders with K_d spanning 10^-14 to 10^-5 M).

Gates (rewritten 2026-10, see below):
  - docking succeeds on >= 3/4 entries
  - the SATURATION signature holds: predictions span < 2 kcal/mol
    while experiment spans > 10, so absolute dG is meaningless here
  - MAE is large and finite, consistent with that saturation
  - Conformal q95 finite and positive (calibration produced a
    real number for downstream interval use)

This test used to assert Spearman >= 0.8, i.e. that docking ranks
these binders usably. That passed only because the bundled
desthiobiotin K_d was wrong by ~6 orders of magnitude (5e-5 M where
the literature and the entry's own cited source say ~1e-11 M). With
the corrected value the measured Spearman is ~0.4 and the claim does
not hold — which is exactly what Vina's tight-binder saturation
predicts, and what the reliability table already says for this class
(`do_not_trust_absolute`). The gate now pins the saturation itself,
which is a real and reproducible property, instead of a ranking
claim that was an artefact of bad reference data.

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

from cellsim.cache import Cache  # noqa: E402
from cellsim.uq import run_calibration  # noqa: E402


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
        assert r.spearman_rho is not None     # reported, deliberately not gated

        # The real, reproducible finding: Vina saturates. These four
        # compounds span ~12 kcal/mol of experimental affinity, and the
        # predictions collapse into a fraction of that range.
        preds = [p.dG_pred_kcalmol for p in r.points
                 if p.dG_pred_kcalmol is not None]
        expts = [p.dG_expt_kcalmol for p in r.points
                 if p.dG_pred_kcalmol is not None]
        pred_span = max(preds) - min(preds)
        expt_span = max(expts) - min(expts)
        print(f"  predicted span = {pred_span:.2f} kcal/mol over an "
              f"experimental span of {expt_span:.2f}")
        assert expt_span > 10.0, (
            f"this set should span a wide affinity range; got {expt_span:.1f}")
        assert pred_span < 2.0, (
            f"predictions span {pred_span:.2f} kcal/mol — Vina is expected to "
            "saturate on this class; if this widened, the scoring path changed")
        assert r.mae_kcalmol is not None
        # Saturation against a 12 kcal/mol range forces a large absolute error.
        assert 3.0 < r.mae_kcalmol < 30.0, (
            f"MAE {r.mae_kcalmol:.2f} kcal/mol is inconsistent with the "
            "saturation this class is known for")
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
