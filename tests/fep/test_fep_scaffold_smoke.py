"""Layer 1.3 FEP scaffold smoke.

Regression-pins the openmmtools-alchemy primitives: if openmmtools
or alchemy breaks / drifts / removes the AbsoluteAlchemicalFactory
API, this gate fails and we see it on PR before the biologist
does.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.fep import alchemical_state_smoke  # noqa: E402


def test_alchemical_factory_builds_on_lj_fluid():
    r = alchemical_state_smoke()
    print(r.summary())
    assert r.ok, r.reason
    assert r.n_alchemical_particles == 32, (
        f"expected 32 particles; got {r.n_alchemical_particles}")
    assert r.openmmtools_version is not None


if __name__ == "__main__":
    try:
        test_alchemical_factory_builds_on_lj_fluid()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
