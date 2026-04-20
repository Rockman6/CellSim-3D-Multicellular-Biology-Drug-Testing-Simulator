"""Cache integration — parametrize_smiles must short-circuit identical physics.

Gates:
  - cold call runs end-to-end (ok=True)
  - warm call (cache hit) also ok=True
  - physical-field outputs (inchi_key, charges, ff_version, elements,
    bonds, partial_charges) identical across cold / warm
  - warm wall time < 1/50 of cold wall time (AM1-BCC is ~20 s cold,
    cache hit is ms; using 50× as a very conservative gate that
    should never false-fail even on a loaded CI runner)

Run:
    conda activate cellsim
    python tests/chem/test_parametrize_cache.py
"""

from __future__ import annotations

import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cache import Cache  # noqa: E402
from src.chem import parametrize_smiles  # noqa: E402


def test_parametrize_cache_roundtrip():
    with tempfile.TemporaryDirectory(prefix="cellsim-pcache-") as tmp:
        cache = Cache(Path(tmp) / "c.sqlite")
        smi = "CC(=O)OC1=CC=CC=C1C(=O)O"  # aspirin

        t0 = time.time()
        r_cold = parametrize_smiles(smi, cache=cache)
        cold_s = time.time() - t0
        assert r_cold.ok, f"cold failed: {r_cold.reason}"
        print(f"[cache] cold: wall={cold_s:.2f} s")

        t0 = time.time()
        r_warm = parametrize_smiles(smi, cache=cache)
        warm_s = time.time() - t0
        assert r_warm.ok, f"warm failed: {r_warm.reason}"
        print(f"[cache] warm: wall={warm_s:.4f} s  "
              f"speedup={cold_s / max(warm_s, 0.0001):.0f}x")

        # Physical-field equality.
        assert r_cold.inchi_key == r_warm.inchi_key
        assert r_cold.ff_version == r_warm.ff_version
        assert r_cold.charge_method == r_warm.charge_method
        assert r_cold.elements == r_warm.elements
        # JSON round-trip turns bond tuples into lists; compare
        # structurally.
        assert [list(b) for b in r_cold.bonds] == \
               [list(b) for b in r_warm.bonds]
        assert r_cold.partial_charges_e == r_warm.partial_charges_e
        # positions can differ by rounding at last digit after
        # JSON-round-trip; allow a tolerance.
        for a, b in zip(r_cold.positions_nm, r_warm.positions_nm):
            for i in range(3):
                assert abs(a[i] - b[i]) < 1e-9

        # Warm path must be dramatically faster than cold.
        assert warm_s * 50 < cold_s, (
            f"warm path not fast enough: "
            f"cold={cold_s:.2f}s warm={warm_s:.4f}s "
            f"(expected warm < cold/50)")
        cache.close()


if __name__ == "__main__":
    try:
        test_parametrize_cache_roundtrip()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
