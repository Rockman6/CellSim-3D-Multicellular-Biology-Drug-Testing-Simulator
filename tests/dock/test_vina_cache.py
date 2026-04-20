"""Cache integration — dock_ligand short-circuits on identical inputs.

Gates:
  - cold call runs Vina end-to-end (ok=True, n_poses > 0)
  - warm call (cache hit) also ok=True
  - top-pose ΔG identical across cold / warm
  - pose count preserved
  - warm wall-time < 1/20 of cold (conservative; Vina cold ~5 s,
    cache hit ~10 ms on laptop, so gate 20× is very loose)
  - ligand_inchi_key, receptor_hash preserved on the reinflated result

Run:
    conda activate cellsim
    python tests/dock/test_vina_cache.py
"""

from __future__ import annotations

import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cache import Cache  # noqa: E402
from src.dock import dock_ligand  # noqa: E402

RECEPTOR = REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb"
LIGAND = "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"


def test_vina_cache_roundtrip():
    assert RECEPTOR.exists(), f"missing {RECEPTOR}"
    with tempfile.TemporaryDirectory(prefix="cellsim-vcache-") as tmp:
        cache = Cache(Path(tmp) / "c.sqlite")

        t0 = time.time()
        r_cold = dock_ligand(
            RECEPTOR, LIGAND,
            center_A=(11.12, 1.68, -10.75),
            box_size_A=(20.0, 20.0, 20.0),
            exhaustiveness=16, num_modes=3, seed=1, cpu=2,
            cache=cache)
        cold_s = time.time() - t0
        assert r_cold.ok, f"cold failed: {r_cold.reason}"
        assert r_cold.poses, "no poses cold"
        print(f"[vcache] cold: ok={r_cold.ok} poses={len(r_cold.poses)}  "
              f"top_dG={r_cold.best_kcalmol:+.2f}  wall={cold_s:.2f} s")

        t0 = time.time()
        r_warm = dock_ligand(
            RECEPTOR, LIGAND,
            center_A=(11.12, 1.68, -10.75),
            box_size_A=(20.0, 20.0, 20.0),
            exhaustiveness=16, num_modes=3, seed=1, cpu=2,
            cache=cache)
        warm_s = time.time() - t0
        assert r_warm.ok, f"warm failed: {r_warm.reason}"
        assert len(r_warm.poses) == len(r_cold.poses)
        assert r_warm.best_kcalmol == r_cold.best_kcalmol
        assert r_warm.ligand_inchi_key == r_cold.ligand_inchi_key
        assert r_warm.receptor_hash == r_cold.receptor_hash
        print(f"[vcache] warm: ok={r_warm.ok} poses={len(r_warm.poses)}  "
              f"top_dG={r_warm.best_kcalmol:+.2f}  "
              f"wall={warm_s:.4f} s  "
              f"speedup={cold_s / max(warm_s, 0.0001):.0f}x")

        # Warm path must be dramatically faster than cold.
        assert warm_s * 20 < cold_s, (
            f"warm path not fast enough: cold={cold_s:.2f}s "
            f"warm={warm_s:.4f}s (expected warm < cold/20)")
        cache.close()


if __name__ == "__main__":
    try:
        test_vina_cache_roundtrip()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
