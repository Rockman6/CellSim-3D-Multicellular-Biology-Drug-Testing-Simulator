"""Cache smoke — round-trip write + read, stable hashes, invalidate.

Verifies src/cache/:
  - compound_hash is stable for canonical SMILES input:
      canonical SMILES "CC(=O)O" (acetic acid) and the atom-reordered
      "OC(C)=O" must hash identically (InChI Key is the same).
  - receptor_hash on 1UBQ is stable and 16 hex chars.
  - Cache.put then Cache.get round-trips a dict value.
  - invalidate_method drops all matching rows.
  - stats() returns the expected n_entries.

Run:
    conda activate cellsim
    python tests/cache/test_cache_smoke.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cache import (  # noqa: E402
    Cache, compound_hash, receptor_hash,
)
from src.cache.hashing import method_key  # noqa: E402


UBQ = REPO_ROOT / "benchmarks" / "md" / "1ubq.pdb"


def test_compound_hash_stable_under_reorder():
    h1 = compound_hash("CC(=O)O")
    h2 = compound_hash("OC(C)=O")        # same molecule, different order
    h3 = compound_hash("C(=O)(C)O")      # another ordering
    assert h1 is not None
    assert h1 == h2 == h3, (
        f"hashes diverged across SMILES orderings: "
        f"{h1} {h2} {h3}")
    print(f"[cache] acetic-acid hash = {h1}  (stable across 3 orderings)")

    # Different compound → different hash.
    h_other = compound_hash("CCO")
    assert h_other != h1
    print(f"[cache] ethanol      hash = {h_other}")


def test_receptor_hash_stable():
    assert UBQ.exists(), f"missing {UBQ}"
    h1 = receptor_hash(UBQ)
    h2 = receptor_hash(UBQ)
    assert h1 is not None and len(h1) == 16
    assert h1 == h2
    print(f"[cache] 1UBQ receptor hash = {h1}")


def test_roundtrip_and_invalidate():
    with tempfile.TemporaryDirectory(prefix="cellsim-cache-") as tmp:
        path = Path(tmp) / "cache.sqlite"
        with Cache(path, cellsim_commit="dead0000") as c:
            # write
            smi = "CC(=O)OC1=CC=CC=C1C(=O)O"
            key = compound_hash(smi)
            method = method_key("parametrize.am1bcc",
                                 "openff-2.1.0",
                                 {"seed": 1})
            payload = {"charges": [-0.3, 0.1, 0.2], "n_atoms": 21}
            c.put(key, method, payload)
            assert c.has(key, method)

            # read
            hit = c.get(key, method)
            assert hit is not None
            assert hit.value["charges"] == payload["charges"]
            assert hit.value["n_atoms"] == 21
            assert hit.cellsim_commit == "dead0000"
            print(f"[cache] round-trip ok  key={key[:8]}…  method={method}")

            # stats
            s = c.stats()
            assert s["n_entries"] == 1
            assert method in s["entries_by_method"]
            print(f"[cache] stats: {s}")

            # invalidate by method prefix
            dropped = c.invalidate_method("parametrize")
            assert dropped == 1
            assert not c.has(key, method)
            s2 = c.stats()
            assert s2["n_entries"] == 0
            print(f"[cache] invalidate_method('parametrize') "
                  f"dropped {dropped} row(s)")


if __name__ == "__main__":
    try:
        test_compound_hash_stable_under_reorder()
        test_receptor_hash_stable()
        test_roundtrip_and_invalidate()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
