"""Layer 1.3 prep smoke test — verifies the Meeko prep pipelines.

This is a prep-only test: it produces a receptor PDBQT from 1UBQ and
a ligand PDBQT from aspirin. Full end-to-end docking (with RMSD-vs-
crystal gates) lands in test_redocking.py once a cocrystal PDB is
bundled.

Requires the `cellsim` conda env (vina, meeko, rdkit).

Run:
    conda activate cellsim
    python tests/dock/test_prep_smoke.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.vina import _prep_receptor_pdbqt, _prep_ligand_pdbqt  # noqa: E402

UBQ_PATH = REPO_ROOT / "benchmarks" / "md" / "1ubq.pdb"


def test_prep_pipelines():
    assert UBQ_PATH.exists(), f"missing test PDB: {UBQ_PATH}"
    with tempfile.TemporaryDirectory(prefix="cellsim-dock-smoke-") as tmp:
        workdir = Path(tmp)
        print(f"[prep] preparing aspirin ligand", flush=True)
        lig = _prep_ligand_pdbqt(
            "CC(=O)OC1=CC=CC=C1C(=O)O", workdir, seed=1)
        size = lig.stat().st_size
        assert size > 500, f"ligand PDBQT too small ({size} B)"
        print(f"  ligand PDBQT ok ({size} B)", flush=True)

        print(f"[prep] preparing 1UBQ receptor", flush=True)
        rec = _prep_receptor_pdbqt(UBQ_PATH, workdir)
        size = rec.stat().st_size
        assert size > 10_000, f"receptor PDBQT too small ({size} B)"
        print(f"  receptor PDBQT ok ({size} B)", flush=True)


if __name__ == "__main__":
    try:
        test_prep_pipelines()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
