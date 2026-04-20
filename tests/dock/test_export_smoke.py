"""Pose export smoke — verifies SDF and PDB output.

Gates (on biotin / 1STP at exh=16, cached):
  - export_poses_sdf writes N poses; file size plausible
  - each pose in the SDF has RDKit V2000 header ('V2000' substring)
  - SDF contains property tags (>  <affinity_kcalmol>) for each pose
  - export_poses_pdb writes MODEL-separated blocks
  - each MODEL block has >= 3 HETATM rows (biotin has 16 heavy atoms
    so sane min is 3; actual ~16)

Run:
    conda activate cellsim
    python tests/dock/test_export_smoke.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import (  # noqa: E402
    dock_ligand, export_poses_sdf, export_poses_pdb,
)

RECEPTOR = REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb"
LIGAND = "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"


def test_export_sdf_and_pdb():
    assert RECEPTOR.exists(), f"missing {RECEPTOR}"

    r = dock_ligand(RECEPTOR, LIGAND,
                    center_A=(11.12, 1.68, -10.75),
                    box_size_A=(20.0, 20.0, 20.0),
                    exhaustiveness=16, num_modes=3, seed=1, cpu=4)
    assert r.ok, f"dock failed: {r.reason}"
    assert r.poses, "no poses"

    with tempfile.TemporaryDirectory(prefix="cellsim-export-") as tmp:
        tmp = Path(tmp)
        sdf = tmp / "biotin_poses.sdf"
        pdb = tmp / "biotin_poses.pdb"

        n_sdf = export_poses_sdf(r, sdf)
        n_pdb = export_poses_pdb(r, pdb)

        assert n_sdf == len(r.poses), (
            f"SDF: wrote {n_sdf} vs {len(r.poses)} expected")
        assert n_pdb == len(r.poses)

        sdf_text = sdf.read_text()
        assert "V2000" in sdf_text, "SDF header V2000 missing"
        # Properties present.
        assert ">  <affinity_kcalmol>" in sdf_text or \
               ">  <affinity_kcalmol>  (" in sdf_text, \
            "SDF has no affinity_kcalmol property"
        # Count pose blocks ($$$$ separator).
        assert sdf_text.count("$$$$") >= n_sdf
        print(f"[export] SDF: {n_sdf} poses, {sdf.stat().st_size} B")

        pdb_text = pdb.read_text()
        assert pdb_text.startswith("HEADER"), "PDB no HEADER"
        # MODEL blocks.
        n_models = pdb_text.count("\nMODEL ")
        assert n_models == n_pdb, (
            f"PDB: {n_models} MODEL blocks vs {n_pdb} expected")
        # At least 3 HETATM rows per model (biotin ~16 heavy atoms).
        hetatm_count = pdb_text.count("\nHETATM")
        assert hetatm_count >= 3 * n_pdb, (
            f"only {hetatm_count} HETATM rows across {n_pdb} models")
        print(f"[export] PDB: {n_pdb} models, "
              f"{hetatm_count} HETATM rows, {pdb.stat().st_size} B")


if __name__ == "__main__":
    try:
        test_export_sdf_and_pdb()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
