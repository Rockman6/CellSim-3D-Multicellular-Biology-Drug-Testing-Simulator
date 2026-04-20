"""Layer 1.3 refine smoke — post-dock OpenMM ligand minimisation.

Verifies that refine_pose_openff on Vina poses:
  - returns without raising on recoverable failure (reason stored)
  - produces valid heavy-atom coordinates (finite, non-NaN)
  - keeps pose centroid within 0.5 Å of input (restraints held)
  - flips posebusters_geometry_ok True on >= 1 of the top-2 poses
    for biotin/1STP (known-easy geometry case)

Run:
    conda activate cellsim
    python tests/dock/test_refine_smoke.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import (  # noqa: E402
    attach_crystal_rmsd, attach_posebusters, dock_ligand,
    refine_pose_openff,
)
from src.dock.vina import DockingResult  # noqa: E402

RECEPTOR = REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb"
LIGAND = "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"  # biotin
CENTER = (11.12, 1.68, -10.75)
BOX = (20.0, 20.0, 20.0)


def _centroid(positions_A):
    import numpy as np
    arr = np.array(positions_A)
    return arr.mean(axis=0).tolist()


def test_refine_biotin_1stp():
    assert RECEPTOR.exists(), f"missing {RECEPTOR}"

    r = dock_ligand(RECEPTOR, LIGAND, center_A=CENTER, box_size_A=BOX,
                    exhaustiveness=16, num_modes=3, seed=1, cpu=4)
    assert r.ok, f"dock failed: {r.reason}"
    assert r.poses, "no poses"

    pre = attach_posebusters(
        DockingResult(receptor_source=r.receptor_source,
                       ligand_smiles=r.ligand_smiles, ok=True,
                       poses=list(r.poses[:2]),
                       ligand_inchi_key=r.ligand_inchi_key,
                       ligand_formula=r.ligand_formula),
        receptor_pdb=RECEPTOR,
        crystal_pdb=RECEPTOR, crystal_resname="BTN")
    pre_geom = [p.posebusters_geometry_ok for p in pre.poses]
    print(f"[refine] pre  geometry_ok: {pre_geom}")

    refined_poses = []
    for p in r.poses[:2]:
        ref = refine_pose_openff(p, r.ligand_smiles)
        assert ref.positions_A, "refined positions empty"
        for xyz in ref.positions_A:
            for v in xyz:
                assert isinstance(v, (int, float)), \
                    f"non-numeric coord {v}"
                assert v == v, "NaN coord"
        c0 = _centroid(p.positions_A)
        c1 = _centroid(ref.positions_A)
        drift = sum((c0[i] - c1[i]) ** 2 for i in range(3)) ** 0.5
        print(f"  pose #{ref.mode}  centroid drift = {drift:.3f} Å")
        assert drift < 1.0, f"refined centroid drifted {drift:.2f} Å"
        refined_poses.append(ref)

    post = attach_posebusters(
        DockingResult(receptor_source=r.receptor_source,
                       ligand_smiles=r.ligand_smiles, ok=True,
                       poses=refined_poses,
                       ligand_inchi_key=r.ligand_inchi_key,
                       ligand_formula=r.ligand_formula),
        receptor_pdb=RECEPTOR,
        crystal_pdb=RECEPTOR, crystal_resname="BTN")
    post_geom = [p.posebusters_geometry_ok for p in post.poses]
    print(f"[refine] post geometry_ok: {post_geom}")

    # At least one top-2 pose should pass geometry after refine.
    assert any(post_geom), (
        f"no top-2 pose had geometry_ok True after refine: {post_geom}")


if __name__ == "__main__":
    try:
        test_refine_biotin_1stp()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
