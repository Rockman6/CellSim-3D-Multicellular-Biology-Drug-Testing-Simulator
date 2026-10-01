"""Guard the docked-pose reconstruction (Meeko atom-order fix).

Meeko REORDERS atoms during ligand prep, so mapping a pose's PDBQT-order
coordinates onto a `MolFromSmiles` template by atom index builds a
topologically SCRAMBLED molecule — physically impossible 3-5 Å "bonds"
— which silently corrupted both PoseBusters validity and symmetry-aware
RMSD (they read a scrambled mol). The fix reconstructs each pose via
Meeko's reverse conversion and stores it on `pose.rdkit_molblock`.

This test docks a ligand KNOWN to be reordered by Meeko (biotin) and
asserts the reconstructed molecule is chemically sane. Without the fix,
~12/17 heavy-atom bonds are impossible and this fails loudly.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

PDB = REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb"
BIOTIN = "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"
# fpocket/crystal biotin pocket centre for 1stp (from mini_bench.yaml).
CENTER = (11.12, 1.68, -10.75)
BOX = (20.0, 20.0, 20.0)


def _dock_top_pose():
    from src.dock import dock_ligand
    r = dock_ligand(PDB, BIOTIN, center_A=CENTER, box_size_A=BOX,
                    exhaustiveness=8, num_modes=5, seed=1, cpu=2)
    assert r.ok, r.reason
    return r.poses[0]


def test_pose_carries_reconstructed_molblock():
    pose = _dock_top_pose()
    assert pose.rdkit_molblock, (
        "pose has no rdkit_molblock — Meeko reverse reconstruction "
        "silently failed; PoseBusters + RMSD would fall back to the "
        "scrambling template path")


def test_reconstructed_mol_has_sane_bond_geometry():
    """The molecule PoseBusters/RMSD actually see must have physically
    plausible heavy-atom bond lengths (no atom-order scramble)."""
    import numpy as np
    from rdkit import Chem
    from src.dock.validity import _build_pose_mol

    pose = _dock_top_pose()
    mol = _build_pose_mol(pose, BIOTIN)
    assert mol is not None and mol.GetNumConformers() > 0
    conf = mol.GetConformer()
    bad = []
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        d = float(np.linalg.norm(
            np.array(conf.GetAtomPosition(b.GetBeginAtomIdx()))
            - np.array(conf.GetAtomPosition(b.GetEndAtomIdx()))))
        # Real covalent heavy-atom bonds are ~1.2-1.9 Å (C-S ~1.81).
        if not (1.1 <= d <= 1.95):
            bad.append((round(d, 2), a1.GetSymbol(), a2.GetSymbol()))
    assert not bad, (
        f"{len(bad)} physically impossible heavy-atom bond(s) — the "
        f"pose→mol atom mapping is scrambled: {bad}")


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
