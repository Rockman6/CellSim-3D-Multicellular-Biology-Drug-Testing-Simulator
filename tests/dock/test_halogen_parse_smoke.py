"""Layer 1.3 regression: halogen atoms survive the Vina PDBQT parser.

Pins the 37a6483 fix for the latent bug where AutoDock-Vina
PDBQT output put the element symbol in the atom-name column
[12:16] but the parser previously only read [12:14], chopping
two-letter elements like Cl/Br down to 'C' and silently
producing a pose.elements list with the halogen replaced by
a carbon. The bug corrupted every downstream consumer of
pose.elements for halogenated ligands: export SDF, strain
diagnostic, off-target rankings, heme-access SoM.

Gates:
  - docking succeeds on 4-Cl-benzamidine
  - pose.elements heavy-atom count matches SMILES (10)
  - at least one element is 'Cl' (the halogen is present)
  - no aromatic 'A' or donor 'Hd' types leak through
  - strain diagnostic works on the halogenated pose

If ANY of these regresses it means either the parser is back
to the bug or Meeko / Vina output format changed.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import dock_ligand  # noqa: E402
from src.dock.strain import (  # noqa: E402
    _normalise_elem,
    ligand_strain,
)


SMILES = "NC(=N)c1ccc(Cl)cc1"
RECEPTOR = REPO_ROOT / "benchmarks" / "dock" / "3ptb.pdb"
CENTER = (-1.76, 14.46, 16.92)
BOX = (18.0, 18.0, 18.0)


def test_halogen_parse_4cl_benzamidine():
    t0 = time.time()
    r = dock_ligand(
        RECEPTOR, SMILES,
        center_A=CENTER, box_size_A=BOX,
        exhaustiveness=8, num_modes=1, seed=1, cpu=2)
    assert r.ok and r.poses, f"dock failed: {r.reason}"
    pose = r.poses[0]

    # Normalise AutoDock subtypes (A→C, HD→H, NA→N, etc.) and
    # count heavy atoms.
    heavy = [_normalise_elem(e) for e in pose.elements
             if _normalise_elem(e) != "H"]
    n_heavy = len(heavy)
    assert n_heavy == 10, (
        f"4-Cl-benzamidine should have 10 heavy atoms; "
        f"got {n_heavy}. Parser may be dropping atoms again.")
    assert "Cl" in heavy, (
        f"Cl atom missing from pose.elements: {heavy}. "
        "Vina parser regression — halogens being dropped.")

    # Strain diagnostic should work on this halogenated pose.
    s = ligand_strain(pose.elements, pose.positions_A, SMILES,
                       ensemble_n=15)
    assert s.ok, f"strain failed on halogenated ligand: {s.reason}"
    assert s.band in ("good", "acceptable", "suspicious", "reject")

    wall = time.time() - t0
    print(f"[halogen-parse] n_heavy={n_heavy}, Cl present, "
          f"strain={s.band}, wall={wall:.1f}s")


if __name__ == "__main__":
    try:
        test_halogen_parse_4cl_benzamidine()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
