"""Layer 1.3 load_smi SDF + .smi dispatch smoke.

Gates that load_smi can read both:
  - a plain .smi file (SMILES + name per line)
  - an SDF file with molecule titles

This is the wet-lab-ingest test: biologists often get compound
lists as SDF out of ChemDraw / synthesis-workflow tools, and
CellSim should eat those directly instead of demanding a manual
smiles conversion.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.batch import load_smi  # noqa: E402


SDF_TWO_MOLECULES = """aspirin
     RDKit          2D

 13 13  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -2.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8971    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8971    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    3.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  2  0
  1  4  1  0
  2  5  1  0
  5  6  1  0
  6  7  2  0
  6  8  1  0
  7  9  1  0
  8 10  2  0
  9 11  2  0
 11 10  1  0
  7 12  1  0
 12 13  2  0
M  END
> <comment>
aspirin, acetylsalicylic acid

$$$$
benzamidine
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971   -2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1962    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1962   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
  4  5  2  0
  4  6  1  0
  5  7  1  0
  6  8  2  0
  7  9  2  0
  9  8  1  0
M  END
$$$$
"""


def test_load_smi_handles_sdf():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sdf = td / "compounds.sdf"
        sdf.write_text(SDF_TWO_MOLECULES)
        rows = load_smi(sdf)
        assert len(rows) == 2, f"expected 2 rows; got {len(rows)}"
        names = [n for _, n in rows]
        assert "aspirin" in names, names
        assert "benzamidine" in names, names
        # RDKit canonical SMILES for aspirin / benzamidine must
        # round-trip (quick sanity check).
        smis = dict(zip(names, [s for s, _ in rows]))
        from rdkit import Chem
        for name, smi in smis.items():
            mol = Chem.MolFromSmiles(smi)
            assert mol is not None, f"bad SMILES for {name}: {smi}"
        print(f"  loaded {len(rows)} from SDF: {names}")


def test_load_smi_handles_plain_smi():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        smi = td / "compounds.smi"
        smi.write_text(
            "CC(=O)OC1=CC=CC=C1C(=O)O\taspirin\n"
            "NC(=N)c1ccccc1 benzamidine\n"
            "# comment line\n"
            "c1ccccc1\n"  # bare smiles, no name
        )
        rows = load_smi(smi)
        assert len(rows) == 3
        assert rows[0][1] == "aspirin"
        assert rows[1][1] == "benzamidine"
        # Bare SMILES gets its own SMILES as the name.
        assert rows[2][1] == "c1ccccc1"
        print(f"  loaded {len(rows)} from .smi: "
              f"{[n for _, n in rows]}")


def test_load_smi_sdf_uses_indexed_name_on_missing_title():
    # SDF record without a title → should get deterministic
    # sdf_001 etc. name.
    no_title = SDF_TWO_MOLECULES.replace(
        "aspirin\n", "\n", 1).replace("benzamidine\n", "\n", 1)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sdf = td / "untitled.sdf"
        sdf.write_text(no_title)
        rows = load_smi(sdf)
        names = [n for _, n in rows]
        assert all(n.startswith("sdf_") for n in names), names
        print(f"  untitled SDF fallback names: {names}")


if __name__ == "__main__":
    tests = [
        test_load_smi_handles_sdf,
        test_load_smi_handles_plain_smi,
        test_load_smi_sdf_uses_indexed_name_on_missing_title,
    ]
    for t in tests:
        try:
            t()
            print(f"  {t.__name__} PASS")
        except AssertionError as e:
            print(f"  {t.__name__} FAIL: {e}")
            sys.exit(1)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  {t.__name__} ERROR: {e}")
            sys.exit(2)
    print("PASS")
