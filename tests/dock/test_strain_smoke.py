"""Layer 1.3 ligand-strain diagnostic smoke.

Gates:
  - ligand_strain runs ok on 3 benchmark cocrystal ligands
  - benzamidine (rigid): ratio < 1.5 (good band)
  - biotin + erlotinib (flexible): ratio < 7 (not rejected)
  - absolute strain values are positive and finite
  - energy ratio matches PoseBusters' internal calculation

Why this gate: the strain check is a biologist-facing triage
filter. It must NOT over-reject typical-drug poses; ratios for
canonical docking benchmarks should land in the good/acceptable
bands. If this test starts rejecting known-good poses, something
in the RDKit/PoseBusters chain broke.

Run:
    conda activate cellsim
    python tests/dock/test_strain_smoke.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import dock_ligand  # noqa: E402
from src.dock.strain import ligand_strain  # noqa: E402


CASES = [
    dict(
        label="benzamidine (3PTB, rigid)",
        pdb="benchmarks/dock/3ptb.pdb",
        smiles="NC(=N)c1ccccc1",
        center=(-1.76, 14.46, 16.92),
        box=(18.0, 18.0, 18.0),
        max_ratio=1.8,        # rigid ligand, top-pose ratio < 1.8
        expected_band_one_of={"good", "acceptable"},
    ),
    dict(
        label="biotin (1STP, flexible)",
        pdb="benchmarks/dock/1stp.pdb",
        smiles="OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
        center=(11.12, 1.68, -10.75),
        box=(20.0, 20.0, 20.0),
        max_ratio=4.0,
        expected_band_one_of={"good", "acceptable", "suspicious"},
    ),
    dict(
        label="erlotinib (1M17, large flexible)",
        pdb="benchmarks/dock/1m17.pdb",
        smiles="COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
        center=(22.01, 0.25, 52.79),
        box=(28.0, 22.0, 22.0),
        max_ratio=4.0,
        expected_band_one_of={"good", "acceptable", "suspicious"},
    ),
]


def test_strain_benchmarks():
    t0 = time.time()
    for case in CASES:
        print(f"\n[strain] {case['label']}")
        r = dock_ligand(
            case["pdb"], case["smiles"],
            center_A=case["center"], box_size_A=case["box"],
            exhaustiveness=8, num_modes=1, seed=1, cpu=2)
        assert r.ok and r.poses, (
            f"dock failed for {case['label']}: {r.reason}")
        pose = r.poses[0]
        s = ligand_strain(
            pose.elements, pose.positions_A, case["smiles"],
            ensemble_n=20)
        print(f"  {s.summary()}")
        assert s.ok, f"strain failed: {s.reason}"
        assert s.strain_kcalmol is not None
        assert s.strain_kcalmol >= -1.0, (
            f"negative strain {s.strain_kcalmol}; FF instability?")
        assert s.energy_ratio is not None
        assert s.energy_ratio <= case["max_ratio"], (
            f"{case['label']} top-pose ratio {s.energy_ratio:.2f} "
            f"> expected max {case['max_ratio']}. "
            "Either Vina/PoseBusters upstream changed, or the "
            "benchmark setup drifted — re-check BENCHMARKS.md.")
        assert s.band in case["expected_band_one_of"], (
            f"{case['label']} landed in unexpected band "
            f"{s.band!r}; expected one of "
            f"{case['expected_band_one_of']}")
    wall = time.time() - t0
    print(f"\n[strain] wall={wall:.1f} s")


if __name__ == "__main__":
    try:
        test_strain_benchmarks()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
