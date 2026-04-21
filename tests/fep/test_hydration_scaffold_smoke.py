"""Layer 1.3 ligand_hydration_fep scaffolded-phase smoke.

Pins the first phase of the hydration-FEP pipeline:
  SMILES → OpenFF Sage + AM1-BCC → OpenMM System →
  openmmtools.alchemy.AbsoluteAlchemicalFactory → alchemical System.

Runs NO MD — that's the Phase 2 commit. This gate just verifies
the parametrise → alchemy pipeline imports cleanly and builds a
valid alchemical system across a small spread of chemotypes
(alkane, alcohol, aromatic).

If openff-toolkit or openmmtools drifts, this fails on PR before
a biologist tries to run an actual FEP calculation.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.fep import ligand_hydration_fep  # noqa: E402


CASES = [
    ("C", 5),                  # methane: 1 C + 4 H
    ("CCO", 9),                # ethanol: 2 C + 1 O + 6 H
    ("c1ccccc1", 12),          # benzene: 6 C + 6 H
]


def test_hydration_scaffold_builds():
    for smi, expected_atoms in CASES:
        r = ligand_hydration_fep(smi)
        print(f"  {smi:>10s}: {r.summary()}")
        assert r.ok, f"scaffold failed on {smi}: {r.reason}"
        assert r.phase == "scaffolded"
        assert r.n_alchemical_atoms == expected_atoms, (
            f"{smi}: expected {expected_atoms} atoms; "
            f"got {r.n_alchemical_atoms}")
        assert r.wall_seconds is not None
        assert r.wall_seconds < 60, (
            f"{smi}: scaffold took {r.wall_seconds:.1f}s; "
            "should be sub-minute.")
        # Phase 2 hasn't shipped — these must still be None.
        assert r.dG_hydration_kcalmol is None
        assert r.dG_hydration_ci95_kcalmol is None


if __name__ == "__main__":
    try:
        test_hydration_scaffold_builds()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
