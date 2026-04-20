"""Layer 1.4 DFT smoke — PySCF B3LYP/def2-SVP on small molecules.

Verifies src.quantum.dft on a minimal molecule (methane) where
the DFT numbers are well-known from the literature.

Gates:
  - ok=True, SCF converged
  - total_energy negative and finite
  - HOMO < 0 (bound molecule)
  - LUMO > HOMO
  - gap in (0, 15) eV
  - dipole finite

Uses methane (CH4) — cheap (~2 s), well-characterised DFT reference:
  HOMO ~ -10.5 eV, LUMO ~ +0.7 eV, dipole = 0.0 D (symmetric).

Run:
    conda activate cellsim
    python tests/quantum/test_dft_smoke.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.quantum import dft_single_point  # noqa: E402


def test_dft_methane():
    t0 = time.time()
    r = dft_single_point("C")
    dt = time.time() - t0
    print(f"[dft] methane  wall={dt:.1f} s")
    print(f"  {r.summary()}")

    assert r.ok, f"DFT failed: {r.reason}"
    assert r.total_energy_Hartree is not None and r.total_energy_Hartree < 0
    assert r.total_energy_eV is not None and r.total_energy_eV < 0
    assert r.homo_eV is not None and r.homo_eV < 0.0, \
        f"HOMO {r.homo_eV} not < 0"
    assert r.lumo_eV is not None and r.lumo_eV > r.homo_eV, \
        f"LUMO {r.lumo_eV} not > HOMO {r.homo_eV}"
    assert 0.0 < r.homo_lumo_gap_eV < 15.0, \
        f"gap {r.homo_lumo_gap_eV} out of (0, 15)"
    assert r.dipole_Debye is not None, "dipole missing"
    assert r.dipole_Debye < 0.5, (
        f"methane dipole {r.dipole_Debye:.2f} D should be ~ 0 "
        "(Td symmetry)")
    assert r.n_atoms == 5


def test_dft_acetic_acid():
    r = dft_single_point("CC(=O)O")
    print(f"  {r.summary()}")
    assert r.ok
    assert r.dipole_Debye is not None and r.dipole_Debye > 0.5, (
        "acetic acid should have a real dipole")


if __name__ == "__main__":
    try:
        test_dft_methane()
        test_dft_acetic_acid()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
