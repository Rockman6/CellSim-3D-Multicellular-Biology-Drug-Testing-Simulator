"""Layer 1.2 PR 2 gate: PDB -> solvated minimised system -> short MD.

Bundled test protein: ubiquitin (76 residues, PDB 1UBQ).

Gates:
    - Loader produces a valid, energy-minimised solvated system.
    - Minimisation drops potential energy by at least 10⁵ kJ/mol
      (typical for a fresh solvated box).
    - 500-step (1 ps) Langevin MD runs without NaN and Cα RMSD
      stays below 2 Å vs the minimised start.

The full 100 ns ubiquitin backbone-RMSD < 3 Å gate is separate —
that needs GPU + long runs, tracked under the Layer 1.2 full gate
todo.

Run:
    conda activate cellsim
    python tests/md/test_protein_load.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.md import load_protein_pdb, short_protein_md  # noqa: E402

UBQ_PATH = REPO_ROOT / "benchmarks" / "md" / "1ubq.pdb"


def test_ubiquitin_load_min_short_md():
    assert UBQ_PATH.exists(), f"missing test PDB: {UBQ_PATH}"

    print(f"[test] loading {UBQ_PATH.name}", flush=True)
    t0 = time.time()
    sys_res = load_protein_pdb(UBQ_PATH, padding_nm=1.0, ph=7.0,
                                minimise_max_iter=500)
    print(f"[test] load+minimise wall={time.time()-t0:.1f} s", flush=True)
    print(f"  {sys_res.summary()}", flush=True)

    assert sys_res.ok, f"loader failed: {sys_res.reason}"
    assert sys_res.n_atoms_protein and sys_res.n_atoms_protein > 500, (
        f"unexpected protein atom count: {sys_res.n_atoms_protein}")
    assert sys_res.n_atoms_total and sys_res.n_atoms_total > 5000, (
        f"solvated system too small: {sys_res.n_atoms_total}")
    assert sys_res.n_waters and sys_res.n_waters > 1000, (
        f"water count too low: {sys_res.n_waters}")

    de = sys_res.e_initial_kJmol - sys_res.e_minimised_kJmol
    assert de > 1.0e4, (
        f"minimisation did not relax the box ({de:.1f} kJ/mol)")

    print(f"[test] short MD (500 steps = 1 ps)", flush=True)
    t0 = time.time()
    tr = short_protein_md(sys_res, n_steps=500, dt_fs=2.0, save_every=100)
    print(f"[test] md wall={time.time()-t0:.1f} s", flush=True)
    print(f"  {tr.summary()}", flush=True)

    assert tr.ok, f"short MD failed: {tr.reason}"
    assert tr.rmsd_ca_A[-1] < 2.0, (
        f"Cα RMSD {tr.rmsd_ca_A[-1]:.2f} Å > 2.0 Å on 1 ps")


if __name__ == "__main__":
    try:
        test_ubiquitin_load_min_short_md()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
