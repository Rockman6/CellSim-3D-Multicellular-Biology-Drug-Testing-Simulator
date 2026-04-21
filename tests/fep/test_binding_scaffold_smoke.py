"""Layer 1.3 FEP Milestone B binding scaffold smoke.

Pins the Hamelberg–Gilson analytical restraint correction against
its closed-form Gaussian integral and verifies the scaffold builder
wires a valid alchemical system end-to-end (PDBFixer + Interchange
+ AbsoluteAlchemicalFactory + CustomCentroidBondForce) on a small
receptor (ubiquitin + methane, ~14k atoms).

Why both checks:
- Formula drift: if someone edits `_harmonic_restraint_free_energy_
  kcalmol` constants (kT, 1 M standard-state volume) the computed
  ΔG_R silently changes — the analytical check traps that.
- API drift: the builder depends on openff-toolkit Topology.from_pdb
  + Interchange.from_smirnoff + openmmtools.alchemy. Any of those
  APIs breaks and this gate catches it before production runs.

Kept out of the headless validators tree on purpose: this test
needs a working conda env (openff, openmmtools, pdbfixer) and is
slow enough (~20 s on CPU) that it doesn't belong in a fast-loop
unit suite. Run it before cutting a release.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.fep.binding import (  # noqa: E402
    _harmonic_restraint_free_energy_kcalmol,
    _build_complex_alchemical_system,
)


def test_harmonic_restraint_correction_matches_closed_form():
    """Cross-check `_harmonic_restraint_free_energy_kcalmol` against
    the Gaussian integral written out by hand."""
    T_K = 298.15
    k_kJ_per_nm2 = 4184.0           # 10 kcal/mol/Å²
    kT_kJ = 0.0083144626 * T_K      # 2.479 kJ/mol
    # Volume of an isotropic 3D Gaussian: (2π·kT/k)^{3/2}
    v_harm_nm3 = (2.0 * math.pi * kT_kJ / k_kJ_per_nm2) ** 1.5
    v_std_nm3 = 1.66054             # 1 M per-molecule volume
    dG_kJ = -kT_kJ * math.log(v_std_nm3 / v_harm_nm3)
    expected_kcalmol = dG_kJ / 4.184

    got = _harmonic_restraint_free_energy_kcalmol(k_kJ_per_nm2)
    assert abs(got - expected_kcalmol) < 1e-6, (
        f"correction drift: got {got:.6f}, "
        f"expected {expected_kcalmol:.6f}")


def test_correction_scales_monotonically_with_k():
    """Looser spring → smaller magnitude correction. If this breaks,
    the sign convention in the formula has flipped."""
    tight = _harmonic_restraint_free_energy_kcalmol(4184.0)
    loose = _harmonic_restraint_free_energy_kcalmol(418.4)
    # Both are negative (standard-state volume > Gaussian volume).
    # Tighter spring → harder confinement → MORE negative correction.
    assert tight < 0 and loose < 0
    assert tight < loose, (
        f"expected tight ({tight:.2f}) more negative than loose "
        f"({loose:.2f}); the spring-k monotonicity is broken")


def test_complex_alchemical_builder_on_ubiquitin():
    """Scaffold build ubiquitin + methane complex — verifies the
    PDBFixer → Topology.from_pdb → Interchange → alchemical factory
    + restraint pipeline."""
    pdb = REPO_ROOT / "benchmarks/md/1ubq.pdb"
    out = _build_complex_alchemical_system(
        "C", pdb, padding_nm=0.8)
    assert out["n_ligand_atoms"] == 5, (
        f"methane should have 5 atoms; got {out['n_ligand_atoms']}")
    # Ubiquitin has 76 residues; after PDBFixer hydrogenation we
    # expect on the order of ~1200 protein atoms (ubiquitin canonical
    # is 1231 atoms per the PDBFixer output).
    assert 1000 <= out["n_protein_atoms"] <= 1500, (
        f"unexpected protein atom count {out['n_protein_atoms']}")
    assert out["n_total_atoms"] > 10_000, (
        "solvated box should be > 10k atoms at 0.8 nm padding")
    # Ligand indices contiguous at end of protein.
    assert out["ligand_indices"] == list(
        range(out["n_protein_atoms"],
              out["n_protein_atoms"] + out["n_ligand_atoms"]))
    assert len(out["anchor_indices"]) >= 5, (
        "restraint anchor should have at least 5 Cα atoms even on "
        "a small protein — fallthrough to the nearest-10 pick "
        "isn't kicking in")


if __name__ == "__main__":
    try:
        test_harmonic_restraint_correction_matches_closed_form()
        print("[PASS] restraint correction closed-form")
        test_correction_scales_monotonically_with_k()
        print("[PASS] restraint k monotonicity")
        test_complex_alchemical_builder_on_ubiquitin()
        print("[PASS] complex builder on 1ubq + methane")
        print("3/3 PASS")
    except AssertionError as e:
        print(f"[FAIL] {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"[ERROR] {e}")
        sys.exit(2)
