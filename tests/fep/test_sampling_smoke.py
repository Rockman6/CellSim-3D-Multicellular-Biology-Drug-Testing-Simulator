"""Layer 1.3 FEP sampling + MBAR smoke.

Verifies the Langevin + MBAR pipeline on the smallest possible
test: decoupling neutral methane in vacuum should return
ΔG ≈ 0 (no self-interactions change with charge/LJ scaling for
a single saturated hydrocarbon).

The test is sub-minute on CPU and pins that:
  - AlchemicalState.from_system succeeds on our unsuffixed
    alchemical region
  - sample_alchemical_windows runs MD + evaluates u_kn without
    crashing
  - pymbar.MBAR returns a finite free-energy difference
  - that value lies within a tight tolerance of zero for the
    methane vacuum reference (physics-sanity check)
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


def _build_methane_vacuum_alchemical():
    import numpy as np
    from openff.toolkit.topology import Molecule
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.units import unit as offunit
    from openmm import unit as ommunit
    from openmmtools import alchemy

    mol = Molecule.from_smiles("C")
    mol.generate_conformers(n_conformers=1)
    mol.assign_partial_charges("am1bcc")
    top = mol.to_topology()
    ff = ForceField("openff-2.1.0.offxml")
    sys_ = ff.create_openmm_system(
        top, charge_from_molecules=[mol])

    region = alchemy.AlchemicalRegion(
        alchemical_atoms=list(range(sys_.getNumParticles())))
    factory = alchemy.AbsoluteAlchemicalFactory()
    alch = factory.create_alchemical_system(sys_, region)

    pos_nm = np.asarray(mol.conformers[0].m_as(offunit.nanometer))
    return (alch, top.to_openmm(),
            pos_nm * ommunit.nanometer)


def test_methane_vacuum_decouple_is_near_zero():
    from src.fep.sampling import sample_alchemical_windows

    alch, top, positions = _build_methane_vacuum_alchemical()
    r = sample_alchemical_windows(
        alch, top, positions,
        n_windows=6, n_equilibration_steps=100,
        n_production_steps=500, sample_stride=50,
        seed=1)
    print(f"  {r.summary()}")
    assert r.ok, r.reason
    assert r.dG_kcalmol is not None
    # Neutral saturated methane vacuum decoupling should give
    # ΔG ≈ 0 (no intramolecular-interaction change with λ).
    # Loose tolerance because the sampling protocol is smoke-
    # size (6 × 10 samples); a real calc would be ≥ 50 ps per
    # window.
    assert abs(r.dG_kcalmol) < 1.0, (
        f"methane vacuum ΔG = {r.dG_kcalmol:.2f} kcal/mol; "
        "expected ≈ 0. Something broke in the sampling + MBAR "
        "pipeline (see AlchemicalState wiring or pymbar call).")


if __name__ == "__main__":
    try:
        test_methane_vacuum_decouple_is_near_zero()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
