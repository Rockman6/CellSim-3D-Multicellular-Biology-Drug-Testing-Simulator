"""Layer 1.3 — alchemical free-energy scaffold (non-AI, physics).

Why this module exists
----------------------
The EGFR calibration bundle (`benchmarks/dock/egfr_calibration.yaml`)
documents a hard, honest result: Vina's empirical scoring function
does NOT reliably rank-order kinase ATP-site inhibitors (Spearman
-0.49). The ΔG is good enough for wide-dynamic-range triage on
non-kinase pockets (streptavidin, trypsin) but not for close
analogs in a congeneric kinase series — a workflow at the heart of
oncology drug discovery.

The physics-legitimate fix is alchemical free-energy perturbation
(FEP): integrate a softcore-enabled ligand through a
λ = 0 → λ = 1 path and recover ΔΔG from the ensemble of short
simulations. This module wraps `openmmtools.alchemy` to do that
non-AI FEP path without pulling in choderalab/perses (which the
conda env doesn't currently ship cleanly — see environment.yml).

What's here now
---------------
- `alchemical_state_smoke()` — sanity check that the openmmtools
  alchemical factory can build a valid AbsoluteAlchemicalFactory
  on a trivial OpenMM test system. This is the "the building
  blocks import and work on this machine" gate.

What's NOT here yet (named open items for subsequent PRs)
---------------------------------------------------------
- `ligand_hydration_fep(smiles)` — compute absolute hydration ΔG
  by a 12-window alchemical transformation (vacuum → TIP3P
  solvation). Calibration: published hydration free energies on
  a small-molecule subset of FreeSolv.

- `relative_binding_fep(smiles_A, smiles_B, receptor_pdb, box)` —
  compute ΔΔG(A→B) for a congeneric-series pair. This is what
  will rescue the EGFR kinase ranking failure.

- `fep_batch_rescore(batch_csv)` — accept a `cellsim dock` output
  CSV and emit a rescored CSV with FEP ΔΔG for the top-K hits.

FEP requires real GPU time (hours, not seconds). Batch-rescore
will run locally on small N and delegate to a cloud H100 for
scale screens; compute path TBD with the Layer 1.7 blind-bench
GPU runner.

Non-AI: every energy term in the softcore potential comes from
an explicit force-field parameter (ff14SB / OpenFF Sage) and the
MBAR / BAR free-energy estimator is closed-form from the
alchemical-state samples — no learned surrogate anywhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class FEPSmokeResult:
    ok: bool
    reason: str = ""
    openmmtools_version: Optional[str] = None
    n_alchemical_particles: Optional[int] = None

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] FEP smoke  {self.reason}"
        return (
            f"[OK]   FEP smoke  openmmtools="
            f"{self.openmmtools_version}  "
            f"alchemical_particles={self.n_alchemical_particles}  "
            "— alchemy primitives importable and factory works")


def alchemical_state_smoke() -> FEPSmokeResult:
    """Verify that openmmtools.alchemy can build an alchemical
    factory on a canonical test system. Does NOT run dynamics —
    that's a later commit. This is the 'env is sane for FEP' gate.

    Uses `openmmtools.testsystems.LennardJonesFluid` as the minimal
    non-trivial OpenMM system (pure LJ particles in a box). An
    alchemical region annotates one particle as softcore; the
    factory build step is the load-bearing verification.
    """
    try:
        import openmmtools
        from openmmtools import alchemy, testsystems
    except ImportError as e:
        return FEPSmokeResult(
            ok=False, reason=f"openmmtools import failed: {e}")

    try:
        tsys = testsystems.LennardJonesFluid(
            nparticles=32, reduced_density=0.5)
        region = alchemy.AlchemicalRegion(
            alchemical_atoms=[0], name="ligand")
        factory = alchemy.AbsoluteAlchemicalFactory()
        alch_system = factory.create_alchemical_system(
            reference_system=tsys.system,
            alchemical_regions=region)
        # Sanity: the resulting system should have the same atom
        # count; one extra CustomNonbondedForce is added for
        # softcore terms.
        n_atoms = alch_system.getNumParticles()
    except Exception as e:
        return FEPSmokeResult(
            ok=False,
            reason=f"AbsoluteAlchemicalFactory build failed: "
                   f"{str(e)[:200]}",
            openmmtools_version=openmmtools.__version__)

    return FEPSmokeResult(
        ok=True,
        openmmtools_version=openmmtools.__version__,
        n_alchemical_particles=n_atoms)


@dataclass
class HydrationFEPResult:
    """Envelope for an absolute-hydration ΔG FEP run.

    The full calculation runs a multi-window alchemical MD
    protocol vacuum → TIP3P and estimates ΔG_hyd via MBAR. This
    envelope is the stable return type for the CLI; it can carry
    a partial / scaffolding-only result too.
    """

    smiles: str
    ok: bool
    reason: str = ""
    dG_hydration_kcalmol: Optional[float] = None
    dG_hydration_ci95_kcalmol: Optional[float] = None
    n_windows: Optional[int] = None
    n_alchemical_atoms: Optional[int] = None
    wall_seconds: Optional[float] = None
    phase: Optional[str] = None     # 'scaffolded' | 'sampled'

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] ΔG_hyd {self.smiles}  {self.reason}"
        if self.dG_hydration_kcalmol is None:
            return (
                f"[OK]   ΔG_hyd {self.smiles}  "
                f"phase={self.phase}  "
                f"n_alchemical_atoms={self.n_alchemical_atoms}  "
                "(MD sampling not yet implemented)")
        return (
            f"[OK]   ΔG_hyd {self.smiles}  "
            f"{self.dG_hydration_kcalmol:+.2f} ± "
            f"{self.dG_hydration_ci95_kcalmol:.2f} kcal/mol  "
            f"wall={self.wall_seconds:.1f} s")


def ligand_hydration_fep(
    smiles: str,
    *,
    n_windows: int = 12,
    parametrize_timeout_s: int = 120,
) -> HydrationFEPResult:
    """Set up an absolute-hydration FEP for `smiles`.

    Phase 1 (this commit — 'scaffolded'): build the OpenFF-
    parametrised OpenMM system for the isolated ligand, wrap it
    with openmmtools.alchemy.AbsoluteAlchemicalFactory, and
    return the alchemical system metadata. Verifies the chem →
    openmm → alchemy pipeline end-to-end without running any MD.

    Phase 2 (next commit — 'sampled'): pair the scaffold with a
    solvated (TIP3P) alchemical system, run n_windows of short
    Langevin MD per leg, MBAR the free-energy difference. Gated
    against FreeSolv ΔG_hyd values on a small held-out subset.

    Returning phase='scaffolded' ok=True on phase-1 success is
    intentional: it lets the CLI + smoke tests treat the scaffold
    as a valid intermediate deliverable without silently claiming
    a ΔG number we haven't computed yet.
    """
    import time
    t0 = time.time()

    result = HydrationFEPResult(
        smiles=smiles, ok=False, phase="scaffolded",
        n_windows=n_windows)

    # Step 1: build an OpenMM System for the bare ligand using
    # OpenFF Sage. Done inline (not via src.chem.parametrize)
    # because parametrize returns structured data — charges,
    # positions, elements — not the live OpenMM System object an
    # alchemical factory needs. Going direct keeps this pipeline
    # short and keeps the non-FEP ADMET path unchanged.
    try:
        from openff.toolkit.topology import Molecule
        from openff.toolkit.typing.engines.smirnoff import \
            ForceField
    except ImportError as e:
        result.reason = f"openff import failed: {e}"
        result.wall_seconds = time.time() - t0
        return result

    try:
        mol = Molecule.from_smiles(smiles)
        mol.generate_conformers(n_conformers=1)
        # AM1-BCC charges — same method src.chem.parametrize
        # uses for docking. Cached internally by the toolkit.
        mol.assign_partial_charges("am1bcc")
        top = mol.to_topology()
        ff = ForceField("openff-2.1.0.offxml")
        vac_system = ff.create_openmm_system(
            top, charge_from_molecules=[mol])
    except Exception as e:
        result.reason = (
            f"OpenFF parametrise failed: {str(e)[:200]}")
        result.wall_seconds = time.time() - t0
        return result

    # Step 2: build an alchemical region covering every ligand
    # atom (absolute hydration FEP — whole molecule decouples
    # over λ).
    try:
        from openmmtools import alchemy
    except ImportError as e:
        result.reason = f"openmmtools import failed: {e}"
        result.wall_seconds = time.time() - t0
        return result

    n_atoms = vac_system.getNumParticles()
    try:
        vac_region = alchemy.AlchemicalRegion(
            alchemical_atoms=list(range(n_atoms)),
            name=f"ligand_{smiles[:16]}")
        factory = alchemy.AbsoluteAlchemicalFactory()
        vac_alch = factory.create_alchemical_system(
            reference_system=vac_system,
            alchemical_regions=vac_region)
    except Exception as e:
        result.reason = (
            f"alchemical factory (vacuum) failed: "
            f"{str(e)[:200]}")
        result.wall_seconds = time.time() - t0
        return result

    # Phase-2a: solvate via openff.interchange + packmol and
    # rebuild the alchemical system on the solvated topology.
    # Hydration free energy is:
    #   ΔG_hyd = ΔG(decouple in solvent) − ΔG(decouple in vacuum)
    # so both legs are required. Using Interchange sidesteps the
    # openff-toolkit 0.16 Modeller route that was failing with
    # a pint NoneType error — solvate_topology builds the water
    # box in OpenFF-native units, then Interchange.from_smirnoff
    # parametrises ligand (Sage) + water (tip3p.offxml) together.
    try:
        from openff.interchange import Interchange
        from openff.interchange.components._packmol import (
            solvate_topology,
        )
        from openff.units import unit as offunit

        solv_top = solvate_topology(
            topology=top,
            nacl_conc=0.0 * offunit.molar,
            padding=0.9 * offunit.nanometer,
        )
        solv_ff = ForceField("openff-2.1.0.offxml", "tip3p.offxml")
        ichg = Interchange.from_smirnoff(
            force_field=solv_ff, topology=solv_top,
            charge_from_molecules=[mol])
        solv_system = ichg.to_openmm(
            combine_nonbonded_forces=True)

        # Alchemical region covers only the ligand atoms (first
        # n_atoms of the solvated topology; packmol places the
        # solute at the front).
        solv_region = alchemy.AlchemicalRegion(
            alchemical_atoms=list(range(n_atoms)),
            name=f"ligand_{smiles[:16]}_solvated")
        _solv_alch = factory.create_alchemical_system(
            reference_system=solv_system,
            alchemical_regions=solv_region)
    except Exception as e:
        result.reason = (
            f"solvated-system build failed: {str(e)[:200]}")
        # Vacuum leg is still valid — be honest about partial
        # progress instead of failing the whole call.
        result.ok = True
        result.n_alchemical_atoms = n_atoms
        result.phase = "scaffolded_vacuum_only"
        result.wall_seconds = time.time() - t0
        return result

    # Phase-2a success: both legs built. Phase-2b (MD + MBAR +
    # FreeSolv gate) follows in a dedicated commit.
    result.ok = True
    result.n_alchemical_atoms = n_atoms
    result.phase = "scaffolded_both_legs"
    result.wall_seconds = time.time() - t0
    return result


__all__ = [
    "alchemical_state_smoke",
    "FEPSmokeResult",
    "ligand_hydration_fep",
    "HydrationFEPResult",
]


def main(argv: Optional[list] = None) -> int:
    """CLI: `cellsim fep-hyd <SMILES>` — run the hydration FEP
    scaffold and print the result summary.

    Currently runs Phase 1 only (SMILES → alchemical system); the
    MD sampling phase is the named next commit.
    """
    import argparse
    import json as _json
    import sys

    ap = argparse.ArgumentParser(
        description="Build an absolute-hydration FEP alchemical "
                    "system for a SMILES. MD sampling is not yet "
                    "implemented; this is the scaffold-phase "
                    "Layer 1.3 FEP tool.")
    ap.add_argument("smiles")
    ap.add_argument("--n-windows", type=int, default=12,
                    help="alchemical λ-windows for the future "
                         "MD sampling phase (default 12)")
    ap.add_argument("--json", action="store_true",
                    help="emit the full result as a JSON object "
                         "(no summary line) for programmatic "
                         "consumption — useful when scripting "
                         "fep-hyd across many compounds.")
    args = ap.parse_args(argv)

    r = ligand_hydration_fep(args.smiles, n_windows=args.n_windows)
    if args.json:
        from dataclasses import asdict
        print(_json.dumps(asdict(r), indent=2, default=str))
    else:
        print(r.summary())
    return 0 if r.ok else 1


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main())
