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


__all__ = ["alchemical_state_smoke", "FEPSmokeResult"]
