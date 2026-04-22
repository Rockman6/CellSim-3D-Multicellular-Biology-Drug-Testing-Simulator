"""Layer 1.3 — alchemical free-energy perturbation (non-AI, physics).

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
λ = 0 → λ = 1 path and recover ΔG / ΔΔG from the ensemble of short
simulations. This module wraps `openmmtools.alchemy` + `pymbar` to
do that non-AI FEP path without pulling in choderalab/perses
(evaluated and not required — openmmtools primitives cover the
Milestone A/B scope).

Shipped primitives (Milestone A — hydration)
--------------------------------------------
- `alchemical_state_smoke()` — openmmtools-import sanity gate.
- `compute_hydration_dg(smiles)` — full ΔG_hyd via 11-window
  MBAR on vacuum + TIP3P decoupling; end-to-end, no callbacks.
  Calibration: FreeSolv-12 at production parameters (M5 Max run
  in flight).
- `ligand_hydration_fep(smiles)` — the earlier scaffold-phase
  primitive, kept for the smoke CI gate that pins the scaffold
  architecture. Production users want `compute_hydration_dg`.

Shipped primitives (Milestone B — binding, see `src.fep.binding`)
-----------------------------------------------------------------
- `compute_absolute_binding_dg(smiles, pdb)` — absolute ΔG_bind
  via DDM (complex decouple − solvent decouple + Hamelberg-Gilson
  restraint correction). Two builder paths:
    * amber14 (default) — OpenMM amber14-all + tip3pfb for
      protein/solvent, SMIRNOFF ligand via SMIRNOFFTemplateGenerator.
    * smirnoff (fallback) — pure Interchange.from_smirnoff.
- `compute_relative_binding_ddg(smi_A, smi_B, pdb)` — ΔΔG(A→B)
  via two independent absolute-ΔG runs + subtraction (avoids a
  perses-style atom-mapper; corrections cancel).

Biologist CLIs (see `scripts/cellsim`)
--------------------------------------
- `cellsim fep-hyd <smi>` — single-compound hydration scaffold.
- `cellsim fep-freesolv <yaml>` — FreeSolv-12 gate runner.
- `cellsim fep-binding validate <yaml>` — sub-second YAML +
  SMILES hygiene check (RDKit + Lipinski + stereo + pKa / pocket
  sanity + wall-time estimate with local-GPU detection).
- `cellsim fep-binding bench <yaml>` — batch runner with
  crash-proof incremental CSV + Ctrl-C handler + `--resume`
  + per-compound ETA.
- `cellsim fep-binding dg / ddg` — single-compound /-pair variants.
- `cellsim fep-report <path>` — autodiscovers latest run, infers
  YAML kind, emits PASS/FAIL/partial verdict + parity PNG +
  normalised CSV for downstream tooling.
- `cellsim bench-all` — one-line-per-YAML regression dashboard.

FEP at production scale requires GPU hours; the scaffold-mode
paths (build + minimise, no MD) run in seconds on CPU and are
the subject of 50+ CI smoke tests that guard every entry point.

Non-AI: every energy term comes from an explicit force-field
parameter (AMBER14 / OpenFF Sage + AM1-BCC / TIP3P) and the MBAR
free-energy estimator is closed-form from the alchemical-state
samples — no learned surrogate anywhere.
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
        # No name on the region → unsuffixed lambda_sterics /
        # lambda_electrostatics globals, which is what
        # AlchemicalState.from_system expects downstream.
        vac_region = alchemy.AlchemicalRegion(
            alchemical_atoms=list(range(n_atoms)))
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
            alchemical_atoms=list(range(n_atoms)))
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


@dataclass
class HydrationDGResult:
    """Envelope for a full hydration ΔG run.

    ΔG_hyd = −( ΔG(decouple in solvent) − ΔG(decouple in vacuum) )

    Positive ΔG_hyd = transfer from gas to water is unfavourable
    (hydrophobic solute). Negative = favourable.
    """

    smiles: str
    ok: bool
    reason: str = ""
    dG_hydration_kcalmol: Optional[float] = None
    dG_vacuum_decouple_kcalmol: Optional[float] = None
    dG_solvent_decouple_kcalmol: Optional[float] = None
    uncertainty_kcalmol: Optional[float] = None
    n_windows: Optional[int] = None
    wall_seconds: Optional[float] = None
    # GHMC acceptance per window per leg — prof's checklist
    # requirement. < 0.70 means the timestep is too large and the
    # ΔG cannot be trusted. Passed through for the report layer.
    ghmc_acceptance_vac: list = field(default_factory=list)
    ghmc_acceptance_solv: list = field(default_factory=list)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] ΔG_hyd {self.smiles}  {self.reason}"
        return (
            f"[OK]   ΔG_hyd {self.smiles}  "
            f"{self.dG_hydration_kcalmol:+.2f} ± "
            f"{self.uncertainty_kcalmol:.2f} kcal/mol  "
            f"(ΔG_vac={self.dG_vacuum_decouple_kcalmol:+.2f}, "
            f"ΔG_solv={self.dG_solvent_decouple_kcalmol:+.2f})  "
            f"wall {self.wall_seconds:.1f} s")


def _build_alchemical_legs(smiles: str, *,
                            softcore_alpha: float = 0.5):
    """Build (vac_alch, solv_alch, top_openmm_vac, solv_positions,
    vac_positions) for a SMILES. `softcore_alpha` controls how
    "flat" the decoupled softcore potential is; larger (e.g. 1.0)
    is more forgiving on unphysical endpoint transitions but
    harder for MBAR to connect states. Default 0.5 follows
    openmmtools convention."""
    import numpy as np
    from openff.toolkit.topology import Molecule
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import (
        solvate_topology,
    )
    from openff.units import unit as offunit
    from openmm import unit as ommunit
    from openmmtools import alchemy

    mol = Molecule.from_smiles(smiles)
    mol.generate_conformers(n_conformers=1)
    mol.assign_partial_charges("am1bcc")
    top = mol.to_topology()

    ff_vac = ForceField("openff-2.1.0.offxml")
    vac_system = ff_vac.create_openmm_system(
        top, charge_from_molecules=[mol])
    n_atoms = vac_system.getNumParticles()
    factory = alchemy.AbsoluteAlchemicalFactory()
    vac_region = alchemy.AlchemicalRegion(
        alchemical_atoms=list(range(n_atoms)),
        softcore_alpha=softcore_alpha)
    vac_alch = factory.create_alchemical_system(
        vac_system, vac_region)

    # 1.2 nm padding → box ≥ 2.4 nm so the default 1.0 nm non-
    # bonded cutoff is ≤ half the box (OpenMM PBC requirement).
    solv_top = solvate_topology(
        topology=top, nacl_conc=0.0 * offunit.molar,
        padding=1.2 * offunit.nanometer)
    ff_solv = ForceField("openff-2.1.0.offxml", "tip3p.offxml")
    ichg = Interchange.from_smirnoff(
        force_field=ff_solv, topology=solv_top,
        charge_from_molecules=[mol])
    solv_system = ichg.to_openmm(
        combine_nonbonded_forces=True)
    solv_region = alchemy.AlchemicalRegion(
        alchemical_atoms=list(range(n_atoms)),
        softcore_alpha=softcore_alpha)
    solv_alch = factory.create_alchemical_system(
        solv_system, solv_region)

    vac_pos_nm = np.asarray(
        mol.conformers[0].m_as(offunit.nanometer), dtype=float)
    vac_positions = vac_pos_nm * ommunit.nanometer
    solv_pos_nm = np.asarray(
        ichg.positions.m_as(offunit.nanometer), dtype=float)
    solv_positions = solv_pos_nm * ommunit.nanometer

    return (vac_alch, solv_alch,
            top.to_openmm(), solv_top.to_openmm(),
            vac_positions, solv_positions,
            n_atoms)


def compute_hydration_dg(
    smiles: str,
    *,
    n_windows: int = 11,
    n_production_steps: int = 2000,
    n_equilibration_steps: int = 500,
    sample_stride: int = 100,
    seed: int = 1,
    softcore_alpha: float = 0.5,
) -> HydrationDGResult:
    """End-to-end absolute hydration free energy.

    Runs Langevin+MBAR alchemical decoupling on both the bare
    ligand (vacuum) and the TIP3P-solvated system, then combines:

        ΔG_hyd = −( ΔG_solv_decouple − ΔG_vacuum_decouple )

    At the default parameters (11 × 2000 steps = 4 ps per window,
    44 ps total MD per leg, CPU) this runs in roughly 5-10 min
    per compound. Accuracy is not publication-grade — that needs
    ≥ 100 ps per window with REST2 or replica exchange — but the
    sign and rough magnitude are right. The Milestone A gate
    against a FreeSolv subset will run with longer windows on a
    rented GPU via workflow_dispatch.
    """
    import time

    t0 = time.time()
    result = HydrationDGResult(
        smiles=smiles, ok=False, n_windows=n_windows)

    try:
        (vac_alch, solv_alch,
         vac_top, solv_top,
         vac_pos, solv_pos, _n) = _build_alchemical_legs(
            smiles, softcore_alpha=softcore_alpha)
    except Exception as e:
        result.reason = f"alchemical build failed: {str(e)[:200]}"
        result.wall_seconds = time.time() - t0
        return result

    from src.fep.sampling import sample_alchemical_windows

    try:
        vac_r = sample_alchemical_windows(
            vac_alch, vac_top, vac_pos,
            n_windows=n_windows,
            n_equilibration_steps=n_equilibration_steps,
            n_production_steps=n_production_steps,
            sample_stride=sample_stride, seed=seed)
        if not vac_r.ok:
            result.reason = f"vacuum leg failed: {vac_r.reason}"
            result.wall_seconds = time.time() - t0
            return result
        solv_r = sample_alchemical_windows(
            solv_alch, solv_top, solv_pos,
            n_windows=n_windows,
            n_equilibration_steps=n_equilibration_steps,
            n_production_steps=n_production_steps,
            sample_stride=sample_stride, seed=seed)
        if not solv_r.ok:
            result.reason = f"solvent leg failed: {solv_r.reason}"
            result.wall_seconds = time.time() - t0
            return result
    except Exception as e:
        result.reason = f"sampling failed: {str(e)[:200]}"
        result.wall_seconds = time.time() - t0
        return result

    # ΔG_hyd = −(ΔG_solv_decouple − ΔG_vac_decouple). A positive
    # value means gas → water transfer is unfavourable (hydro-
    # phobic). For methane the experimental value is +2.0
    # kcal/mol.
    import math
    ddg = solv_r.dG_kcalmol - vac_r.dG_kcalmol
    result.dG_hydration_kcalmol = -ddg
    result.dG_vacuum_decouple_kcalmol = vac_r.dG_kcalmol
    result.dG_solvent_decouple_kcalmol = solv_r.dG_kcalmol
    result.uncertainty_kcalmol = math.sqrt(
        (vac_r.dG_uncertainty_kcalmol or 0.0) ** 2
        + (solv_r.dG_uncertainty_kcalmol or 0.0) ** 2)
    result.ghmc_acceptance_vac = list(vac_r.ghmc_acceptance)
    result.ghmc_acceptance_solv = list(solv_r.ghmc_acceptance)
    result.ok = True
    result.wall_seconds = time.time() - t0
    return result


# Lazy re-export of binding primitives: eager `from src.fep.binding
# import ...` at module-init pulls in binding.py before anyone can
# run `python -m src.fep.binding` as __main__, producing a
# RuntimeWarning about the module being in sys.modules prior to
# execution. Using __getattr__ + a local cache defers the import
# until a caller actually accesses `from src.fep import
# compute_absolute_binding_dg` (which is still supported).
_LAZY_BINDING_ATTRS = {
    "BindingDGResult",
    "RelativeBindingDDGResult",
    "compute_absolute_binding_dg",
    "compute_relative_binding_ddg",
}


def __getattr__(name: str):
    if name in _LAZY_BINDING_ATTRS:
        from src.fep import binding as _b
        val = getattr(_b, name)
        globals()[name] = val   # cache; single import after first use
        return val
    raise AttributeError(
        f"module 'src.fep' has no attribute {name!r}")


__all__ = [
    "alchemical_state_smoke",
    "FEPSmokeResult",
    "ligand_hydration_fep",
    "HydrationFEPResult",
    "compute_hydration_dg",
    "HydrationDGResult",
    "compute_absolute_binding_dg",
    "compute_relative_binding_ddg",
    "BindingDGResult",
    "RelativeBindingDDGResult",
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
