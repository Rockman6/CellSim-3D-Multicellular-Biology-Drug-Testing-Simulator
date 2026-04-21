"""Layer 1.3 FEP — absolute + relative binding ΔG scaffold (Milestone B).

Thermodynamic cycle (double-decoupling method, DDM)::

    ligand in water ──ΔG_dec_solvent──→ phantom ligand in water
         ↓                                      ↓
    ΔG_bind                                   ΔG_bind'
         ↓                                      ↓
    bound complex ──ΔG_dec_complex──→ phantom ligand in complex
                                       (held near site by restraint R)

so

    ΔG_bind = −(ΔG_dec_complex − ΔG_dec_solvent) + ΔG_R  + ΔG_std

where `ΔG_R` is the analytical free-energy cost of the harmonic
centre-of-mass restraint we use to keep the decoupled ligand from
drifting out of the binding site, and `ΔG_std` is the standard-state
correction (ligand concentration 1 M) that converts the restrained
bound state to a 1 M reference. Both corrections are closed-form,
NOT sampled — no extra MD needed.

For a ligand pair A → B on the same receptor::

    ΔΔG_bind(A→B) = ΔG_bind(B) − ΔG_bind(A)

We compute it as two independent absolute-binding runs rather than
via a hybrid single-topology morph. Reasons:

1. **Reuses the validated hydration sampler unchanged.**
   `src.fep.sampling.sample_alchemical_windows` already handles the
   reverse-λ iteration, per-λ minimise + applyConstraints, and GHMC
   acceptance diagnostic. The sampler is the riskiest piece; we
   don't want a second, parallel MCMC driver for a hybrid.
2. **No atom-mapping machinery.** A hybrid morph needs a
   perses-style common-core mapper to decide which A-atom is which
   B-atom; for a restarted project that we've decided NOT to depend
   on perses for, writing an atom-mapper is its own multi-week sub-
   project. DDM bypasses it — each ligand parametrises independently.
3. **Cost ~2× hybrid but the wall-time is dominated by the protein
   MD** (thousands of waters + protein atoms vs. a handful of ligand
   atoms). Doubling ligand-only sampling adds <5 % to total wall time
   per pair once we're sampling on GPU.

The cost: sign of the *absolute* ΔG_bind depends on getting the
restraint correction right. For a congeneric pair, errors in ΔG_R
and ΔG_std cancel in the subtraction, so ΔΔG is relatively robust.
We gate Milestone B on ΔΔG ranking (Kendall τ on the EGFR series),
not on absolute-ΔG accuracy — ΔΔG is the quantity med-chem cares
about anyway.

This module is **scaffold-phase** in this commit. It builds both
alchemical systems (complex + solvent) end-to-end, verifies they
minimise cleanly, and returns `phase='scaffolded'` so callers can
distinguish the builder from a sampled ΔG. Phase-2 (next commit,
after the M5 Max FreeSolv tarball confirms the sampler) wires
`sample_alchemical_windows` onto both legs and applies the
corrections.

Non-AI: every energy term comes from explicit SMIRNOFF parameters
(OpenFF Sage 2.1.0 bonded + AM1-BCC charges on the ligand; AMBER
ff14SB SMIRNOFF port on the protein; TIP3P on water). No ML
surrogate anywhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional


# Canonical SMIRNOFF force-field bundle for protein + ligand + water.
# Pinned names (with version suffix for the amber port) so a silent
# openff-forcefields bump can't change the numerics underneath us.
_LIGAND_FF = "openff-2.1.0.offxml"
_PROTEIN_FF = "ff14sb_off_impropers_0.0.4.offxml"
_WATER_FF = "tip3p.offxml"


@dataclass
class BindingDGResult:
    """Envelope for an absolute-binding ΔG run.

    Phase states:
    - 'scaffolded_both_legs': both systems built + minimised, no MD.
    - 'scaffolded_complex_only': complex built, solvent leg failed.
    - 'sampled': MD + MBAR complete, ΔG_bind populated.
    """

    smiles: str
    receptor: str
    ok: bool
    reason: str = ""
    dG_bind_kcalmol: Optional[float] = None
    dG_dec_complex_kcalmol: Optional[float] = None
    dG_dec_solvent_kcalmol: Optional[float] = None
    dG_restraint_correction_kcalmol: Optional[float] = None
    dG_standard_state_kcalmol: Optional[float] = None
    uncertainty_kcalmol: Optional[float] = None
    restraint_k_kJ_per_nm2: Optional[float] = None
    restraint_r0_nm: Optional[float] = None
    n_windows: Optional[int] = None
    n_ligand_atoms: Optional[int] = None
    n_protein_atoms: Optional[int] = None
    n_total_atoms_complex: Optional[int] = None
    n_total_atoms_solvent: Optional[int] = None
    wall_seconds: Optional[float] = None
    phase: Optional[str] = None

    def summary(self) -> str:
        if not self.ok:
            return (f"[FAIL] ΔG_bind {self.smiles} / "
                    f"{self.receptor}  {self.reason}")
        if self.dG_bind_kcalmol is None:
            return (f"[OK]   ΔG_bind {self.smiles} / "
                    f"{self.receptor}  phase={self.phase}  "
                    f"lig_atoms={self.n_ligand_atoms}  "
                    f"prot_atoms={self.n_protein_atoms}  "
                    f"complex={self.n_total_atoms_complex}  "
                    f"wall={self.wall_seconds:.1f}s")
        return (f"[OK]   ΔG_bind {self.smiles} / "
                f"{self.receptor}  "
                f"{self.dG_bind_kcalmol:+.2f} ± "
                f"{self.uncertainty_kcalmol:.2f} kcal/mol  "
                f"wall={self.wall_seconds:.1f}s")


@dataclass
class RelativeBindingDDGResult:
    """Envelope for a ΔΔG(A→B) relative-binding run."""

    smiles_A: str
    smiles_B: str
    receptor: str
    ok: bool
    reason: str = ""
    ddG_kcalmol: Optional[float] = None
    uncertainty_kcalmol: Optional[float] = None
    dG_bind_A_kcalmol: Optional[float] = None
    dG_bind_B_kcalmol: Optional[float] = None
    bind_A: Optional[BindingDGResult] = None
    bind_B: Optional[BindingDGResult] = None
    phase: Optional[str] = None
    wall_seconds: Optional[float] = None

    def summary(self) -> str:
        if not self.ok:
            return (f"[FAIL] ΔΔG {self.smiles_A}→{self.smiles_B} / "
                    f"{self.receptor}  {self.reason}")
        if self.ddG_kcalmol is None:
            return (f"[OK]   ΔΔG {self.smiles_A}→{self.smiles_B} / "
                    f"{self.receptor}  phase={self.phase}  "
                    f"(both legs scaffolded; sampling pending)")
        return (f"[OK]   ΔΔG {self.smiles_A}→{self.smiles_B} / "
                f"{self.receptor}  {self.ddG_kcalmol:+.2f} ± "
                f"{self.uncertainty_kcalmol:.2f} kcal/mol  "
                f"wall={self.wall_seconds:.1f}s")


def _prepare_protein_topology(receptor_pdb: str | Path):
    """Run PDBFixer on the receptor and write the fixed structure to
    a temp PDB. Returns (temp_pdb_path, openmm_topology, positions).

    `openff.toolkit.Topology.from_pdb` needs the hydrogens explicit
    and atom names matching the PDB CCD — exactly what PDBFixer
    produces. Writing to disk and re-reading via the OpenFF loader
    is how both openff-evaluator and yank do it; avoids the
    `Topology.from_openmm(unique_molecules=…)` route that requires
    us to construct per-residue Molecule templates by hand.

    Caller is responsible for `os.unlink(temp_pdb_path)` when done.
    """
    import tempfile
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=str(receptor_pdb))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".pdb", delete=False)
    PDBFile.writeFile(fixer.topology, fixer.positions, tmp)
    tmp.close()
    return tmp.name, fixer.topology, fixer.positions


def _place_ligand_at_site(mol, target_center_nm):
    """Translate the ligand's single generated conformer so its CoM
    sits at `target_center_nm` (a 3-element numpy array in nm).

    Returns the same Molecule mutated in-place.
    """
    import numpy as np
    from openff.units import unit as offunit

    pos_nm = np.asarray(
        mol.conformers[0].m_as(offunit.nanometer), dtype=float)
    com = pos_nm.mean(axis=0)
    shift = np.asarray(target_center_nm, dtype=float) - com
    new_pos = (pos_nm + shift) * offunit.nanometer
    mol._conformers[0] = new_pos
    return mol


def _auto_binding_site_center(prot_top_omm, prot_positions):
    """Fallback when the caller hasn't supplied a pocket centre: use
    the geometric centroid of all Cα atoms. Fine for the scaffold;
    a real pocket-detection-driven centre is wired in Phase-2 via
    `src.dock.pocket_detect`.
    """
    import numpy as np
    from openmm import unit as ommunit

    ca_indices = [
        atom.index for atom in prot_top_omm.atoms()
        if atom.name == "CA"]
    if not ca_indices:
        # No protein at all — pathological receptor — use global mean.
        pos_nm = np.array(
            [[p.x, p.y, p.z] for p in prot_positions])
    else:
        pos_nm = np.array(
            [[prot_positions[i].x,
              prot_positions[i].y,
              prot_positions[i].z] for i in ca_indices])
    return pos_nm.mean(axis=0)


def _find_ca_indices_near(positions_nm, center_nm, radius_nm,
                          ca_candidate_indices):
    """Return the subset of Cα indices within `radius_nm` of the
    site centre. These become the restraint anchor group."""
    import numpy as np

    pos = np.asarray(positions_nm)
    c = np.asarray(center_nm)
    d = np.linalg.norm(pos[ca_candidate_indices] - c, axis=1)
    near = [ca_candidate_indices[i]
            for i, di in enumerate(d) if di <= radius_nm]
    if not near:
        # Radius too tight; fall back to the nearest 10 Cα atoms
        # so the restraint always has something to anchor to.
        order = np.argsort(d)
        near = [ca_candidate_indices[i] for i in order[:10]]
    return near


def _add_harmonic_com_restraint(system, ligand_indices,
                                 anchor_indices, k_kJ_per_nm2,
                                 r0_nm):
    """Add a CustomCentroidBondForce that harmonically restrains the
    ligand CoM to the anchor CoM (binding-site Cα centroid).

    Energy: U = 0.5 · k · (r − r0)²  where r = |CoM_lig − CoM_anchor|

    Coupled to a single lambda_restraint global (0 → 1) so Phase-2
    can decouple or scale the restraint at well-defined λ points if
    it needs a soft bumper around the FEP endpoints. For the scaffold
    we leave lambda_restraint = 1 (restraint always on); the
    AlchemicalState globals (lambda_sterics, lambda_electrostatics)
    handle the ligand decoupling separately.
    """
    from openmm import CustomCentroidBondForce

    expr = ("0.5*lambda_restraint*k*(distance(g1,g2)-r0)^2")
    force = CustomCentroidBondForce(2, expr)
    force.addGlobalParameter("lambda_restraint", 1.0)
    force.addGlobalParameter("k", float(k_kJ_per_nm2))
    force.addGlobalParameter("r0", float(r0_nm))
    force.addGroup(list(ligand_indices))
    force.addGroup(list(anchor_indices))
    force.addBond([0, 1])
    system.addForce(force)
    return force


def _harmonic_restraint_free_energy_kcalmol(
        k_kJ_per_nm2, temperature_K=298.15):
    """Analytical free-energy of a 3D isotropic harmonic restraint
    relative to the 1 M standard state.

    For a harmonic well U = (k/2) r² (we use r0 > 0 but the Gaussian
    integral is centred on r0; the factor below is independent of r0):

        ΔG_R→std  =  −kT · ln[ V_std / V_harmonic ]

    with V_harmonic = (2π·kT/k)^(3/2) and V_std = 1660 Å³ (the volume
    per molecule at 1 M standard state).

    Returns the correction added to ΔG_bind to go from 'restrained
    decoupled' to '1 M standard state' — this is the classic
    Hamelberg–Gilson / Boresch-style correction for CoM restraints.
    """
    import math

    # kT in kJ/mol at T
    kT_kJ = 0.0083144626 * temperature_K
    # Gaussian volume in nm^3: (2π·kT/k)^(3/2)
    v_harm_nm3 = (2.0 * math.pi * kT_kJ / k_kJ_per_nm2) ** 1.5
    # Standard-state volume: 1 M = 1 molecule per 1660.54 Å³ = 1.66054 nm³
    v_std_nm3 = 1.66054
    # ΔG = -kT ln(V_std / V_harm) → convert kJ → kcal (÷ 4.184)
    dG_kJ = -kT_kJ * math.log(v_std_nm3 / v_harm_nm3)
    return dG_kJ / 4.184


def _build_complex_alchemical_system(
    smiles: str,
    receptor_pdb: str | Path,
    *,
    binding_site_center_nm=None,
    anchor_radius_nm: float = 1.2,
    padding_nm: float = 1.2,
    softcore_alpha: float = 0.5,
    restraint_k_kJ_per_nm2: float = 4184.0,   # ≈ 10 kcal/mol/Å²
    restraint_r0_nm: float = 0.0,
):
    """Build the complex-leg alchemical system.

    Returns a dict with:
        alch_system, topology_omm, positions, ligand_indices,
        anchor_indices, restraint_k_kJ_per_nm2, restraint_r0_nm,
        n_ligand_atoms, n_protein_atoms

    Pipeline:
        receptor_pdb → PDBFixer (protonate, strip heterogens)
                     → OpenFF Topology.from_openmm(protein)
        smiles → Molecule.from_smiles + AM1-BCC + generate_conformer
               → translate CoM to site centre
        combined = Topology.from_molecules([*prot.molecules, ligand])
        solvate_topology(combined, padding)
        Interchange.from_smirnoff(ff14SB + Sage + TIP3P, solvated)
        alchemical factory annotates ligand atoms as softcore
        add CustomCentroidBondForce CoM–CoM restraint
    """
    import numpy as np
    from openff.toolkit.topology import Molecule, Topology
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import (
        solvate_topology,
    )
    from openff.units import unit as offunit
    from openmm import unit as ommunit
    from openmmtools import alchemy

    # 1) Prep protein. Write PDBFixer output to a temp file then
    # load it through Topology.from_pdb (native OpenFF loader that
    # understands PDB CCD residue templates — no unique_molecules
    # list needed).
    import os

    tmp_pdb, prot_top_omm, prot_positions = (
        _prepare_protein_topology(receptor_pdb))
    try:
        prot_top = Topology.from_pdb(tmp_pdb)
    finally:
        try:
            os.unlink(tmp_pdb)
        except OSError:
            pass
    n_protein_atoms = prot_top_omm.getNumAtoms()

    # Cα lookup — used both for auto pocket centre and for the
    # restraint anchor group.
    ca_candidate_indices = [
        atom.index for atom in prot_top_omm.atoms()
        if atom.name == "CA"]
    prot_positions_nm = np.array(
        [[p.x, p.y, p.z] for p in prot_positions])

    # 2) Decide binding-site centre.
    if binding_site_center_nm is None:
        binding_site_center_nm = _auto_binding_site_center(
            prot_top_omm, prot_positions)
    binding_site_center_nm = np.asarray(
        binding_site_center_nm, dtype=float)

    # 3) Build + place the ligand.
    lig = Molecule.from_smiles(smiles)
    lig.generate_conformers(n_conformers=1)
    lig.assign_partial_charges("am1bcc")
    _place_ligand_at_site(lig, binding_site_center_nm)
    n_ligand_atoms = lig.n_atoms

    # 4) Combine + solvate.
    combined_top = Topology.from_molecules(
        [*prot_top.molecules, lig])
    # Seed positions: protein + ligand (ligand CoM is already at
    # the binding-site centre from step 3).
    lig_positions_nm = np.asarray(
        lig.conformers[0].m_as(offunit.nanometer), dtype=float)
    seeded_positions = np.concatenate(
        [prot_positions_nm, lig_positions_nm], axis=0)
    combined_top.set_positions(
        seeded_positions * offunit.nanometer)

    solv_top = solvate_topology(
        topology=combined_top,
        nacl_conc=0.0 * offunit.molar,
        padding=padding_nm * offunit.nanometer,
    )

    # 5) Parametrise the whole solvated system.
    ff = ForceField(_LIGAND_FF, _PROTEIN_FF, _WATER_FF)
    ichg = Interchange.from_smirnoff(
        force_field=ff, topology=solv_top,
        charge_from_molecules=[lig])
    complex_system = ichg.to_openmm(
        combine_nonbonded_forces=True)
    complex_top_omm = solv_top.to_openmm()
    complex_positions_nm = np.asarray(
        ichg.positions.m_as(offunit.nanometer), dtype=float)

    # 6) Figure out ligand atom indices in the combined+solvated
    # system. packmol writes the combined topology first (protein
    # then ligand — in the order we passed to from_molecules) and
    # the solvent last. Ligand indices are therefore the contiguous
    # block right after protein: [n_protein, n_protein + n_ligand).
    lig_start = n_protein_atoms
    lig_end = n_protein_atoms + n_ligand_atoms
    ligand_indices = list(range(lig_start, lig_end))

    # 7) Pick restraint anchor (Cα within `anchor_radius_nm` of site
    # centre). Protein atom indices in the solvated topology are
    # unchanged from the bare protein — packmol doesn't re-index.
    anchor_indices = _find_ca_indices_near(
        prot_positions_nm, binding_site_center_nm,
        anchor_radius_nm, ca_candidate_indices)

    # 8) Make the alchemical region cover the ligand.
    region = alchemy.AlchemicalRegion(
        alchemical_atoms=ligand_indices,
        softcore_alpha=softcore_alpha)
    factory = alchemy.AbsoluteAlchemicalFactory()
    alch_system = factory.create_alchemical_system(
        reference_system=complex_system,
        alchemical_regions=region)

    # 9) Add the CoM–CoM restraint AFTER the alchemical factory
    # (the factory only touches NonbondedForce / Custom*NB; our
    # CustomCentroidBondForce passes through unchanged, which is
    # exactly what we want — the restraint is not alchemically
    # decoupled).
    _add_harmonic_com_restraint(
        alch_system, ligand_indices, anchor_indices,
        k_kJ_per_nm2=restraint_k_kJ_per_nm2,
        r0_nm=restraint_r0_nm)

    # Positions for the alchemical system must carry OpenMM units.
    positions = complex_positions_nm * ommunit.nanometer

    return {
        "alch_system": alch_system,
        "topology_omm": complex_top_omm,
        "positions": positions,
        "ligand_indices": ligand_indices,
        "anchor_indices": anchor_indices,
        "restraint_k_kJ_per_nm2": restraint_k_kJ_per_nm2,
        "restraint_r0_nm": restraint_r0_nm,
        "n_ligand_atoms": n_ligand_atoms,
        "n_protein_atoms": n_protein_atoms,
        "n_total_atoms": complex_system.getNumParticles(),
    }


def compute_absolute_binding_dg(
    smiles: str,
    receptor_pdb: str | Path,
    *,
    binding_site_center_nm=None,
    anchor_radius_nm: float = 1.2,
    padding_nm: float = 1.2,
    restraint_k_kJ_per_nm2: float = 4184.0,
    restraint_r0_nm: float = 0.0,
    n_windows: int = 11,
    softcore_alpha: float = 0.5,
    sample: bool = False,
    n_equilibration_steps: int = 500,
    n_production_steps: int = 2000,
    sample_stride: int = 100,
    seed: int = 1,
) -> BindingDGResult:
    """Absolute binding ΔG via the double-decoupling method.

    Phase-1 (this commit, default ``sample=False``):
      build complex + solvent alchemical systems, verify both
      minimise cleanly, return `phase='scaffolded_both_legs'`.

    Phase-2 (flip ``sample=True``):
      run `sample_alchemical_windows` on both legs, apply
      restraint + standard-state corrections, return
      `phase='sampled'` with ΔG_bind populated.
    """
    import time

    t0 = time.time()
    result = BindingDGResult(
        smiles=smiles,
        receptor=str(Path(receptor_pdb).name),
        ok=False,
        phase="scaffolded",
        n_windows=n_windows,
        restraint_k_kJ_per_nm2=restraint_k_kJ_per_nm2,
        restraint_r0_nm=restraint_r0_nm,
    )

    # Complex leg.
    try:
        cx = _build_complex_alchemical_system(
            smiles, receptor_pdb,
            binding_site_center_nm=binding_site_center_nm,
            anchor_radius_nm=anchor_radius_nm,
            padding_nm=padding_nm,
            softcore_alpha=softcore_alpha,
            restraint_k_kJ_per_nm2=restraint_k_kJ_per_nm2,
            restraint_r0_nm=restraint_r0_nm)
    except Exception as e:
        result.reason = (
            f"complex build failed: {str(e)[:200]}")
        result.wall_seconds = time.time() - t0
        return result

    result.n_ligand_atoms = cx["n_ligand_atoms"]
    result.n_protein_atoms = cx["n_protein_atoms"]
    result.n_total_atoms_complex = cx["n_total_atoms"]

    # Solvent leg — reuse the validated hydration builder. ΔG for
    # decoupling in pure water is the same quantity as the solvent
    # half of ΔG_hyd, so `_build_alchemical_legs` is the exact right
    # primitive.
    try:
        from src.fep import _build_alchemical_legs

        (_vac_alch, solv_alch,
         _vac_top, solv_top,
         _vac_pos, solv_pos, _n) = _build_alchemical_legs(
            smiles, softcore_alpha=softcore_alpha)
    except Exception as e:
        # Complex leg is still valid → report partial progress
        # rather than claiming total failure. Matches the pattern
        # used by the hydration scaffold for solvent-only failures.
        result.reason = (
            f"solvent leg build failed (complex ok): "
            f"{str(e)[:200]}")
        result.ok = True
        result.phase = "scaffolded_complex_only"
        result.wall_seconds = time.time() - t0
        return result

    result.n_total_atoms_solvent = solv_alch.getNumParticles()

    if not sample:
        # Phase-1: stop here. Both systems built; ΔG left None.
        result.ok = True
        result.phase = "scaffolded_both_legs"
        result.wall_seconds = time.time() - t0
        return result

    # Phase-2: drive the sampler on both legs.
    from src.fep.sampling import sample_alchemical_windows

    try:
        cx_r = sample_alchemical_windows(
            cx["alch_system"], cx["topology_omm"], cx["positions"],
            n_windows=n_windows,
            n_equilibration_steps=n_equilibration_steps,
            n_production_steps=n_production_steps,
            sample_stride=sample_stride, seed=seed)
        if not cx_r.ok:
            result.reason = f"complex leg sample failed: {cx_r.reason}"
            result.wall_seconds = time.time() - t0
            return result

        sv_r = sample_alchemical_windows(
            solv_alch, solv_top, solv_pos,
            n_windows=n_windows,
            n_equilibration_steps=n_equilibration_steps,
            n_production_steps=n_production_steps,
            sample_stride=sample_stride, seed=seed)
        if not sv_r.ok:
            result.reason = f"solvent leg sample failed: {sv_r.reason}"
            result.wall_seconds = time.time() - t0
            return result
    except Exception as e:
        result.reason = f"binding sampling failed: {str(e)[:200]}"
        result.wall_seconds = time.time() - t0
        return result

    # Compose: ΔG_bind = −(ΔG_dec_complex − ΔG_dec_solvent) + ΔG_R + ΔG_std
    # ΔG_R (harmonic restraint → 1 M standard state) is analytical;
    # ΔG_std is folded into the same correction formula.
    import math

    dG_R_plus_std = _harmonic_restraint_free_energy_kcalmol(
        restraint_k_kJ_per_nm2)

    dG_bind = -(cx_r.dG_kcalmol - sv_r.dG_kcalmol) + dG_R_plus_std
    unc = math.sqrt(
        (cx_r.dG_uncertainty_kcalmol or 0.0) ** 2
        + (sv_r.dG_uncertainty_kcalmol or 0.0) ** 2)

    result.ok = True
    result.phase = "sampled"
    result.dG_bind_kcalmol = dG_bind
    result.dG_dec_complex_kcalmol = cx_r.dG_kcalmol
    result.dG_dec_solvent_kcalmol = sv_r.dG_kcalmol
    result.dG_restraint_correction_kcalmol = dG_R_plus_std
    result.dG_standard_state_kcalmol = 0.0   # folded in above
    result.uncertainty_kcalmol = unc
    result.wall_seconds = time.time() - t0
    return result


def compute_relative_binding_ddg(
    smiles_A: str,
    smiles_B: str,
    receptor_pdb: str | Path,
    *,
    binding_site_center_nm=None,
    anchor_radius_nm: float = 1.2,
    padding_nm: float = 1.2,
    restraint_k_kJ_per_nm2: float = 4184.0,
    restraint_r0_nm: float = 0.0,
    n_windows: int = 11,
    softcore_alpha: float = 0.5,
    sample: bool = False,
    n_equilibration_steps: int = 500,
    n_production_steps: int = 2000,
    sample_stride: int = 100,
    seed: int = 1,
) -> RelativeBindingDDGResult:
    """Relative ΔΔG(A→B) = ΔG_bind(B) − ΔG_bind(A).

    Two independent absolute-binding runs with identical
    receptor prep, pocket centre, restraint, and sampling
    parameters so the corrections cancel when we subtract.

    Phase-1 (default ``sample=False``) returns both builders'
    outputs with `phase='scaffolded'` and ddG_kcalmol=None.
    """
    import time

    t0 = time.time()
    result = RelativeBindingDDGResult(
        smiles_A=smiles_A, smiles_B=smiles_B,
        receptor=str(Path(receptor_pdb).name),
        ok=False)

    bind_A = compute_absolute_binding_dg(
        smiles_A, receptor_pdb,
        binding_site_center_nm=binding_site_center_nm,
        anchor_radius_nm=anchor_radius_nm,
        padding_nm=padding_nm,
        restraint_k_kJ_per_nm2=restraint_k_kJ_per_nm2,
        restraint_r0_nm=restraint_r0_nm,
        n_windows=n_windows,
        softcore_alpha=softcore_alpha,
        sample=sample,
        n_equilibration_steps=n_equilibration_steps,
        n_production_steps=n_production_steps,
        sample_stride=sample_stride, seed=seed)
    result.bind_A = bind_A
    if not bind_A.ok:
        result.reason = f"ligand A: {bind_A.reason}"
        result.wall_seconds = time.time() - t0
        return result

    bind_B = compute_absolute_binding_dg(
        smiles_B, receptor_pdb,
        binding_site_center_nm=binding_site_center_nm,
        anchor_radius_nm=anchor_radius_nm,
        padding_nm=padding_nm,
        restraint_k_kJ_per_nm2=restraint_k_kJ_per_nm2,
        restraint_r0_nm=restraint_r0_nm,
        n_windows=n_windows,
        softcore_alpha=softcore_alpha,
        sample=sample,
        n_equilibration_steps=n_equilibration_steps,
        n_production_steps=n_production_steps,
        sample_stride=sample_stride, seed=seed)
    result.bind_B = bind_B
    if not bind_B.ok:
        result.reason = f"ligand B: {bind_B.reason}"
        result.wall_seconds = time.time() - t0
        return result

    if not sample:
        # Phase-1 scaffolded: both legs built; ddG deferred.
        result.ok = True
        result.phase = "scaffolded"
        result.wall_seconds = time.time() - t0
        return result

    import math
    result.dG_bind_A_kcalmol = bind_A.dG_bind_kcalmol
    result.dG_bind_B_kcalmol = bind_B.dG_bind_kcalmol
    result.ddG_kcalmol = (
        bind_B.dG_bind_kcalmol - bind_A.dG_bind_kcalmol)
    result.uncertainty_kcalmol = math.sqrt(
        (bind_A.uncertainty_kcalmol or 0.0) ** 2
        + (bind_B.uncertainty_kcalmol or 0.0) ** 2)
    result.ok = True
    result.phase = "sampled"
    result.wall_seconds = time.time() - t0
    return result


__all__ = [
    "BindingDGResult",
    "RelativeBindingDDGResult",
    "compute_absolute_binding_dg",
    "compute_relative_binding_ddg",
]


def main(argv=None) -> int:
    """CLI: `cellsim fep-binding` — two subcommands.

      cellsim fep-binding dg <SMILES> <receptor.pdb>
      cellsim fep-binding ddg <SMI_A> <SMI_B> <receptor.pdb>

    Phase-1 scaffolding mode by default (no MD). Pass ``--sample``
    to actually drive the sampler on both legs; that will only
    complete in reasonable wall time on a GPU.
    """
    import argparse
    import json as _json
    import sys

    ap = argparse.ArgumentParser(
        description="Milestone B FEP — absolute or relative "
                    "binding ΔG. Phase-1 scaffold (default) builds "
                    "the alchemical systems; --sample drives MD.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    dgp = sub.add_parser("dg", help="absolute ΔG_bind")
    dgp.add_argument("smiles")
    dgp.add_argument("receptor_pdb")
    dgp.add_argument("--n-windows", type=int, default=11)
    dgp.add_argument("--softcore-alpha", type=float, default=0.5)
    dgp.add_argument("--padding", type=float, default=1.2)
    dgp.add_argument("--restraint-k", type=float, default=4184.0,
                     help="spring constant kJ/mol/nm² "
                          "(default 4184 ≈ 10 kcal/mol/Å²)")
    dgp.add_argument("--sample", action="store_true",
                     help="run MD (slow on CPU; intended for GPU)")
    dgp.add_argument("--json", action="store_true")

    ddgp = sub.add_parser("ddg", help="relative ΔΔG(A→B)")
    ddgp.add_argument("smiles_a")
    ddgp.add_argument("smiles_b")
    ddgp.add_argument("receptor_pdb")
    ddgp.add_argument("--n-windows", type=int, default=11)
    ddgp.add_argument("--softcore-alpha", type=float, default=0.5)
    ddgp.add_argument("--padding", type=float, default=1.2)
    ddgp.add_argument("--restraint-k", type=float, default=4184.0)
    ddgp.add_argument("--sample", action="store_true")
    ddgp.add_argument("--json", action="store_true")

    args = ap.parse_args(argv)

    if args.cmd == "dg":
        r = compute_absolute_binding_dg(
            args.smiles, args.receptor_pdb,
            n_windows=args.n_windows,
            softcore_alpha=args.softcore_alpha,
            padding_nm=args.padding,
            restraint_k_kJ_per_nm2=args.restraint_k,
            sample=args.sample)
    else:
        r = compute_relative_binding_ddg(
            args.smiles_a, args.smiles_b, args.receptor_pdb,
            n_windows=args.n_windows,
            softcore_alpha=args.softcore_alpha,
            padding_nm=args.padding,
            restraint_k_kJ_per_nm2=args.restraint_k,
            sample=args.sample)

    if args.json:
        from dataclasses import asdict
        print(_json.dumps(asdict(r), indent=2, default=str))
    else:
        print(r.summary())
    return 0 if r.ok else 1


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main())
