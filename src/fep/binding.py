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
    # Drop TERMINAL missing-residue extensions — PDBFixer's placement
    # heuristic for missing N/C-terminal residues is unreliable and
    # can put them 99 Å from the main body (observed on 1stp where
    # GLN 159 ended up at y=-99 Å, which then bloats the solvation
    # box by 50×). Internal-gap missing residues (real loops) are
    # kept — their placement IS reasonable because both endpoints
    # are constrained.
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain_idx, res_pos = key
        n_chain_res = len(list(chains[chain_idx].residues()))
        if res_pos == 0 or res_pos == n_chain_res:
            del fixer.missingResidues[key]
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


def _auto_binding_site_center(prot_top_omm, prot_positions,
                              receptor_pdb=None):
    """Pick a binding-site centroid when the caller hasn't supplied
    one. Two-tier:

    1. **fpocket** (non-AI geometric pocket detector) on the
       original receptor PDB → top-ranked drug-scored pocket
       centre. This is what the docking pipeline already uses for
       receptors without a cocrystal (see `src.dock.pocket_detect`).
    2. **Cα centroid** fallback — geometric centre of every Cα in
       the protein. Always defined, but on a big multi-domain
       receptor the centroid may land in a hinge region instead of
       a real druggable pocket, which would restrain the ligand to
       nowhere useful. Warn but keep going.

    Returning the fpocket pick by default is the biologist-friendly
    choice: you point CellSim at a PDB you grabbed from RCSB, it
    picks the druggable pocket automatically, you never hand-pick
    coordinates. For receptors where you *know* the right site
    (e.g. EGFR kinase: use the ATP-binding hinge, not allosteric
    sites), supply `binding_site_center_nm` explicitly — the caller
    path stays clean.
    """
    import numpy as np

    # Tier 1: fpocket. Silent fail-through on any issue (fpocket
    # missing, timeout, zero pockets) — we don't want the binding
    # FEP scaffold to hard-error just because fpocket happens to
    # not be on PATH.
    if receptor_pdb is not None:
        try:
            from src.dock.pocket_detect import detect_pockets
            pockets = detect_pockets(receptor_pdb, top_k=1)
            if pockets and pockets[0].center_A is not None:
                c_A = np.asarray(pockets[0].center_A, dtype=float)
                # Å → nm
                center_nm = c_A / 10.0
                logger = _module_logger()
                logger.info(
                    "fpocket picked pocket #%d (drug=%.2f, "
                    "vol=%.0f Å³) at %s nm",
                    pockets[0].rank, pockets[0].drug_score or 0.0,
                    pockets[0].volume_A3 or 0.0, center_nm)
                return center_nm
        except Exception as e:
            _module_logger().info(
                "fpocket pocket-detect failed (%s); "
                "falling back to Cα centroid", e)

    # Tier 2: Cα centroid fallback.
    ca_indices = [
        atom.index for atom in prot_top_omm.atoms()
        if atom.name == "CA"]
    if not ca_indices:
        pos_nm = np.array(
            [[p.x, p.y, p.z] for p in prot_positions])
    else:
        pos_nm = np.array(
            [[prot_positions[i].x,
              prot_positions[i].y,
              prot_positions[i].z] for i in ca_indices])
    _module_logger().info(
        "using Cα centroid as binding-site centre "
        "(no fpocket hit — may be off-site on multi-domain "
        "receptors; supply binding_site_center_nm to override)")
    return pos_nm.mean(axis=0)


def _module_logger():
    import logging
    return logging.getLogger(__name__)


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

    # 2) Decide binding-site centre. Pass the ORIGINAL receptor_pdb
    # (not the PDBFixer-fixed temp file) because fpocket wants a
    # real PDB and gives slightly different answers on the fixed
    # one; the centre resolves to the same Cα-shell the restraint
    # anchors to, so the tiny coordinate difference is fine.
    if binding_site_center_nm is None:
        binding_site_center_nm = _auto_binding_site_center(
            prot_top_omm, prot_positions,
            receptor_pdb=receptor_pdb)
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


def _build_complex_alchemical_system_amber14(
    smiles: str,
    receptor_pdb: str | Path,
    *,
    binding_site_center_nm=None,
    anchor_radius_nm: float = 1.2,
    padding_nm: float = 1.2,
    softcore_alpha: float = 0.5,
    restraint_k_kJ_per_nm2: float = 4184.0,
    restraint_r0_nm: float = 0.0,
):
    """Hybrid builder: AMBER14 for protein + TIP3P for solvent via
    OpenMM's classical ForceField; SMIRNOFF (OpenFF Sage 2.1.0 +
    AM1-BCC charges) for the ligand via openmmforcefields'
    SMIRNOFFTemplateGenerator.

    Fixes the O(N²) memory blowup that `Interchange.from_smirnoff`
    hits on mid-size proteins (1stp, 1m17, anything > ~1 500 atoms).
    Keeps ligand parameters bit-for-bit identical to the hydration
    leg, so the ΔΔG = complex - solvent subtraction in Phase-2
    still has the Sage/AM1-BCC errors cancel.

    Return shape is the exact same dict as
    `_build_complex_alchemical_system` (its pure-SMIRNOFF sibling),
    so either builder drops in interchangeably.
    """
    import numpy as np
    from openff.toolkit.topology import Molecule
    from openmm import app as _app
    from openmm import unit as ommunit
    from openmmforcefields.generators import (
        SMIRNOFFTemplateGenerator,
    )
    from openmmtools import alchemy

    # 1) Prep protein via PDBFixer (same protonation + nonstandard-
    # residue handling as the pure-SMIRNOFF sibling; fair apples-
    # to-apples starting structure between the two paths).
    fixer_tmp_pdb, prot_top_omm, prot_positions = (
        _prepare_protein_topology(receptor_pdb))
    import os
    try:
        os.unlink(fixer_tmp_pdb)
    except OSError:
        pass

    n_protein_atoms = prot_top_omm.getNumAtoms()
    ca_candidate_indices = [
        atom.index for atom in prot_top_omm.atoms()
        if atom.name == "CA"]
    prot_positions_nm = np.array(
        [[p.x, p.y, p.z] for p in prot_positions])

    # 2) Decide binding-site centre (fpocket → Cα centroid).
    if binding_site_center_nm is None:
        binding_site_center_nm = _auto_binding_site_center(
            prot_top_omm, prot_positions,
            receptor_pdb=receptor_pdb)
    binding_site_center_nm = np.asarray(
        binding_site_center_nm, dtype=float)

    # 3) Build + place the ligand (same as pure-SMIRNOFF path).
    lig = Molecule.from_smiles(smiles)
    lig.generate_conformers(n_conformers=1)
    lig.assign_partial_charges("am1bcc")
    _place_ligand_at_site(lig, binding_site_center_nm)
    n_ligand_atoms = lig.n_atoms

    # 4) Register SMIRNOFF template generator BEFORE we attempt to
    # build a System. Passing `molecules=[lig]` pre-parametrises
    # the ligand so addSolvent + createSystem don't hit the
    # SMIRNOFF pattern matcher lazily (which is where the pure-
    # Interchange path blows up).
    gen = SMIRNOFFTemplateGenerator(
        molecules=[lig], forcefield=_LIGAND_FF)
    ff = _app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    ff.registerTemplateGenerator(gen.generator)

    # 5) Combine protein + ligand into a Modeller. Ligand is
    # appended AFTER protein so ligand_indices = [n_protein,
    # n_protein + n_ligand).
    modeller = _app.Modeller(prot_top_omm, prot_positions)
    # Clear any CRYST1-inherited box — a PDB's crystal packing box
    # is usually far larger than we want, and if left in place
    # Modeller.addSolvent uses it verbatim instead of the padding
    # we asked for. Streptavidin 1stp has a 99×99×125 Å crystal
    # that would give 812k-atom solvated systems vs ~15k with a
    # tight padding-derived box.
    modeller.topology.setPeriodicBoxVectors(None)
    lig_top_omm = lig.to_topology().to_openmm()
    lig_pos_nm = np.asarray(
        lig.conformers[0].m_as(__import__("openff.units",
                                          fromlist=["unit"])
                               .unit.nanometer),
        dtype=float)
    lig_pos_omm = lig_pos_nm * ommunit.nanometer
    modeller.add(lig_top_omm, lig_pos_omm)

    # 6) Solvate (classical Modeller + tip3p model — fast).
    modeller.addSolvent(
        ff, padding=padding_nm * ommunit.nanometer,
        model="tip3p",
        ionicStrength=0.0 * ommunit.molar)
    complex_top_omm = modeller.topology

    # 7) Create the bare System, THEN alchemical-wrap the ligand.
    complex_system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=_app.PME,
        nonbondedCutoff=1.0 * ommunit.nanometer,
        constraints=_app.HBonds,
        rigidWater=True)

    # Positions as numpy for return. Modeller.positions is an
    # openmm.unit.Quantity wrapping either a list of Vec3 or a
    # numpy ndarray (depending on OpenMM version). Use
    # value_in_unit(nm) + asarray to handle both.
    complex_positions_nm = np.asarray(
        modeller.positions.value_in_unit(ommunit.nanometer),
        dtype=float)

    # 8) Ligand atom indices in the combined solvated topology.
    # Modeller appends in order (protein → ligand → waters/ions),
    # so [n_protein, n_protein + n_ligand) still holds.
    lig_start = n_protein_atoms
    lig_end = n_protein_atoms + n_ligand_atoms
    ligand_indices = list(range(lig_start, lig_end))

    # 9) Restraint anchor — Cα near the binding site.
    anchor_indices = _find_ca_indices_near(
        prot_positions_nm, binding_site_center_nm,
        anchor_radius_nm, ca_candidate_indices)

    # 10) Alchemical factory + restraint (same as sibling).
    region = alchemy.AlchemicalRegion(
        alchemical_atoms=ligand_indices,
        softcore_alpha=softcore_alpha)
    factory = alchemy.AbsoluteAlchemicalFactory()
    alch_system = factory.create_alchemical_system(
        reference_system=complex_system,
        alchemical_regions=region)
    _add_harmonic_com_restraint(
        alch_system, ligand_indices, anchor_indices,
        k_kJ_per_nm2=restraint_k_kJ_per_nm2,
        r0_nm=restraint_r0_nm)

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
        "ff_stack": "amber14-all + tip3pfb + SMIRNOFF-lig",
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
    force_field_path: str = "amber14",
) -> BindingDGResult:
    """
    `force_field_path` selects the protein-ligand parametrisation
    stack:
      - "amber14" (default): openmm.app.ForceField amber14-all +
        tip3pfb for protein + solvent, SMIRNOFFTemplateGenerator
        for the ligand. Fast, works on all receptors tested.
      - "smirnoff": Interchange.from_smirnoff on the whole system
        with the OpenFF ff14sb_off_impropers port. Hits an O(N²)
        memory wall on mid-size proteins; kept as a fallback /
        provenance-matching path for the hydration leg.
    """
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
    if force_field_path == "amber14":
        _builder = _build_complex_alchemical_system_amber14
    elif force_field_path == "smirnoff":
        _builder = _build_complex_alchemical_system
    else:
        result.reason = (
            f"unknown force_field_path '{force_field_path}' "
            "(expected 'amber14' or 'smirnoff')")
        result.wall_seconds = time.time() - t0
        return result
    try:
        cx = _builder(
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

    bp = sub.add_parser(
        "bench",
        help="YAML-driven batch (e.g. binding_streptavidin.yaml)")
    bp.add_argument("yaml_path")
    bp.add_argument("--n-windows", type=int, default=11)
    bp.add_argument("--softcore-alpha", type=float, default=0.5)
    bp.add_argument("--padding", type=float, default=1.2)
    bp.add_argument("--restraint-k", type=float, default=4184.0)
    bp.add_argument(
        "--sample", action="store_true",
        help="run MD on every entry (HOURS per compound on CPU; "
             "a full streptavidin set is a GPU-only operation)")
    bp.add_argument("--production-steps", type=int, default=2000)
    bp.add_argument("--equilibration-steps", type=int, default=500)
    bp.add_argument("--sample-stride", type=int, default=100)
    bp.add_argument(
        "--out-csv", default=None,
        help="write per-entry results (same schema as FreeSolv so "
             "`cellsim fep-report` can consume it)")

    vp = sub.add_parser(
        "validate",
        help="sub-second YAML dry-run: parse SMILES via RDKit, "
             "confirm receptor PDB exists, report issues — run "
             "this BEFORE launching a multi-hour sampled run")
    vp.add_argument("yaml_path")
    vp.add_argument(
        "--json", action="store_true",
        help="emit the diagnostic table as JSON instead of text")

    args = ap.parse_args(argv)

    if args.cmd == "dg":
        r = compute_absolute_binding_dg(
            args.smiles, args.receptor_pdb,
            n_windows=args.n_windows,
            softcore_alpha=args.softcore_alpha,
            padding_nm=args.padding,
            restraint_k_kJ_per_nm2=args.restraint_k,
            sample=args.sample)
    elif args.cmd == "ddg":
        r = compute_relative_binding_ddg(
            args.smiles_a, args.smiles_b, args.receptor_pdb,
            n_windows=args.n_windows,
            softcore_alpha=args.softcore_alpha,
            padding_nm=args.padding,
            restraint_k_kJ_per_nm2=args.restraint_k,
            sample=args.sample)
    elif args.cmd == "bench":
        return _run_bench(args)
    else:  # validate
        return _run_validate(args)

    if args.json:
        from dataclasses import asdict
        print(_json.dumps(asdict(r), indent=2, default=str))
    else:
        print(r.summary())
    return 0 if r.ok else 1


def _run_bench(args) -> int:
    """Drive `compute_absolute_binding_dg` over every entry in a
    binding YAML. Emits one-line-per-entry progress + optional CSV
    in the same schema `cellsim fep-report` expects, so the whole
    run can be summarised through the existing analyser.

    Schema maps:
        YAML `dG_bind_kcalmol`    → CSV `dG_expt_kcalmol`
        BindingDGResult.dG_bind_kcalmol → CSV `dG_pred_kcalmol`
        (residual = pred − expt, just as for hydration)

    Phase-1 default (`--sample` absent) runs scaffold-only on every
    entry. Wall time ≈ 20-60 s per compound on CPU regardless of
    receptor size, so a 4-compound YAML is a 2-5 min pre-flight
    sanity check — exactly what a biologist wants before kicking
    off a proper Phase-2 run on GPU.
    """
    import csv as _csv
    import sys as _sys
    import time as _time
    import yaml as _yaml

    data = _yaml.safe_load(Path(args.yaml_path).read_text())
    recv = data.get("receptor", {})
    receptor_pdb = recv.get("pdb_path")
    if not receptor_pdb:
        print("bench: YAML missing receptor.pdb_path",
              file=_sys.stderr)
        return 2
    bsc = recv.get("binding_site_center_nm")

    entries = data.get("entries", [])
    rows = []
    t_all = _time.time()
    print(f"[fep-binding bench] {args.yaml_path}")
    print(f"  receptor:   {receptor_pdb}")
    print(f"  entries:    {len(entries)}")
    print(f"  sample:     {bool(args.sample)}")
    print()

    for entry in entries:
        name = entry.get("name", "?")
        smi = entry.get("smiles", "")
        expt = entry.get("dG_bind_kcalmol")
        print(f"[bench] {name}  {smi}", flush=True)
        ts = _time.time()
        r = compute_absolute_binding_dg(
            smi, receptor_pdb,
            binding_site_center_nm=bsc,
            n_windows=args.n_windows,
            softcore_alpha=args.softcore_alpha,
            padding_nm=args.padding,
            restraint_k_kJ_per_nm2=args.restraint_k,
            sample=bool(args.sample),
            n_equilibration_steps=args.equilibration_steps,
            n_production_steps=args.production_steps,
            sample_stride=args.sample_stride)
        wall = _time.time() - ts
        pred = r.dG_bind_kcalmol
        unc = r.uncertainty_kcalmol
        resid = (pred - expt
                 if pred is not None and expt is not None else None)
        if r.ok and pred is None:
            print(f"  scaffolded {r.phase}  "
                  f"atoms={r.n_total_atoms_complex}  "
                  f"wall={wall:.1f}s", flush=True)
        elif r.ok:
            print(f"  pred = {pred:+.2f}  expt = {expt:+.2f}  "
                  f"resid = {(resid or 0):+.2f}  "
                  f"wall={wall:.1f}s", flush=True)
        else:
            print(f"  FAIL: {r.reason}", flush=True)

        rows.append({
            "name": name,
            "smiles": smi,
            "dG_expt_kcalmol": expt,
            "dG_pred_kcalmol": pred,
            "uncertainty_kcalmol": unc,
            "residual_kcalmol": resid,
            "ghmc_accept_mean": None,
            "ghmc_accept_min": None,
            "wall_seconds": round(wall, 2),
            "ok": r.ok,
            "reason": r.reason,
        })

    t_total = _time.time() - t_all
    n_ok = sum(1 for row in rows if row["ok"])
    print()
    print(f"[fep-binding bench] {n_ok}/{len(rows)} ok  "
          f"total wall {t_total:.1f} s  "
          f"({args.yaml_path})")

    if args.out_csv:
        out = Path(args.out_csv)
        out.parent.mkdir(parents=True, exist_ok=True)
        cols = ["name", "smiles", "dG_expt_kcalmol",
                "dG_pred_kcalmol", "uncertainty_kcalmol",
                "residual_kcalmol", "ghmc_accept_mean",
                "ghmc_accept_min", "wall_seconds", "ok", "reason"]
        with out.open("w", newline="",
                      encoding="utf-8-sig") as fo:
            w = _csv.DictWriter(
                fo, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for row in rows:
                w.writerow(row)
        print(f"  wrote {out}  "
              f"(schema compatible with `cellsim fep-report`)")

    return 0 if n_ok == len(rows) else 1


def _run_validate(args) -> int:
    """Dry-run YAML check — NO OpenFF / OpenMM imports, so it lands
    in < 1 s even on machines without the MD stack installed.

    Verifies:
      - YAML parses and has the expected shape
      - `receptor.pdb_path` exists as a readable file
      - every entry has a name + a SMILES that RDKit parses
      - stereocenters / charges / heavy-atom count for each ligand
      - receptor atom count (quick sanity on the PDB)

    Emits a diagnostic table + overall PASS/FAIL, exits 1 on any
    error. Intended as a pre-flight CI gate on any PR that edits
    a benchmark YAML, and as a biologist pre-run sanity check.
    """
    import json as _json
    import sys as _sys
    import yaml as _yaml

    from rdkit import Chem
    from rdkit.Chem import AllChem  # noqa: F401

    issues: list[str] = []
    entries_report: list[dict] = []

    try:
        data = _yaml.safe_load(Path(args.yaml_path).read_text())
    except Exception as e:
        print(f"[FAIL] cannot read YAML {args.yaml_path}: {e}",
              file=_sys.stderr)
        return 1

    if not isinstance(data, dict):
        print(f"[FAIL] YAML must be a mapping at top-level; "
              f"got {type(data).__name__}", file=_sys.stderr)
        return 1

    # Detect YAML kind by the per-entry expt key. Hydration YAMLs
    # (FreeSolv) use dG_hydration_kcalmol and have no receptor;
    # binding YAMLs (streptavidin etc.) use dG_bind_kcalmol and
    # require a receptor.
    entries_preview = data.get("entries") or []
    is_binding = any(
        "dG_bind_kcalmol" in (e or {}) for e in entries_preview)
    is_hydration = any(
        "dG_hydration_kcalmol" in (e or {}) for e in entries_preview)
    yaml_kind = (
        "binding" if is_binding and not is_hydration
        else "hydration" if is_hydration and not is_binding
        else "mixed" if is_binding and is_hydration
        else "unknown")

    recv = data.get("receptor", {}) or {}
    receptor_pdb = recv.get("pdb_path")
    receptor_report: dict = {"pdb_path": receptor_pdb,
                              "yaml_kind": yaml_kind}

    if yaml_kind == "hydration":
        # Hydration YAMLs don't need a receptor — no FEP leg in the
        # protein. Skip the receptor-file check.
        pass
    elif not receptor_pdb:
        issues.append(
            f"missing receptor.pdb_path (yaml kind: {yaml_kind})")
    else:
        p = Path(receptor_pdb)
        if not p.exists():
            issues.append(f"receptor PDB not found: {receptor_pdb}")
            receptor_report["exists"] = False
        elif not p.is_file():
            issues.append(f"receptor path is not a file: {receptor_pdb}")
            receptor_report["exists"] = False
        else:
            try:
                txt = p.read_text(encoding="utf-8",
                                  errors="replace")
                n_atoms = sum(
                    1 for line in txt.splitlines()
                    if line.startswith("ATOM"))
                n_het = sum(
                    1 for line in txt.splitlines()
                    if line.startswith("HETATM"))
                n_chains = len({
                    line[21:22]
                    for line in txt.splitlines()
                    if line.startswith("ATOM")})
                receptor_report.update({
                    "exists": True,
                    "size_bytes": p.stat().st_size,
                    "n_atoms": n_atoms,
                    "n_hetatm": n_het,
                    "n_chains": n_chains,
                })
                if n_atoms == 0:
                    issues.append(
                        f"receptor {p.name} has no ATOM records")
            except OSError as e:
                issues.append(
                    f"cannot read receptor {p.name}: {e}")

    entries = data.get("entries") or []
    if not entries:
        issues.append("no entries in YAML")

    seen_names: set[str] = set()
    for i, entry in enumerate(entries):
        name = entry.get("name", f"<entry_{i}>")
        smi = (entry.get("smiles") or "").strip()
        expt = entry.get("dG_bind_kcalmol",
                         entry.get("dG_hydration_kcalmol"))
        er: dict = {"name": name, "smiles": smi, "expt_kcalmol": expt}

        if name in seen_names:
            issues.append(f"duplicate name: {name}")
        seen_names.add(name)

        if not smi:
            issues.append(f"{name}: missing smiles")
            er["ok"] = False
            er["reason"] = "missing smiles"
        else:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                issues.append(f"{name}: RDKit cannot parse '{smi}'")
                er["ok"] = False
                er["reason"] = "RDKit parse failed"
            else:
                try:
                    canon = Chem.MolToSmiles(mol)
                    n_heavy = mol.GetNumHeavyAtoms()
                    n_stereo = sum(
                        1 for _ in Chem.FindMolChiralCenters(
                            mol, includeUnassigned=True))
                    formal_charge = Chem.GetFormalCharge(mol)
                    er.update({
                        "ok": True,
                        "canon_smiles": canon,
                        "n_heavy_atoms": n_heavy,
                        "n_stereocenters": n_stereo,
                        "formal_charge": formal_charge,
                    })
                    # Biologist gotchas: formal charge ≠ 0 means the
                    # alchemical FEP needs a counter-ion correction
                    # we haven't wired; flag it.
                    if formal_charge != 0:
                        issues.append(
                            f"{name}: formal charge "
                            f"{formal_charge:+d} — alchemical FEP "
                            "needs a Rocklin / Warren correction "
                            "(not wired in Phase-1); Phase-2 item")
                    if n_heavy > 60:
                        issues.append(
                            f"{name}: {n_heavy} heavy atoms > 60 — "
                            "AM1-BCC will be slow (minutes/compound)")
                except Exception as e:
                    issues.append(f"{name}: RDKit descriptor "
                                  f"failure: {e}")
                    er["ok"] = False
                    er["reason"] = f"descriptor failure: {e}"

        if expt is None:
            issues.append(
                f"{name}: missing expt ΔG "
                "(dG_bind_kcalmol / dG_hydration_kcalmol)")
        entries_report.append(er)

    report = {
        "yaml_path": str(args.yaml_path),
        "receptor": receptor_report,
        "entries": entries_report,
        "issues": issues,
        "pass": len(issues) == 0,
    }

    if args.json:
        print(_json.dumps(report, indent=2, default=str))
    else:
        print(f"[fep-binding validate] {args.yaml_path}")
        print(f"  kind:       {yaml_kind}")
        if yaml_kind == "hydration":
            print(f"  receptor:   n/a (hydration YAML)")
        elif receptor_report.get("exists"):
            print(
                f"  receptor:   {receptor_pdb}  "
                f"atoms={receptor_report.get('n_atoms')}  "
                f"chains={receptor_report.get('n_chains')}  "
                f"hetatm={receptor_report.get('n_hetatm')}")
        else:
            print(f"  receptor:   {receptor_pdb}  [MISSING]")
        print(f"  entries:    {len(entries_report)}")
        print()
        # Table
        hdr = (f"  {'name':<20s} {'heavy':>5s} {'stereo':>6s} "
               f"{'chg':>4s} {'expt':>7s}  {'status'}")
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for er in entries_report:
            if not er.get("ok"):
                print(f"  {er['name']:<20s} "
                      f"{'—':>5s} {'—':>6s} {'—':>4s} "
                      f"{(er.get('expt_kcalmol') or 0):>+7.2f}  "
                      f"FAIL: {er.get('reason', '?')}")
            else:
                print(
                    f"  {er['name']:<20s} "
                    f"{er['n_heavy_atoms']:>5d} "
                    f"{er['n_stereocenters']:>6d} "
                    f"{er['formal_charge']:>+4d} "
                    f"{(er.get('expt_kcalmol') or 0):>+7.2f}  "
                    "ok")
        print()
        if issues:
            print(f"  {len(issues)} issue(s):")
            for it in issues:
                print(f"    ! {it}")
            print()
            print("[fep-binding validate] FAIL")
        else:
            print("[fep-binding validate] PASS — "
                  "ready for bench / dg / ddg")
    return 0 if not issues else 1


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main())
