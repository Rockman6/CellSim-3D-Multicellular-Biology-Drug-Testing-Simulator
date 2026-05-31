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

Two interchangeable force-field stacks are available via
`compute_absolute_binding_dg(..., force_field_path=...)`:

- **"amber14"** (default): OpenMM classical ForceField
  `amber14-all.xml` + `amber14/tip3pfb.xml` for the protein and
  solvent, with `openmmforcefields.SMIRNOFFTemplateGenerator`
  providing the OpenFF Sage 2.1.0 + AM1-BCC ligand template.
  Fast and robust across all probed receptors (1ubq / 1stp /
  1m17 scaffold in under 15 s on CPU).
- **"smirnoff"**: pure `openff.interchange.Interchange.from_smirnoff`
  on `openff-2.1.0 + ff14sb_off_impropers + tip3p`. Kept for
  users who want provenance-identical Sage-bonded parameters
  between the hydration and binding legs; slower but viable
  after the terminal-residue bounding-box fix landed in
  `_prepare_protein_topology`.

Ligand parameters are identical between the two paths (Sage +
AM1-BCC), so the protein-leg stack is the only difference; ΔΔG
subtractions within the same receptor cancel either way.

Phase state returned on `BindingDGResult`:
- `'scaffolded_both_legs'` — systems built end-to-end, no MD.
- `'sampled'` — MD + MBAR complete, ΔG_bind populated. Pass
  `sample=True` to reach this state; expect hours per compound
  on GPU.

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
    # GHMC acceptance per window per leg — prof's checklist
    # requirement. < 0.70 signals timestep too large and ΔG
    # untrustworthy. Passed through for the fep-report gate.
    ghmc_acceptance_complex: list = field(default_factory=list)
    ghmc_acceptance_solvent: list = field(default_factory=list)

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


# Platform throughput anchors for the wall-time estimator.
# Units: steps × atoms per second, measured empirically and used
# as throughput[platform] / n_atoms → effective steps/sec for a
# system of size n_atoms. The anchors assume a 2 fs timestep with
# GHMC + HMR constraints (standard CellSim profile); large atomic
# systems drop linearly in steps/sec with atom count.
#
# Anchors calibrated against:
#  - CPU:    24-atom methane in vacuum, single-core measured at
#            ~80 k steps/sec → 2 M steps·atom/sec
#  - Metal:  ~2 k-atom FreeSolv M5 Max run, ~13 k steps/sec
#            effective on the full Phase-2 sampler → ~26 M
#  - CUDA:   ~150 M is the paper-quality anchor for an H100
#            running dhfr-class 23k-atom systems at ~6 k steps/sec
#
# These are anchors with ±2-3× error bars. Reported as such in
# the CLI output so biologists don't treat them as precise.
_PLATFORM_STEP_ATOMS_PER_SECOND = {
    "cpu": 2_000_000,
    "metal_m5max": 26_000_000,
    "cuda_h100": 150_000_000,
}

# Wall-time amortisation for per-λ overheads (minimisation,
# context switches, MBAR eval). Empirically ~15-30% of pure MD
# wall, so add a 25% surcharge.
_OVERHEAD_SURCHARGE = 1.25


def estimate_sampling_wall_hours(
    n_atoms: int,
    n_windows: int,
    n_production_steps: int,
    n_equilibration_steps: int,
    n_legs: int = 2,
) -> dict[str, float]:
    """Rough estimate of sample_alchemical_windows wall time across
    reference platforms. Returns {'cpu': h, 'metal_m5max': h,
    'cuda_h100': h}. Accuracy: ±2-3× — this is a sanity-check
    preview, not a guarantee.

    Model: total_steps = n_legs × n_windows × (n_equil + n_prod).
    Effective throughput on a given platform scales inversely with
    n_atoms; multiply by a 25% overhead surcharge for per-λ
    minimise + MBAR eval.
    """
    total_steps = n_legs * n_windows * (
        n_equilibration_steps + n_production_steps)
    out: dict[str, float] = {}
    for plat, const in _PLATFORM_STEP_ATOMS_PER_SECOND.items():
        steps_per_sec = const / max(n_atoms, 1)
        seconds = total_steps / max(steps_per_sec, 1e-6)
        seconds *= _OVERHEAD_SURCHARGE
        out[plat] = seconds / 3600.0
    return out


def format_wall_estimate_block(
    est: dict[str, float], *, gate_hours: float = 48.0,
) -> str:
    """Pretty-print the wall estimate as a biologist-readable block.
    Flags anything over `gate_hours` on CPU as "CPU-only is not
    viable" — a 48-hour run on CPU usually means the user forgot
    to flag GPU, not that they want a CPU run of that length."""
    def _fmt(h: float) -> str:
        if h < 1.0:
            return f"{h * 60:.0f} min"
        if h < 48.0:
            return f"{h:.1f} h"
        return f"{h / 24:.1f} d"

    lines = [
        "  ~ Estimated sampling wall time "
        "(±2-3× accuracy; for planning, not guarantees):",
        f"      CPU            : {_fmt(est['cpu'])}",
        f"      Metal (M5 Max) : {_fmt(est['metal_m5max'])}",
        f"      CUDA H100      : {_fmt(est['cuda_h100'])}",
    ]
    if est["cpu"] > gate_hours:
        lines.append(
            "  ! CPU-only is not viable for this config "
            "(> 48 h). Use Metal or CUDA.")
    return "\n".join(lines)


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


# Charged-ligand finite-size PBC correction (Rocklin 2013 J Chem
# Phys 139:184103) is NOT yet implemented. The full correction is
# 5 terms (EMP + NET + RIP + DSC + ANA) and requires knowing the
# protein's net charge + ion pairing distribution — too much
# auxiliary state to derive from the YAML alone in a way I trust
# to ship publication-grade numbers.
#
# Practical guidance for biologists hitting the validator's
# 'formal_charge ≠ 0' message:
#   - For a CONGENERIC SERIES of equally-charged ligands (e.g.
#     all-anionic carboxylates), the PBC self-interaction error
#     cancels in ΔΔG. Rank-correlation gates (Kendall τ) are still
#     valid; absolute MAE is offset by a constant.
#   - For absolute ΔG_bind on a charged ligand vs neutral protein,
#     extend the Rocklin correction post-hoc per the cited paper.
#     A typical -1 ligand in a 5 nm cubic TIP3P box adds ~5-15
#     kcal/mol depending on protein charge.
#   - When the Phase-3 implementation lands, this comment block
#     gets replaced with an analytical _rocklin_pbc_correction
#     function + autowiring into compute_absolute_binding_dg.


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
    # allow_undefined_stereo=True so SMILES like 2-iminobiotin
    # (C=N exocyclic geometry ambiguous) don't block a bench run.
    # RDKit will pick a deterministic but arbitrary stereo; the
    # validator warns so the biologist knows to pin explicit E/Z
    # if ΔG_bind depends on it.
    lig = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
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
    # allow_undefined_stereo=True so SMILES like 2-iminobiotin
    # (C=N exocyclic geometry ambiguous) don't block a bench run.
    # RDKit will pick a deterministic but arbitrary stereo; the
    # validator warns so the biologist knows to pin explicit E/Z
    # if ΔG_bind depends on it.
    lig = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
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
    use_replica_exchange: bool = False,
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
            sample_stride=sample_stride, seed=seed,
            use_replica_exchange=use_replica_exchange)
        if not cx_r.ok:
            result.reason = f"complex leg sample failed: {cx_r.reason}"
            result.wall_seconds = time.time() - t0
            return result

        sv_r = sample_alchemical_windows(
            solv_alch, solv_top, solv_pos,
            n_windows=n_windows,
            n_equilibration_steps=n_equilibration_steps,
            n_production_steps=n_production_steps,
            sample_stride=sample_stride, seed=seed,
            use_replica_exchange=use_replica_exchange)
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
    result.ghmc_acceptance_complex = list(cx_r.ghmc_acceptance)
    result.ghmc_acceptance_solvent = list(sv_r.ghmc_acceptance)
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
    force_field_path: str = "amber14",
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
        sample_stride=sample_stride, seed=seed,
        force_field_path=force_field_path)
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
        sample_stride=sample_stride, seed=seed,
        force_field_path=force_field_path)
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
    """CLI: `cellsim fep-binding` — four subcommands.

      cellsim fep-binding validate <yaml>
          Sub-second YAML + SMILES hygiene check (RDKit + Lipinski +
          stereo + pKa + unit check + wall-time estimate with
          local-GPU detection). Run this BEFORE any multi-hour
          sampled run.

      cellsim fep-binding dg <SMILES> <receptor.pdb>
          Absolute ΔG_bind for a single compound. --sample drives
          MD; without --sample, scaffold-only phase returns.

      cellsim fep-binding ddg <SMI_A> <SMI_B> <receptor.pdb>
          Relative ΔΔG(A→B) via two independent absolute-binding
          runs + subtraction.

      cellsim fep-binding bench <yaml>
          YAML-driven batch with crash-proof incremental CSV,
          --resume for restarting from a crash, per-compound ETA,
          Ctrl-C handler with partial-CSV preservation.

    Scaffold mode (no --sample) runs in seconds per compound on
    CPU; sampled mode requires GPU hours for real-benchmark YAMLs.
    """
    import argparse
    import json as _json
    import sys

    ap = argparse.ArgumentParser(
        description="Milestone B FEP — absolute and relative "
                    "binding ΔG. Typical workflow: (1) validate — "
                    "sub-second YAML + SMILES hygiene check; "
                    "(2) bench — build & optionally sample on GPU; "
                    "(3) dg / ddg — single-compound variants.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    dgp = sub.add_parser("dg", help="absolute ΔG_bind")
    dgp.add_argument("smiles")
    dgp.add_argument("receptor_pdb")
    dgp.add_argument(
        "--n-windows", type=int, default=11,
        help="alchemical λ-windows per leg (default 11)")
    dgp.add_argument(
        "--softcore-alpha", type=float, default=0.5,
        help="softcore LJ smoothing (default 0.5)")
    dgp.add_argument(
        "--padding", type=float, default=1.2,
        help="water padding around complex in nm (default 1.2)")
    dgp.add_argument("--restraint-k", type=float, default=4184.0,
                     help="CoM harmonic spring in kJ/mol/nm² "
                          "(default 4184 ≈ 10 kcal/mol/Å²)")
    dgp.add_argument("--sample", action="store_true",
                     help="run MD (slow on CPU; intended for GPU)")
    dgp.add_argument("--force-field-path", default="amber14",
                     choices=["amber14", "smirnoff"],
                     help="amber14 (default; fast, AMBER14 protein "
                          "+ tip3pfb + SMIRNOFF ligand) or "
                          "smirnoff (pure Interchange, slower but "
                          "provenance-identical with hydration)")
    dgp.add_argument("--json", action="store_true")

    ddgp = sub.add_parser("ddg", help="relative ΔΔG(A→B)")
    ddgp.add_argument("smiles_a")
    ddgp.add_argument("smiles_b")
    ddgp.add_argument("receptor_pdb")
    ddgp.add_argument(
        "--n-windows", type=int, default=11,
        help="alchemical λ-windows per leg (default 11)")
    ddgp.add_argument(
        "--softcore-alpha", type=float, default=0.5,
        help="softcore LJ smoothing (default 0.5)")
    ddgp.add_argument(
        "--padding", type=float, default=1.2,
        help="water padding around complex in nm (default 1.2)")
    ddgp.add_argument(
        "--restraint-k", type=float, default=4184.0,
        help="CoM harmonic spring in kJ/mol/nm² (default 4184)")
    ddgp.add_argument(
        "--sample", action="store_true",
        help="run MD (slow on CPU; intended for GPU)")
    ddgp.add_argument("--force-field-path", default="amber14",
                      choices=["amber14", "smirnoff"],
                      help="protein/ligand parametrisation path "
                           "(default amber14)")
    ddgp.add_argument("--json", action="store_true")

    bp = sub.add_parser(
        "bench",
        help="YAML-driven batch (e.g. binding_streptavidin.yaml)")
    bp.add_argument(
        "yaml_path",
        help="path to a binding YAML under benchmarks/fep/ — "
             "run `cellsim fep-binding validate` on it first")
    bp.add_argument(
        "--n-windows", type=int, default=11,
        help="number of alchemical λ-windows per leg "
             "(default 11 — matches Milestone A production)")
    bp.add_argument(
        "--softcore-alpha", type=float, default=0.5,
        help="softcore LJ smoothing parameter (default 0.5; "
             "raise to 1.0 for less-kinked PMF on hard systems)")
    bp.add_argument(
        "--padding", type=float, default=1.2,
        help="water padding around complex in nm (default 1.2; "
             "a 1.0 nm nonbonded cutoff wants padding ≥ 1.0)")
    bp.add_argument(
        "--restraint-k", type=float, default=4184.0,
        help="CoM harmonic restraint spring in kJ/mol/nm² "
             "(default 4184 ≈ 10 kcal/mol/Å², standard Boresch k)")
    bp.add_argument(
        "--sample", action="store_true",
        help="run MD on every entry (HOURS per compound on CPU; "
             "a full streptavidin set is a GPU-only operation)")
    bp.add_argument(
        "--production-steps", type=int, default=2000,
        help="GHMC production steps per window (default 2000 ≈ "
             "4 ps at 2 fs; Milestone A uses 25 000 = 50 ps)")
    bp.add_argument(
        "--equilibration-steps", type=int, default=500,
        help="GHMC equilibration steps per window (default 500 = "
             "1 ps; Milestone A uses 2 500 = 5 ps)")
    bp.add_argument(
        "--sample-stride", type=int, default=100,
        help="sample interval in production steps (default 100 = "
             "200 samples per 2 000-step window)")
    bp.add_argument(
        "--force-field-path", default="amber14",
        choices=["amber14", "smirnoff"],
        help="protein/ligand parametrisation path (default "
             "amber14; 'smirnoff' for the pure-Interchange "
             "fallback that's slower but provenance-matches "
             "hydration's solvent leg)")
    bp.add_argument(
        "--out-csv", default=None,
        help="write per-entry results (same schema as FreeSolv so "
             "`cellsim fep-report` can consume it)")
    bp.add_argument(
        "--resume", action="store_true",
        help="if --out-csv already exists, skip compounds already "
             "present with a non-empty dG_pred_kcalmol. Intended "
             "for restarting a crashed multi-hour run without "
             "losing completed compounds.")
    bp.add_argument(
        "--replica-exchange", action="store_true",
        help="use openmmtools.multistate.ReplicaExchangeSampler "
             "(Hamiltonian replica exchange between adjacent "
             "λ-states) instead of the hand-rolled independent-"
             "replica path. The canonical literature fix for MBAR "
             "overlap failures on tight intramolecular charge "
             "networks (acetic_acid, acetamide, biotin); see "
             "BENCHMARKS.md § 'Replica exchange unblocks tight "
             "binders'. Costs ~same wall time as the hand-rolled "
             "path at equal sampling budget but converges where "
             "hand-rolled returns MBAR overlap failure.")

    vp = sub.add_parser(
        "validate",
        help="sub-second YAML dry-run: parse SMILES via RDKit, "
             "confirm receptor PDB exists, report issues — run "
             "this BEFORE launching a multi-hour sampled run")
    vp.add_argument(
        "yaml_path",
        help="path to a benchmark YAML (binding_*.yaml or "
             "freesolv_*.yaml under benchmarks/fep/)")
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
            sample=args.sample,
            force_field_path=args.force_field_path)
    elif args.cmd == "ddg":
        r = compute_relative_binding_ddg(
            args.smiles_a, args.smiles_b, args.receptor_pdb,
            n_windows=args.n_windows,
            softcore_alpha=args.softcore_alpha,
            padding_nm=args.padding,
            restraint_k_kJ_per_nm2=args.restraint_k,
            sample=args.sample,
            force_field_path=args.force_field_path)
    elif args.cmd == "bench":
        return _run_bench(args)
    else:  # validate
        return _run_validate(args)

    if args.json:
        from dataclasses import asdict
        print(_json.dumps(asdict(r), indent=2, default=str))
    else:
        print(r.summary())
        # Biologist preview: when a scaffold-only run completes
        # (no --sample), print how long a sampled run of this config
        # would take on the reference platforms. Prevents "I'll just
        # try --sample overnight" surprises that turn into 5-day CPU
        # burns.
        if (not getattr(args, "sample", False)
                and args.cmd == "dg"
                and isinstance(r, BindingDGResult)
                and r.ok
                and r.n_total_atoms_complex):
            est = estimate_sampling_wall_hours(
                n_atoms=r.n_total_atoms_complex,
                n_windows=args.n_windows,
                n_production_steps=25000,
                n_equilibration_steps=2500)
            print()
            print("  (reference: 11 × 25 000 prod + 2 500 equil "
                  "× 2 legs — Milestone-A-tier sampling)")
            print(format_wall_estimate_block(est))
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

    # --resume: if the output CSV already exists and contains
    # compounds with a non-empty dG_pred_kcalmol, seed `rows` with
    # those prior entries and skip them when iterating. The biologist
    # restarting a crashed 6-hour EGFR run recovers without losing
    # the completed compounds.
    completed_names: set[str] = set()
    if args.resume and args.out_csv and Path(args.out_csv).exists():
        try:
            with Path(args.out_csv).open(
                    "r", encoding="utf-8-sig", newline="") as fi:
                reader = _csv.DictReader(fi)
                for prior in reader:
                    pred = (prior.get("dG_pred_kcalmol") or "").strip()
                    if pred and pred != "None":
                        # Normalise numeric columns (CSV dict has all
                        # strings).
                        def _f(key):
                            v = (prior.get(key) or "").strip()
                            if not v or v == "None":
                                return None
                            try:
                                return float(v)
                            except ValueError:
                                return None
                        rows.append({
                            "name": prior.get("name", "").strip(),
                            "smiles": prior.get("smiles", "").strip(),
                            "dG_expt_kcalmol": _f("dG_expt_kcalmol"),
                            "dG_pred_kcalmol": _f("dG_pred_kcalmol"),
                            "uncertainty_kcalmol":
                                _f("uncertainty_kcalmol"),
                            "residual_kcalmol": _f("residual_kcalmol"),
                            "ghmc_accept_mean":
                                _f("ghmc_accept_mean"),
                            "ghmc_accept_min": _f("ghmc_accept_min"),
                            "wall_seconds": _f("wall_seconds"),
                            "ok": (prior.get("ok", "").strip().lower()
                                   in ("true", "1")),
                            "reason": prior.get("reason", "").strip(),
                        })
                        completed_names.add(
                            prior.get("name", "").strip())
        except OSError as e:
            print(f"bench --resume: cannot read {args.out_csv}: {e}",
                  file=_sys.stderr)
            return 2

    t_all = _time.time()
    print(f"[fep-binding bench] {args.yaml_path}")
    print(f"  receptor:   {receptor_pdb}")
    print(f"  entries:    {len(entries)}")
    print(f"  sample:     {bool(args.sample)}")
    if completed_names:
        print(f"  resumed:    {len(completed_names)} compounds "
              f"already in {args.out_csv} — skipping "
              f"{sorted(completed_names)}")
    print()

    interrupted = False
    for entry in entries:
        name = entry.get("name", "?")
        smi = entry.get("smiles", "")
        expt = entry.get("dG_bind_kcalmol")
        if name in completed_names:
            print(f"[bench] {name}  (skipped, resumed from CSV)",
                  flush=True)
            continue
        print(f"[bench] {name}  {smi}", flush=True)
        ts = _time.time()
        try:
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
                sample_stride=args.sample_stride,
                force_field_path=args.force_field_path,
                use_replica_exchange=getattr(
                    args, "replica_exchange", False))
        except KeyboardInterrupt:
            # Graceful exit: the crash-proof CSV (previous commit)
            # already has every completed compound. Just announce
            # the halt + partial-CSV location so biologist knows
            # there's something to resume from.
            interrupted = True
            break
        wall = _time.time() - ts
        pred = r.dG_bind_kcalmol
        unc = r.uncertainty_kcalmol
        resid = (pred - expt
                 if pred is not None and expt is not None else None)
        # Progress ETA: based on the compounds completed so far in
        # THIS bench (skipped/resumed ones don't count — no new
        # wall observed). Format 'n/m done; eta ~Xm|h remaining'.
        # Biologist can step away and estimate when to come back.
        active_rows = [
            row for row in rows
            if row.get("name") not in completed_names]
        n_done_now = len(active_rows) + 1  # includes current
        n_remaining = max(
            0,
            len(entries) - len(completed_names) - n_done_now)
        elapsed_now = _time.time() - t_all
        if n_done_now >= 1 and n_remaining > 0:
            mean_per = elapsed_now / n_done_now
            eta_s = mean_per * n_remaining
            if eta_s < 90:
                eta_tag = f"eta ~{eta_s:.0f}s"
            elif eta_s < 3600:
                eta_tag = f"eta ~{eta_s/60:.1f}m"
            else:
                eta_tag = f"eta ~{eta_s/3600:.1f}h"
            progress = (
                f"  [{n_done_now}/"
                f"{len(entries) - len(completed_names)} active; "
                f"{eta_tag}]")
        else:
            progress = ""

        if r.ok and pred is None:
            print(f"  scaffolded {r.phase}  "
                  f"atoms={r.n_total_atoms_complex}  "
                  f"wall={wall:.1f}s{progress}", flush=True)
        elif r.ok:
            print(f"  pred = {pred:+.2f}  expt = {expt:+.2f}  "
                  f"resid = {(resid or 0):+.2f}  "
                  f"wall={wall:.1f}s{progress}", flush=True)
        else:
            print(f"  FAIL: {r.reason}{progress}", flush=True)

        # Minimum across both legs is the conservative per-compound
        # acceptance; if either leg was under-sampled the ΔG bind
        # doesn't stand. Matches hydration's same-schema convention.
        all_accepts = (
            list(r.ghmc_acceptance_complex)
            + list(r.ghmc_acceptance_solvent))
        ghmc_mean = (
            sum(all_accepts) / len(all_accepts)
            if all_accepts else None)
        ghmc_min = min(all_accepts) if all_accepts else None

        rows.append({
            "name": name,
            "smiles": smi,
            "dG_expt_kcalmol": expt,
            "dG_pred_kcalmol": pred,
            "uncertainty_kcalmol": unc,
            "residual_kcalmol": resid,
            "ghmc_accept_mean": ghmc_mean,
            "ghmc_accept_min": ghmc_min,
            "wall_seconds": round(wall, 2),
            "ok": r.ok,
            "reason": r.reason,
        })
        # Write CSV after EVERY compound finishes, not just at the
        # end. On a multi-hour bench a biologist might Ctrl-C, lose
        # power, or hit a segfault; they still want the CSV with
        # completed compounds. Overwrite-whole-file is O(n) in rows
        # per compound but bench is typically < 20 compounds so this
        # is < 1 ms overhead per compound — unmeasurable vs the 10+
        # min sample time.
        if args.out_csv:
            _write_bench_csv(args.out_csv, rows)

    t_total = _time.time() - t_all
    n_ok = sum(1 for row in rows if row["ok"])
    print()
    if interrupted:
        n_remaining = (
            len(entries) - len(completed_names) - len(
                [row for row in rows
                 if row.get("name") not in completed_names]))
        print(f"[fep-binding bench] INTERRUPTED after "
              f"{n_ok}/{len(rows)} completed "
              f"({n_remaining} of YAML's entries skipped).  "
              f"total wall {t_total:.1f} s  "
              f"({args.yaml_path})")
    else:
        print(f"[fep-binding bench] {n_ok}/{len(rows)} ok  "
              f"total wall {t_total:.1f} s  "
              f"({args.yaml_path})")

    if args.out_csv:
        # Final re-write ensures even a 0-entry YAML produces the
        # header (no compounds → incremental writer never fired).
        _write_bench_csv(args.out_csv, rows)
        print(f"  wrote {args.out_csv}  "
              f"(schema compatible with `cellsim fep-report`)")
        if interrupted:
            print(f"  to continue: cellsim fep-binding bench "
                  f"{args.yaml_path} --out-csv {args.out_csv} "
                  "--resume [other args]")
    elif interrupted and rows:
        # No --out-csv was set. The in-memory results are lost
        # on process exit; at least tell the biologist what
        # happened so they re-run with the flag.
        print(f"  partial results for {len(rows)} compounds "
              "NOT saved (no --out-csv flag). Re-run with e.g. "
              "--out-csv bench_run.csv so the next Ctrl-C "
              "preserves completed compounds.")

    # Interrupted → exit 130 (conventional Ctrl-C code, >1 for
    # CI to detect).
    if interrupted:
        return 130
    return 0 if n_ok == len(rows) else 1


def _write_bench_csv(out_csv, rows) -> None:
    """Atomically write the bench CSV (overwrite whole file).
    Called after every compound finishes so a Ctrl-C / power-off
    mid-bench still leaves the partial CSV on disk."""
    import csv as _csv
    out = Path(out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["name", "smiles", "dG_expt_kcalmol",
            "dG_pred_kcalmol", "uncertainty_kcalmol",
            "residual_kcalmol", "ghmc_accept_mean",
            "ghmc_accept_min", "wall_seconds", "ok", "reason"]
    # Write to a temp sibling file first then rename — atomic on
    # POSIX, so a crash during the write doesn't corrupt the
    # previous incremental CSV.
    tmp = out.with_suffix(out.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8-sig") as fo:
        w = _csv.DictWriter(
            fo, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    import os
    os.replace(tmp, out)


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

    # Severity split: errors fail the validate exit; warnings
    # print but don't block. A biologist running validate on a
    # drug-like ligand (e.g. lapatinib with 11 rotatable bonds)
    # should get a heads-up, not a refusal to proceed.
    issues: list[str] = []   # errors — make this a hard FAIL
    warnings: list[str] = []  # advisories — print but exit 0
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

    # Unit-sanity check on binding_site_center_nm. A common
    # biologist error is copy-pasting a coordinate from a docking
    # YAML (which uses Ångströms) into a binding FEP YAML (which
    # expects nm). Protein bounding boxes are 3–10 nm per side; a
    # value > 10 nm is almost certainly an Å number that should
    # have been divided by 10. Warn loudly with the suggested
    # nm conversion. Second check: is the coordinate actually
    # inside the protein's ATOM-record bounding box? If not, the
    # restraint would anchor at a phantom location.
    bsc = recv.get("binding_site_center_nm")
    if bsc is not None and yaml_kind != "hydration":
        try:
            bsc_xyz = [float(c) for c in bsc]
            if any(abs(c) > 10.0 for c in bsc_xyz):
                nm_conv = [round(c / 10.0, 4) for c in bsc_xyz]
                warnings.append(
                    f"binding_site_center_nm = {bsc_xyz} has a "
                    f"component > 10 nm — likely Å, not nm. Did "
                    f"you mean {nm_conv}? (Protein bounding boxes "
                    "are rarely > 10 nm. Set to null to let fpocket "
                    "auto-detect the druggable pocket.)")
            receptor_report["binding_site_center_nm"] = bsc_xyz

            # In-bounds check — only if we successfully parsed the
            # receptor PDB atoms above.
            if receptor_pdb and Path(receptor_pdb).is_file():
                try:
                    xs, ys, zs = [], [], []
                    for line in Path(receptor_pdb).read_text(
                            encoding="utf-8",
                            errors="replace").splitlines():
                        if not line.startswith("ATOM"):
                            continue
                        # PDB columns 31-54 are X/Y/Z (8-char floats)
                        try:
                            xs.append(float(line[30:38]))
                            ys.append(float(line[38:46]))
                            zs.append(float(line[46:54]))
                        except ValueError:
                            continue
                    if xs:
                        # Å → nm, pad by 1.5 nm (accommodates
                        # solvent-exposed surface sites).
                        lo = [min(xs)/10 - 1.5, min(ys)/10 - 1.5,
                              min(zs)/10 - 1.5]
                        hi = [max(xs)/10 + 1.5, max(ys)/10 + 1.5,
                              max(zs)/10 + 1.5]
                        off_axis = [
                            i for i, c in enumerate(bsc_xyz)
                            if c < lo[i] or c > hi[i]]
                        if off_axis:
                            ax = ["x", "y", "z"]
                            names = ",".join(ax[i] for i in off_axis)
                            warnings.append(
                                f"binding_site_center_nm {bsc_xyz} "
                                f"is outside the protein's atom-"
                                f"coord bounding box (lo={lo}, "
                                f"hi={hi} nm) on axis {names}. "
                                "Possible wrong-receptor paste or "
                                "wrong frame; verify by opening the "
                                "PDB in PyMOL and checking the site.")
                except OSError:
                    pass
        except (TypeError, ValueError):
            issues.append(
                f"binding_site_center_nm must be a 3-element list "
                f"of floats in nm; got {bsc!r}")

    entries = data.get("entries") or []
    if not entries:
        issues.append("no entries in YAML")

    # Allowed per-entry key names — anything else is a probable
    # typo (dG_bond_kcalmol vs dG_bind_kcalmol, Smiles vs smiles)
    # which would silently skip the field and give the biologist
    # a ΔG of None or a missing SMILES.
    _ALLOWED_ENTRY_KEYS = {
        "name", "smiles",
        "dG_bind_kcalmol", "dG_hydration_kcalmol",
        "IC50_M",
        "expt_unc_kcalmol",
        "source_notes",
    }

    seen_names: set[str] = set()
    for i, entry in enumerate(entries):
        name = entry.get("name", f"<entry_{i}>")
        smi = (entry.get("smiles") or "").strip()
        expt = entry.get("dG_bind_kcalmol",
                         entry.get("dG_hydration_kcalmol"))
        er: dict = {"name": name, "smiles": smi, "expt_kcalmol": expt}

        # Typo detection: any key not in the allowed set is
        # suspicious. Don't FAIL — authors may legitimately extend
        # the schema with custom fields (e.g. series tag, assay
        # ID) — just WARN with the specific key.
        if isinstance(entry, dict):
            extra = set(entry.keys()) - _ALLOWED_ENTRY_KEYS
            for bad_key in sorted(extra):
                # Suggest the closest known key via a tiny
                # Levenshtein-ish heuristic (no stdlib dep).
                def _similar(a, b):
                    a, b = a.lower(), b.lower()
                    if a == b:
                        return 100
                    common = sum(1 for c in a if c in b)
                    return common * 100 // max(len(a), len(b))
                suggestion = max(
                    _ALLOWED_ENTRY_KEYS, key=lambda k: _similar(
                        k, bad_key))
                if _similar(suggestion, bad_key) >= 70:
                    warnings.append(
                        f"{name}: unknown entry key '{bad_key}' — "
                        f"did you mean '{suggestion}'? (Will be "
                        "silently ignored by the bench driver.)")
                else:
                    warnings.append(
                        f"{name}: unknown entry key '{bad_key}' "
                        "(silently ignored; custom fields are OK "
                        "but may hide typos)")

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
                    # Detect ambiguous E/Z bonds using RDKit's
                    # own CIP-aware stereo finder. This is the
                    # same primitive openff-toolkit uses, so if
                    # it returns a bond here the Molecule.from_smiles
                    # path WILL reject it without
                    # allow_undefined_stereo=True. Catches imino
                    # (iminobiotin N=C) without false-positiving on
                    # amide C=O or urea C=N where one end has no
                    # heavy neighbour.
                    n_undef_bonds = 0
                    try:
                        potentials = Chem.FindPotentialStereo(mol)
                        n_undef_bonds = sum(
                            1 for p in potentials
                            if p.type == Chem.StereoType.Bond_Double
                            and p.specified ==
                                Chem.StereoSpecified.Unspecified)
                    except Exception:
                        # Old RDKit without FindPotentialStereo —
                        # skip the check (better than false alarms).
                        pass
                    # Rotatable bond count + HBD/HBA — biologists
                    # already read these in Lipinski terms, and
                    # they predict sampling convergence cost
                    # (flexible ligands need more MD).
                    from rdkit.Chem import rdMolDescriptors as _rmd
                    try:
                        n_rot = _rmd.CalcNumRotatableBonds(mol)
                        n_hbd = _rmd.CalcNumHBD(mol)
                        n_hba = _rmd.CalcNumHBA(mol)
                        tpsa = _rmd.CalcTPSA(mol)
                    except Exception:
                        n_rot = n_hbd = n_hba = 0
                        tpsa = 0.0

                    # Per-compound scaffold + sampled wall estimates.
                    # Constants calibrated from 2026-04-22 1ubq
                    # smoke: 82.5 s for ~1200 MCMC steps × 15k-atom
                    # complex → ~5e-7 s per atom-step on CPU. GPU
                    # conservative 20× speedup (Metal/OpenCL on
                    # M-series).
                    receptor_atoms = (
                        receptor_report.get("n_atoms") or 0)
                    # Scaffold estimate (AM1-BCC + OpenMM parametrise
                    # + packmol): empirically ~3 s/heavy-atom.
                    scaffold_s = max(5.0, 3.0 * n_heavy)
                    # Complex atom count guess: receptor × 10 for
                    # padding 0.8 nm (fits 1ubq 15k, 1stp 25k, 1m17
                    # 110k from today's benches).
                    complex_atoms = max(15_000, receptor_atoms * 10)
                    # Default Milestone-A-ish sampled params:
                    # 11 windows × 25 k prod steps × 2 legs = 550 k
                    # integration steps total per compound. Plus
                    # ~10% equilibration overhead.
                    n_steps = 11 * (25_000 + 2_500) * 2
                    cpu_s_per_atom_step = 5e-7  # from 1ubq smoke
                    gpu_s_per_atom_step = 2.5e-8
                    sampled_cpu_s = (
                        complex_atoms * n_steps * cpu_s_per_atom_step)
                    sampled_gpu_s = (
                        complex_atoms * n_steps * gpu_s_per_atom_step)
                    er.update({
                        "ok": True,
                        "canon_smiles": canon,
                        "n_heavy_atoms": n_heavy,
                        "n_stereocenters": n_stereo,
                        "n_undef_double_bonds": n_undef_bonds,
                        "formal_charge": formal_charge,
                        "n_rotatable_bonds": n_rot,
                        "n_hbd": n_hbd,
                        "n_hba": n_hba,
                        "tpsa_A2": tpsa,
                        "est_scaffold_wall_s": scaffold_s,
                        "est_sampled_wall_cpu_s": sampled_cpu_s,
                        "est_sampled_wall_gpu_s": sampled_gpu_s,
                    })
                    # Charged ligand: PBC self-interaction error in
                    # absolute ΔG. Cancels in ΔΔG within a same-charge
                    # series (so still valid for rank-correlation
                    # gates) but makes absolute MAE off by 5-15
                    # kcal/mol depending on protein charge. Warning,
                    # not error — biologists studying ionic ligand
                    # series shouldn't be blocked, just informed.
                    if formal_charge != 0:
                        warnings.append(
                            f"{name}: formal charge "
                            f"{formal_charge:+d} — absolute ΔG_bind "
                            "carries a 5-15 kcal/mol PBC self-"
                            "interaction error (Rocklin 2013 J Chem "
                            "Phys 139:184103). For ΔΔG within a "
                            "same-charge congeneric series the error "
                            "cancels; for absolute ΔG apply Rocklin "
                            "post-hoc.")
                    if n_heavy > 60:
                        issues.append(
                            f"{name}: {n_heavy} heavy atoms > 60 — "
                            "AM1-BCC will be slow (minutes/compound)")
                    if n_undef_bonds:
                        issues.append(
                            f"{name}: {n_undef_bonds} acyclic "
                            "double bond(s) with undefined E/Z — "
                            "builder will pick arbitrarily; pin "
                            "explicit geometry in SMILES if ΔG "
                            "depends on it (e.g. iminobiotin's "
                            "exocyclic C=N)")
                    # Advisory warning (not a FAIL): rotatable-bond
                    # count suggests sampling difficulty. Lipinski
                    # convention: > 10 is 'very flexible'. For FEP
                    # this means standard 50 ps/window may not
                    # converge; recommend longer production steps.
                    # Goes into `warnings` so exit stays 0 — the
                    # biologist gets a heads-up, not a block.
                    if er.get("n_rotatable_bonds", 0) > 10:
                        warnings.append(
                            f"{name}: {er['n_rotatable_bonds']} "
                            "rotatable bonds > 10 — highly flexible "
                            "ligand. FEP may not converge at default "
                            "50 ps/window. Consider "
                            "--production-steps 100000 (200 ps)")
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
        "warnings": warnings,
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
        # Table (wider: Lipinski-ish descriptors biologists read).
        hdr = (
            f"  {'name':<20s} {'heavy':>5s} {'rot':>4s} "
            f"{'HBD':>3s} {'HBA':>3s} {'TPSA':>5s} "
            f"{'stereo':>6s} {'chg':>4s} "
            f"{'expt':>7s}  {'status'}")
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for er in entries_report:
            if not er.get("ok"):
                print(f"  {er['name']:<20s} "
                      f"{'—':>5s} {'—':>4s} {'—':>3s} {'—':>3s} "
                      f"{'—':>5s} {'—':>6s} {'—':>4s} "
                      f"{(er.get('expt_kcalmol') or 0):>+7.2f}  "
                      f"FAIL: {er.get('reason', '?')}")
            else:
                print(
                    f"  {er['name']:<20s} "
                    f"{er['n_heavy_atoms']:>5d} "
                    f"{er.get('n_rotatable_bonds', 0):>4d} "
                    f"{er.get('n_hbd', 0):>3d} "
                    f"{er.get('n_hba', 0):>3d} "
                    f"{er.get('tpsa_A2', 0):>5.0f} "
                    f"{er['n_stereocenters']:>6d} "
                    f"{er['formal_charge']:>+4d} "
                    f"{(er.get('expt_kcalmol') or 0):>+7.2f}  "
                    "ok")
        print()
        # Hydration wall-time estimate. Constants are empirical
        # from friend's in-flight M5 Max run on binding=freesolv_12
        # (observed 4-6 h GPU for the full 12-compound set → ~25
        # min/compound GPU at production params). Small systems
        # (1500-atom solvated ligand) have large per-step OpenMM
        # fixed overhead that dwarfs the per-atom cost calibrated
        # from a larger-box binding smoke. Base-per-step model:
        #   per_step_s = 0.003 s (GPU) or 0.08 s (CPU)
        # matches the observed small-hydration throughput.
        if yaml_kind in ("hydration", "mixed"):
            ok_entries = [e for e in entries_report if e.get("ok")]
            if ok_entries:
                def _fmt_wall(s):
                    if s < 90:
                        return f"{s:.0f} s"
                    if s < 3600:
                        return f"{s/60:.1f} min"
                    return f"{s/3600:.1f} h"
                total_scaffold = sum(
                    max(5.0, 1.5 * e.get("n_heavy_atoms", 0))
                    for e in ok_entries)
                # 11 windows × 25 k prod + 2.5 k equil × 2 legs.
                n_steps_per_cpd = 11 * (25_000 + 2_500) * 2
                total_cpu = (
                    len(ok_entries) * n_steps_per_cpd * 0.08)
                total_gpu = (
                    len(ok_entries) * n_steps_per_cpd * 0.003)
                print("  estimated wall time (full YAML, 11 "
                      "windows × 25 000 prod steps × 2 legs):")
                print(f"    scaffold only : {_fmt_wall(total_scaffold)}")
                print(f"    sampled (CPU) : {_fmt_wall(total_cpu)}")
                print(f"    sampled (GPU) : {_fmt_wall(total_gpu)}  "
                      "(M5 Max observed 4-6 h on 12 compounds)")
                print()

        # Wall-time budget summary (only meaningful for binding
        # YAMLs — hydration YAMLs skip protein build).
        if yaml_kind in ("binding", "mixed"):
            ok_entries = [e for e in entries_report if e.get("ok")]
            if ok_entries:
                total_scaffold = sum(
                    e.get("est_scaffold_wall_s", 0)
                    for e in ok_entries)
                total_cpu = sum(
                    e.get("est_sampled_wall_cpu_s", 0)
                    for e in ok_entries)
                total_gpu = sum(
                    e.get("est_sampled_wall_gpu_s", 0)
                    for e in ok_entries)

                def _fmt_wall(s):
                    if s < 90:
                        return f"{s:.0f} s"
                    if s < 3600:
                        return f"{s/60:.1f} min"
                    return f"{s/3600:.1f} h"

                # Detect which OpenMM platform is actually usable
                # on THIS machine so the biologist reads the right
                # column. OpenMM ships Reference + CPU everywhere;
                # Metal / OpenCL / CUDA appear only where hardware
                # + drivers are installed.
                try:
                    from openmm import Platform as _Plat
                    _plats = {
                        _Plat.getPlatform(i).getName()
                        for i in range(_Plat.getNumPlatforms())}
                except Exception:
                    _plats = {"CPU"}
                _gpu_name = next(
                    (n for n in ("Metal", "CUDA", "OpenCL")
                     if n in _plats), None)
                _rec_wall = total_gpu if _gpu_name else total_cpu
                _rec_note = (
                    f"THIS machine: {_gpu_name} available → "
                    f"~{_fmt_wall(_rec_wall)} expected"
                    if _gpu_name else
                    f"THIS machine: no GPU detected → use CPU "
                    f"wall (~{_fmt_wall(_rec_wall)}); rent a GPU "
                    "for production runs")

                print("  estimated wall time (full YAML, 11 "
                      "windows × 25 000 prod steps × 2 legs):")
                print(f"    scaffold only : {_fmt_wall(total_scaffold)}")
                print(f"    sampled (CPU) : {_fmt_wall(total_cpu)}  "
                      "(single-threaded, not recommended for real runs)")
                print(f"    sampled (GPU) : {_fmt_wall(total_gpu)}  "
                      "(M-series Metal / CUDA estimate)")
                print(f"    → {_rec_note}")
                print()

        if warnings:
            print(f"  {len(warnings)} warning(s):")
            for w in warnings:
                print(f"    ? {w}")
            print()
        if issues:
            print(f"  {len(issues)} issue(s):")
            for it in issues:
                print(f"    ! {it}")
            print()
            print("[fep-binding validate] FAIL")
        else:
            tag = (" (with warnings)" if warnings else "")
            print(f"[fep-binding validate] PASS{tag} — "
                  "ready for bench / dg / ddg")
    return 0 if not issues else 1


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main())
