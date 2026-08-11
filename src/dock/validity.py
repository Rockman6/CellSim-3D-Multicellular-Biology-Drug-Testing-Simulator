"""PoseBusters physical-validity wrapper for docked poses.

Every biologist looking at a docking result wants to know: "is this
pose physically plausible, or is it structurally broken?"
PoseBusters answers that by running ~20 well-defined physical-
validity tests on a pose:

  Molecule-level (mol_pred alone):
    - sanitisation, stereochemistry, connectivity, molecular formula
    - reasonable bond lengths / angles, no atom clashes, planarity
      of aromatic rings, correct kekulisation, etc.

  Protein-conditional (mol_pred + mol_cond = receptor):
    - no ligand-protein steric clashes
    - no ligand atoms inside the protein

  Cocrystal-conditional (mol_pred + mol_true = crystal ligand):
    - heavy-atom RMSD ≤ 2 Å (our existing `pose_rmsd` metric)

Each test returns True/False; an "overall PB-ok" flag is True iff
every test passes. This is exactly what a biologist needs in a
one-glance summary.

We surface a per-pose `posebusters_ok` boolean on DockingPose and a
per-pose flag dict for a detailed breakdown. The viewer colour-codes
poses by this flag alongside crystal-RMSD.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional
import logging

from .vina import DockingPose, DockingResult

logger = logging.getLogger(__name__)


def _build_pose_mol(pose: DockingPose, smiles: str):
    """Build a chemically-correct 3D RDKit mol for a pose.

    Preferred path: the pose's `rdkit_molblock`, reconstructed via
    Meeko's reverse conversion, which carries the TRUE atom mapping +
    hydrogens. This is required for correctness — Meeko reorders atoms
    during prep, so the legacy template path below (which assumes the
    pose coords are in SMILES order) builds a topologically scrambled
    molecule and every geometry test fails. The template path is kept
    only as a degraded fallback for poses lacking a molblock.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molblock = getattr(pose, "rdkit_molblock", None)
    if molblock:
        m = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=True)
        if m is not None and m.GetNumConformers() > 0:
            # Meeko's reconstruction already includes hydrogens with
            # geometry; only add them if somehow absent.
            if not any(a.GetAtomicNum() == 1 for a in m.GetAtoms()):
                try:
                    m = Chem.AddHs(m, addCoords=True)
                except Exception:
                    pass
            return m
        logger.debug("MolFromMolBlock failed; falling back to template")

    tpl = Chem.MolFromSmiles(smiles)
    if tpl is None:
        return None
    tpl = Chem.RemoveHs(tpl)  # heavy-atom topology in SMILES order

    pose_heavy_xyz = [pose.positions_A[i]
                      for i, e in enumerate(pose.elements)
                      if e.upper() != "H"]
    n_heavy = tpl.GetNumAtoms()
    if n_heavy != len(pose_heavy_xyz):
        logger.debug("pose/template heavy atom mismatch: %d vs %d",
                     len(pose_heavy_xyz), n_heavy)
        return None

    conf = Chem.Conformer(n_heavy)
    for i, (x, y, z) in enumerate(pose_heavy_xyz):
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    tpl.AddConformer(conf, assignId=True)

    # Add hydrogens with geometry so PoseBusters' clash tests are
    # meaningful. `addCoords=True` places Hs at ideal geometry from
    # their parent heavy atoms.
    try:
        tpl_h = Chem.AddHs(tpl, addCoords=True)
    except Exception as e:
        logger.debug("AddHs(addCoords=True) failed: %s", e)
        return tpl  # fall back to heavy-atom-only (worse clash test)
    try:
        Chem.SanitizeMol(tpl_h)
    except Exception:
        pass
    return tpl_h


def _build_crystal_mol(pdb_path: Path, resname: str, smiles: str):
    """Build an RDKit mol for the crystal cocrystal ligand using the
    SMILES template + PDB HETATM coords. Returns None on failure."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Parse PDB, keep only the target HETATM residue.
    lines = [L for L in Path(pdb_path).read_text().splitlines()
             if L.startswith("HETATM") and L[17:20].strip() == resname]
    if not lines:
        return None
    lines = ["HEADER    crystal"] + lines + ["END"]
    ref = Chem.MolFromPDBBlock("\n".join(lines), sanitize=False,
                                removeHs=True)
    if ref is None:
        return None
    tpl = Chem.MolFromSmiles(smiles)
    if tpl is None:
        return None
    try:
        ref = AllChem.AssignBondOrdersFromTemplate(tpl, ref)
    except Exception as e:
        logger.debug("crystal AssignBondOrdersFromTemplate: %s", e)
    try:
        Chem.SanitizeMol(ref)
    except Exception:
        pass
    return ref


def _protein_only_pdb(pdb_path: Path, out_path: Path) -> Path:
    """Write a cleaned protein-only PDB (no heteroatoms, no waters).

    PoseBusters expects the receptor as a clean PDB/mol. We pre-clean
    it here so the "no clashes with protein" tests use only ATOM
    rows, not the cocrystal ligand or solvent.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with Path(pdb_path).open() as fi, out_path.open("w") as fo:
        for line in fi:
            if line.startswith(("HEADER", "TITLE", "REMARK")):
                fo.write(line)
            elif line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname in ("HOH", "WAT", "TIP3",
                               "NA", "CL", "NA+", "CL-",
                               "ZN", "MG", "CA", "MN", "FE"):
                    continue
                fo.write(line)
            elif line.startswith(("TER", "END")):
                fo.write(line)
    return out_path


def attach_posebusters(
    result: DockingResult,
    *,
    receptor_pdb: Optional[str | Path] = None,
    crystal_pdb: Optional[str | Path] = None,
    crystal_resname: Optional[str] = None,
) -> DockingResult:
    """Run PoseBusters on every pose and fill `posebusters_ok`.

    If a receptor is supplied, PoseBusters' protein-conditioned tests
    (no steric clashes, no ligand inside protein) are included. If a
    crystal PDB + resname is supplied, the cocrystal RMSD test is
    included.

    Fails silently (logs + keeps `posebusters_ok = None`) rather than
    raising, so a PoseBusters glitch on one pose never breaks a run.
    """
    try:
        from posebusters import PoseBusters
    except ImportError as e:
        logger.warning("posebusters not installed: %s", e)
        return result

    import tempfile

    # Optionally prepare a protein-only PDB (stripped of waters / ions
    # / heteroatoms).
    protein_mol: Optional[Path] = None
    crystal_mol = None
    with tempfile.TemporaryDirectory(prefix="cellsim-pb-") as tmp:
        tmp = Path(tmp)
        if receptor_pdb:
            protein_mol = _protein_only_pdb(
                Path(receptor_pdb), tmp / "receptor_only.pdb")
        if crystal_pdb and crystal_resname:
            crystal_mol = _build_crystal_mol(
                Path(crystal_pdb), crystal_resname, result.ligand_smiles)

        config = ("redock" if (protein_mol and crystal_mol)
                  else "dock" if protein_mol
                  else "mol")
        buster = PoseBusters(config=config, max_workers=0)

        for pose in result.poses:
            pose_mol = _build_pose_mol(pose, result.ligand_smiles)
            if pose_mol is None:
                pose.posebusters_ok = None
                continue
            try:
                kwargs: dict = {"mol_pred": pose_mol}
                if crystal_mol is not None:
                    kwargs["mol_true"] = crystal_mol
                if protein_mol is not None:
                    kwargs["mol_cond"] = str(protein_mol)
                df = buster.bust(**kwargs)
            except Exception as e:
                logger.debug("PoseBusters failed on pose %d: %s",
                             pose.mode, e)
                pose.posebusters_ok = None
                continue
            # `df` has one row with many True/False test columns.
            flags_series = df.iloc[0].drop(
                labels=["file", "molecule"], errors="ignore")
            flags_dict: dict = {}
            for k, v in flags_series.items():
                # Booleans / strings / numbers — keep anything simple.
                if isinstance(v, bool):
                    flags_dict[str(k)] = bool(v)
            pose.posebusters_flags = flags_dict

            # Pocket-fit: the set of tests that matter for a biologist
            # deciding "is this pose in the right place for wet-lab?"
            pocket_tests = [
                "mol_pred_loaded", "protein-ligand_maximum_distance",
                "minimum_distance_to_protein",
                "minimum_distance_to_organic_cofactors",
                "minimum_distance_to_inorganic_cofactors",
                "minimum_distance_to_waters",
                "volume_overlap_with_protein",
                "volume_overlap_with_organic_cofactors",
                "volume_overlap_with_inorganic_cofactors",
                "volume_overlap_with_waters",
            ]
            # Geometry: covalent + internal quality. Vina often fails
            # these because of approximate sterics; a short MD
            # post-refinement (future) fixes them.
            geometry_tests = [
                "sanitization", "all_atoms_connected", "no_radicals",
                "molecular_formula", "molecular_bonds",
                "double_bond_stereochemistry", "tetrahedral_chirality",
                "bond_lengths", "bond_angles", "internal_steric_clash",
                "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
                "double_bond_flatness", "internal_energy",
            ]

            def _all_pass(keys):
                present = [flags_dict[k] for k in keys if k in flags_dict]
                return bool(present) and all(present)

            pose.posebusters_pocket_ok = _all_pass(pocket_tests)
            pose.posebusters_geometry_ok = _all_pass(geometry_tests)
            pose.posebusters_ok = bool(flags_series.all())
    return result
