"""Export Vina docked poses to standard formats (SDF, PDB).

So biologists can open the poses in PyMOL / ChimeraX / Maestro
without leaving their existing workflow.

SDF is the canonical ligand container (one molecule block per
pose, with embedded properties like ΔG / K_d / crystal-RMSD).
PDB is lower-information but universally viewable.

Usage:
    from src.dock import DockingResult, export_poses_sdf
    r = dock_ligand(...)
    export_poses_sdf(r, "/tmp/biotin_poses.sdf")
    # → one multi-MOL SDF file. Drag-drop into PyMOL:
    #     File → Open → biotin_poses.sdf
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from .vina import DockingResult, DockingPose

logger = logging.getLogger(__name__)


def _pose_to_rdkit_mol(pose: DockingPose, smiles: str):
    """Reconstruct an RDKit Mol with pose coordinates.

    Returns an RDKit Mol (heavy atoms only — matches the ordering
    Meeko preserves through Vina), coordinates set from the pose.
    Returns None if RDKit can't parse the SMILES.
    """
    from rdkit import Chem

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return None
    template = Chem.RemoveHs(template)
    n_heavy = template.GetNumAtoms()

    pose_heavy_xyz = [pose.positions_A[i]
                      for i, e in enumerate(pose.elements)
                      if e.upper() != "H"]
    if len(pose_heavy_xyz) != n_heavy:
        logger.debug("pose heavy %d vs template heavy %d mismatch",
                     len(pose_heavy_xyz), n_heavy)
        return None

    mol = Chem.Mol(template)
    conf = Chem.Conformer(n_heavy)
    for i, (x, y, z) in enumerate(pose_heavy_xyz):
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    return mol


def export_poses_sdf(
    result: DockingResult,
    out_path: str | Path,
    *,
    include_properties: bool = True,
    top_k: Optional[int] = None,
    extra_props: Optional[dict] = None,
) -> int:
    """Write a multi-pose SDF.

    Returns the number of poses successfully written. Never raises
    on recoverable failures — bad poses are skipped with a debug log.

    Each pose block carries these SDF properties when
    `include_properties=True`:
        rank, affinity_kcalmol, affinity_kJmol, Kd_nM,
        crystal_rmsd_A (if known), pocket_ok, geometry_ok,
        pb_all_ok, seed, exhaustiveness, receptor_hash,
        ligand_inchi_key, cellsim_provenance.

    `extra_props` is an optional dict of additional SDF properties
    to attach to every pose (same value on each pose block). Used
    by src.dock.batch to carry the triage verdict + strain band
    into the SDF so biologists opening the file in PyMOL /
    ChimeraX see the pose-level wet-lab decision context.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import SDWriter
    except ImportError as e:
        logger.warning("RDKit missing; cannot export SDF: %s", e)
        return 0

    if not result.ok or not result.poses:
        logger.debug("nothing to export: ok=%s poses=%d",
                     result.ok, len(result.poses))
        return 0

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    poses = result.poses
    if top_k is not None:
        poses = poses[:top_k]

    w = SDWriter(str(out_path))
    written = 0
    for pose in poses:
        mol = _pose_to_rdkit_mol(pose, result.ligand_smiles)
        if mol is None:
            continue
        if include_properties:
            mol.SetProp("_Name",
                         f"{result.ligand_formula or 'lig'}_"
                         f"pose_{pose.mode}")
            mol.SetProp("rank", str(pose.mode))
            mol.SetProp("affinity_kcalmol",
                         f"{pose.affinity_kcalmol:.3f}")
            mol.SetProp("affinity_kJmol",
                         f"{pose.affinity_kJmol:.3f}")
            mol.SetProp("Kd_nM", f"{pose.kd_implied_nM:.3g}")
            if pose.rmsd_vs_reference_A is not None:
                mol.SetProp("crystal_rmsd_A",
                             f"{pose.rmsd_vs_reference_A:.3f}")
            if pose.posebusters_pocket_ok is not None:
                mol.SetProp("pocket_ok",
                             str(pose.posebusters_pocket_ok))
            if pose.posebusters_geometry_ok is not None:
                mol.SetProp("geometry_ok",
                             str(pose.posebusters_geometry_ok))
            if pose.posebusters_ok is not None:
                mol.SetProp("pb_all_ok", str(pose.posebusters_ok))
            if result.seed is not None:
                mol.SetProp("seed", str(result.seed))
            if result.exhaustiveness is not None:
                mol.SetProp("exhaustiveness",
                             str(result.exhaustiveness))
            if result.receptor_hash:
                mol.SetProp("receptor_hash", result.receptor_hash)
            if result.ligand_inchi_key:
                mol.SetProp("ligand_inchi_key", result.ligand_inchi_key)
            mol.SetProp("cellsim_provenance",
                         f"CellSim/{result.exhaustiveness}/"
                         f"seed={result.seed}")
            if extra_props:
                for k, v in extra_props.items():
                    if v is None or v == "":
                        continue
                    mol.SetProp(str(k), str(v))
        w.write(mol)
        written += 1
    w.close()
    logger.info("wrote %d poses to %s", written, out_path)
    return written


def export_poses_pdb(
    result: DockingResult,
    out_path: str | Path,
    *,
    top_k: Optional[int] = None,
) -> int:
    """Write all poses as a single MODEL-separated PDB file.

    Intended for users who don't have RDKit in their downstream
    viewer or just want a lightweight text format.
    """
    if not result.ok or not result.poses:
        return 0

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    poses = result.poses
    if top_k is not None:
        poses = poses[:top_k]

    written = 0
    with out_path.open("w") as f:
        f.write(f"HEADER    CellSim docking poses  "
                f"ligand={result.ligand_formula or '?'}  "
                f"receptor={Path(result.receptor_source).name}  "
                f"seed={result.seed}  "
                f"exh={result.exhaustiveness}\n")
        if result.ligand_inchi_key:
            f.write(f"REMARK  1 InChI-Key {result.ligand_inchi_key}\n")
        for pose in poses:
            f.write(f"MODEL    {pose.mode:>4d}\n")
            f.write(f"REMARK  2 AFFINITY {pose.affinity_kcalmol:+.3f} "
                    f"kcal/mol  K_d≈{pose.kd_implied_nM:.3g} nM\n")
            if pose.rmsd_vs_reference_A is not None:
                f.write(f"REMARK  3 CRYSTAL_RMSD "
                        f"{pose.rmsd_vs_reference_A:.3f} A\n")
            for i, (elem, xyz) in enumerate(
                    zip(pose.elements, pose.positions_A), 1):
                x, y, z = xyz
                f.write(
                    f"HETATM{i:>5d}  {elem:<3s} LIG A   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"  1.00  0.00           "
                    f"{elem:>2s}\n")
            f.write("ENDMDL\n")
            written += 1
        f.write("END\n")
    return written
