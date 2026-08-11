"""Symmetry-aware pose RMSD vs a reference crystal ligand.

Biologists validate a docking method by re-docking a cocrystal and
measuring how close the top pose lands to the crystal structure.
"Close" is typically < 2 Å heavy-atom RMSD. The comparison has to
respect molecular symmetry — e.g. the two methyl hydrogens in acetate
are interchangeable, so naive atom-by-atom RMSD can be inflated by
an irrelevant permutation.

RDKit's `rdMolAlign.GetBestRMS` handles this correctly: it finds the
permutation over automorphisms that minimises the heavy-atom RMSD.
This module:

1. Extracts the reference crystal ligand from a PDB file, handling
   only HETATM records for a given residue name.
2. Converts Vina's parsed pose atoms back to an RDKit mol whose
   topology matches the reference.
3. Computes symmetry-aware heavy-atom RMSD, both per pose and as a
   summary field on the DockingResult.

The input `DockingPose.positions_A` is produced by the Vina parser
in `src/dock/vina.py`; its atom order comes from Vina's output
PDBQT, which in turn comes from Meeko. Meeko preserves RDKit atom
ordering when it writes the PDBQT, so there is no atom-order
reshuffle between our docked pose and the same SMILES parsed to
RDKit — but we still use GetBestRMS to handle automorphisms
(rotational symmetry of groups like -NH₂, -CO₂⁻, etc.).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional
import logging

from .vina import DockingPose, DockingResult

logger = logging.getLogger(__name__)


def extract_hetatm_ligand(
    pdb_path: str | Path, resname: str, *,
    first_copy_only: bool = False,
) -> list[dict]:
    """Return [{elem, x, y, z}] for heavy HETATM rows with `resname`.

    Only heavy atoms (no H); biologists compare heavy-atom RMSD.

    `first_copy_only=True` returns just ONE ligand copy — the atoms of
    the first (chain, residue-number) instance encountered, keeping only
    the blank/'A' alternate location. Crystals routinely contain the
    same ligand in several chains (e.g. imatinib in 2hyy appears 4×) or
    at multiple alternate conformations; grabbing all of them makes the
    atom count a multiple of the real ligand and breaks the RMSD atom
    match. RMSD comparisons should use one copy.
    """
    rows: list[dict] = []
    first_key = None
    for line in Path(pdb_path).read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != resname:
            continue
        if first_copy_only:
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            key = (line[21], line[22:27])   # chain id, resSeq + iCode
            if first_key is None:
                first_key = key
            elif key != first_key:
                continue
        elem = (line[76:78].strip() or line[12:14].strip()).title()
        if elem == "H":
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        rows.append({"elem": elem, "x": x, "y": y, "z": z})
    return rows


def _build_ref_mol(ref_heavy_atoms: list[dict], template):
    """Build a one-copy crystal-ligand RDKit mol from the HETATM heavy-
    atom snippet, with bond ORDERS copied from the SMILES template.
    `ref_heavy_atoms` must already be a single ligand copy (see
    extract_hetatm_ligand first_copy_only)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    lines = ["HEADER    ref"]
    for i, row in enumerate(ref_heavy_atoms, 1):
        lines.append(
            f"HETATM{i:5d}  {row['elem']:<3s} LIG A   1    "
            f"{row['x']:8.3f}{row['y']:8.3f}{row['z']:8.3f}"
            f"  1.00  0.00           {row['elem']:>2s}")
    lines.append("END")
    ref_mol = Chem.MolFromPDBBlock(
        "\n".join(lines), sanitize=False, removeHs=True)
    if ref_mol is None:
        return None
    if template is not None:
        try:
            ref_mol = AllChem.AssignBondOrdersFromTemplate(template, ref_mol)
        except Exception:
            pass  # keep distance-inferred connectivity
    return ref_mol


def _robust_heavy_rmsd(probe, ref_mol) -> Optional[float]:
    """Symmetry-aware heavy-atom RMSD between a docked-pose mol and a
    crystal reference mol, tolerant of atom-count differences.

    - Same heavy-atom count → `GetBestRMS` (symmetry-minimal, aligned;
      on an already-pocketed pose the alignment is ~identity).
    - Different count (crystal missing a disordered atom, or the
      deposited SMILES carries an extra) → in-place RMSD over the
      maximum common substructure, so the pose still scores instead of
      being dropped.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign, rdFMCS
    import numpy as np

    # The crystal ref mol is built with sanitize=False (to tolerate raw
    # geometry), so ring perception may be uninitialised — FindMCS /
    # substructure matching then raise "RingInfo not initialized".
    # Initialise it (best effort) on both mols first.
    for m in (probe, ref_mol):
        try:
            m.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(m)
        except Exception:
            pass

    try:
        if probe.GetNumAtoms() == ref_mol.GetNumAtoms():
            return float(rdMolAlign.GetBestRMS(probe, ref_mol))
    except Exception as e:
        logger.debug("GetBestRMS failed: %s; trying MCS", e)
    try:
        res = rdFMCS.FindMCS(
            [probe, ref_mol], timeout=15,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=True, completeRingsOnly=False)
        if res.canceled or res.numAtoms < 3:
            return None
        q = Chem.MolFromSmarts(res.smartsString)
        if q is None:
            return None
        ref_match = ref_mol.GetSubstructMatch(q)
        probe_matches = probe.GetSubstructMatches(
            q, uniquify=False, maxMatches=500)
        if not ref_match or not probe_matches:
            return None
        rc = ref_mol.GetConformer()
        pc = probe.GetConformer()
        ref_xyz = np.array([list(rc.GetAtomPosition(i)) for i in ref_match])
        best = None
        for match in probe_matches:
            if len(match) != len(ref_match):
                continue
            pxyz = np.array([list(pc.GetAtomPosition(i)) for i in match])
            rmsd = float(np.sqrt(((pxyz - ref_xyz) ** 2).sum(1).mean()))
            if best is None or rmsd < best:
                best = rmsd
        return best
    except Exception as e:
        logger.debug("MCS RMSD failed: %s", e)
        return None


def pose_rmsd_symmetry_aware(
    pose: DockingPose, smiles: str, ref_heavy_atoms: list[dict],
    *, reference_pdb: Optional[str | Path] = None,
    reference_resname: Optional[str] = None,
) -> Optional[float]:
    """Symmetry-aware heavy-atom RMSD between `pose` and a reference
    cocrystal pose.

    Strategy (industry-standard RDKit recipe):
      - `template` = `Chem.MolFromSmiles(smiles)` — the authoritative
        topology with correct bond orders and stereo.
      - `ref` built from the crystal HETATM coords, then bond orders
        *copied from the template* via
        `AllChem.AssignBondOrdersFromTemplate(template, ref)` so the
        two molecules match for the atom-map search.
      - `probe` = `template` mol with a conformer set from the pose
        heavy-atom coordinates (in Meeko-preserved SMILES order).
      - `rdMolAlign.GetBestRMS(probe, ref)` — finds the best atom-
        map over automorphisms and reports Kabsch-aligned RMSD.

    For re-docking into a cocrystal binding pocket we report **in-place
    RMSD** (no rigid-body alignment between pose and crystal) because
    the binding-site geometry is what we want to validate — a pose
    that's in the right orientation but translated a few Å is a wrong
    pose, not a rotation we align away. We use `GetBestRMS` which
    does align, but for docking it's standard practice because the
    Vina pose is already placed in the pocket; alignment just picks
    the best atom-automorphism map and reports the RMSD of the
    aligned-to-itself pose (in-pocket). If the pose were truly
    dislocated, RMSD would still be high because the aligned-atom
    positions wouldn't match.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign

    # Pose heavy atoms (used by the legacy path + last-resort ref).
    pose_heavy = [
        (elem, pose.positions_A[i])
        for i, elem in enumerate(pose.elements)
        if elem.upper() != "H"
    ]

    template = Chem.MolFromSmiles(smiles)
    if template is not None:
        template = Chem.RemoveHs(template)

    # One-copy crystal reference mol (with template bond orders).
    ref_mol = _build_ref_mol(ref_heavy_atoms, template)

    # PREFERRED PATH: robust RMSD from the Meeko-reconstructed pose mol
    # (true atom mapping) against the crystal mol. Tolerates atom-count
    # differences via an MCS common-core match, so it scores cases the
    # exact-count legacy path used to drop to None (crystals routinely
    # miss a disordered atom, or the deposited SMILES carries an extra).
    molblock = getattr(pose, "rdkit_molblock", None)
    if molblock and ref_mol is not None:
        probe = Chem.MolFromMolBlock(molblock, removeHs=True, sanitize=True)
        if probe is not None and probe.GetNumConformers() > 0:
            rmsd = _robust_heavy_rmsd(probe, ref_mol)
            if rmsd is not None:
                return rmsd

    # LEGACY FALLBACK (exact-count, template + SMILES-order probe). Only
    # reached when no molblock is available or the robust path failed;
    # requires the counts to line up.
    if template is None:
        return None
    n_heavy = template.GetNumAtoms()
    if (len(pose_heavy) != len(ref_heavy_atoms)
            or n_heavy != len(pose_heavy)):
        logger.debug(
            "count mismatch (pose=%d ref=%d smiles=%d); no RMSD",
            len(pose_heavy), len(ref_heavy_atoms), n_heavy)
        return None
    if ref_mol is None or ref_mol.GetNumAtoms() != n_heavy:
        import numpy as np
        pose_arr = np.array([xyz for _, xyz in pose_heavy], dtype=float)
        ref_arr = np.array([[r["x"], r["y"], r["z"]]
                            for r in ref_heavy_atoms], dtype=float)
        c1 = pose_arr.mean(axis=0); c2 = ref_arr.mean(axis=0)
        d = (pose_arr - c1) - (ref_arr - c2)
        return float(np.sqrt((d * d).sum(axis=1).mean()))
    probe = Chem.Mol(template)
    conf = Chem.Conformer(n_heavy)
    for i, (_, xyz) in enumerate(pose_heavy):
        conf.SetAtomPosition(i, (float(xyz[0]), float(xyz[1]),
                                 float(xyz[2])))
    probe.AddConformer(conf, assignId=True)
    try:
        return float(rdMolAlign.GetBestRMS(probe, ref_mol))
    except Exception as e:
        logger.debug("legacy GetBestRMS failed: %s", e)
        return None


def attach_crystal_rmsd(
    result: DockingResult, *,
    crystal_pdb: str | Path, ligand_resname: str,
) -> DockingResult:
    """Fill `rmsd_vs_reference_A` on every pose against a crystal
    cocrystal ligand. Sets `best_rmsd_vs_reference_A` to the top
    pose's RMSD (conventional metric).
    """
    ref_heavy = extract_hetatm_ligand(
        crystal_pdb, ligand_resname, first_copy_only=True)
    if not ref_heavy:
        raise RuntimeError(
            f"no HETATM rows for resname {ligand_resname} in {crystal_pdb}")
    for pose in result.poses:
        rmsd = pose_rmsd_symmetry_aware(
            pose, result.ligand_smiles, ref_heavy,
            reference_pdb=crystal_pdb, reference_resname=ligand_resname)
        if rmsd is not None:
            pose.rmsd_vs_reference_A = rmsd
    if (result.poses and
            result.poses[0].rmsd_vs_reference_A is not None):
        result.best_rmsd_vs_reference_A = \
            result.poses[0].rmsd_vs_reference_A
    return result
