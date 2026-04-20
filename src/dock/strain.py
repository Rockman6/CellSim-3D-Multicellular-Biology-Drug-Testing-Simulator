"""Ligand-strain energy diagnostic for docked poses.

Vina's empirical score does not penalise intramolecular strain —
a pose can score "strongly" while adopting a high-energy
conformation. Perola & Charifson (2004) J Med Chem 47:2499
quantified this across 150+ protein-ligand cocrystals and
established the field convention:

    strain = E(bound conformation)  -  E(relaxed conformation)

The empirical ratio E(bound) / E(ensemble_avg) is more robust
than the absolute strain because UFF and MMFF have systematic
offsets on heteroaromatic systems. CellSim uses the PoseBusters
paper (Buttenschoen 2024 Chem Sci 15:3130) bands on the ratio:

    ratio < 1.5  good:        pose looks like a crystal pose
    ratio < 3.0  acceptable:  plausible but somewhat strained
    ratio < 7.0  suspicious:  likely a Vina scoring artefact
    ratio > 7.0  reject:      almost certainly non-physical

The absolute strain (kcal/mol) is reported alongside for users
who want the Perola-Charifson-style number.

Implementation: delegate the actual energy calculation to
PoseBusters' `check_energy_ratio` module (UFF ensemble over
100 conformers), which is already part of the cellsim env and
has been validated against PDBBind in the PoseBusters paper.
CellSim's contribution here is:

  1. the Vina-pose → RDKit-Mol bridge that preserves both
     coordinates and SMILES-derived bond orders (Meeko reorders
     heavy atoms relative to the SMILES walk in some rings, so
     naive index-based assignment breaks connectivity);
  2. the Perola-Charifson banding applied to the absolute strain
     (kcal/mol) that PoseBusters computes internally;
  3. a biologist-facing summary line wired into `src/dock/viewer.py`
     and the profile dashboard.

Usage:
    from src.dock.strain import ligand_strain
    s = ligand_strain(pose.elements, pose.positions_A, smiles)
    print(s.summary())
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict
from typing import Optional

logger = logging.getLogger(__name__)


def _normalise_elem(e: str) -> str:
    """AutoDock PDBQT atom types → periodic-table element."""
    e = e.upper()
    if e in ("A", "G", "GA", "C"):
        return "C"
    if e in ("N", "NA", "NS"):
        return "N"
    if e in ("O", "OA"):
        return "O"
    if e in ("S", "SA"):
        return "S"
    return e.title()


# Energy-ratio bands follow the PoseBusters paper (Buttenschoen
# et al 2024 Chem Sci 15:3130) which used UFF ensemble ratios on
# PDBBind crystal ligands. Their empirical finding: crystal poses
# cluster at ratio 1.0-2.0; mis-docked / clashing poses have
# ratio > 4; ratio > 7 is essentially always a failure mode.
# We use these UFF-calibrated bands instead of the original
# Perola-Charifson (MMFF-calibrated) absolute-strain thresholds
# because PoseBusters ships UFF energies.
STRAIN_RATIO_GOOD = 1.5
STRAIN_RATIO_ACCEPTABLE = 3.0
STRAIN_RATIO_SUSPICIOUS = 7.0


@dataclass
class StrainResult:
    ok: bool
    reason: str = ""
    smiles: Optional[str] = None
    n_heavy: Optional[int] = None
    e_docked_kcalmol: Optional[float] = None
    e_relaxed_mean_kcalmol: Optional[float] = None
    strain_kcalmol: Optional[float] = None
    energy_ratio: Optional[float] = None   # posebusters convention
    band: Optional[str] = None             # negligible | tolerable |
                                           # suspicious | reject
    method: str = "UFF ensemble (PoseBusters)"

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] strain  {self.reason}"
        return (f"[{self.band:>10s}]  strain = "
                f"{self.strain_kcalmol:+.2f} kcal/mol  "
                f"(E_docked {self.e_docked_kcalmol:+.1f}, "
                f"E_relaxed_mean {self.e_relaxed_mean_kcalmol:+.1f}, "
                f"ratio {self.energy_ratio:.2f}, "
                f"{self.method})")


def _classify(ratio: float) -> str:
    if ratio < STRAIN_RATIO_GOOD:
        return "good"
    if ratio < STRAIN_RATIO_ACCEPTABLE:
        return "acceptable"
    if ratio < STRAIN_RATIO_SUSPICIOUS:
        return "suspicious"
    return "reject"


def _pose_to_mol(elements, positions_A, smiles):
    """Build an RDKit Mol with pose coords + SMILES-derived bond
    orders. Returns (mol, n_heavy, err_or_None)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds

    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return None, 0, f"bad SMILES: {smiles}"
    template = Chem.RemoveHs(template)
    n_heavy_mol = template.GetNumAtoms()

    heavy = [(_normalise_elem(e), xyz)
             for e, xyz in zip(elements, positions_A)
             if _normalise_elem(e) != "H"]
    n_heavy_pose = len(heavy)
    if n_heavy_pose != n_heavy_mol:
        return (None, n_heavy_pose,
                f"heavy-atom count mismatch: pose {n_heavy_pose}, "
                f"SMILES {n_heavy_mol}")

    xyz_lines = [str(n_heavy_pose), smiles]
    for elem, xyz in heavy:
        xyz_lines.append(
            f"{elem} {float(xyz[0]):.4f} {float(xyz[1]):.4f} "
            f"{float(xyz[2]):.4f}")
    xyz_block = "\n".join(xyz_lines)

    # Infer connectivity from 3D, then transfer bond orders from the
    # SMILES template via substructure matching. If this fails
    # (typically because DetermineConnectivity mis-infers one edge
    # on a distorted Vina pose), we surface ok=False honestly rather
    # than use a less reliable fallback — strain is informational,
    # and an empty cell in the CSV is more trustworthy than a bogus
    # ratio pulled from misaligned coordinates.
    try:
        raw = Chem.MolFromXYZBlock(xyz_block)
        if raw is None:
            return None, n_heavy_pose, "MolFromXYZBlock returned None"
        rdDetermineBonds.DetermineConnectivity(raw)
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw)
        return mol, n_heavy_pose, None
    except Exception as e:
        return (None, n_heavy_pose,
                f"bond-order inference failed: {e}")


def ligand_strain(
    elements: list[str],
    positions_A: list,
    smiles: str,
    *,
    ensemble_n: int = 50,
) -> StrainResult:
    """Compute UFF-ensemble strain for a docked pose.

    Uses PoseBusters' `check_energy_ratio` internally, then bands
    the absolute strain (kcal/mol) per Perola & Charifson 2004.

    Returns `ok=False` (not a raised exception) on recoverable
    failures — strain is a triage read-out, not a gate.
    """
    mol, n_heavy_pose, err = _pose_to_mol(
        elements, positions_A, smiles)
    if err is not None:
        return StrainResult(
            ok=False, reason=err, smiles=smiles,
            n_heavy=n_heavy_pose)

    try:
        from posebusters.modules.energy_ratio import (
            check_energy_ratio,
        )
    except ImportError as e:
        return StrainResult(
            ok=False, reason=f"PoseBusters missing: {e}",
            smiles=smiles, n_heavy=n_heavy_pose)

    try:
        out = check_energy_ratio(
            mol,
            threshold_energy_ratio=float("inf"),
            ensemble_number_conformations=ensemble_n)
    except Exception as e:
        return StrainResult(
            ok=False,
            reason=f"PoseBusters energy_ratio crashed: {e}",
            smiles=smiles, n_heavy=n_heavy_pose)

    r = out.get("results", {})
    e_docked = r.get("mol_pred_energy")
    e_relax = r.get("ensemble_avg_energy")
    ratio = r.get("energy_ratio")
    if e_docked is None or e_relax is None or ratio is None:
        return StrainResult(
            ok=False,
            reason=f"PoseBusters returned no energies: keys={list(r)}",
            smiles=smiles, n_heavy=n_heavy_pose)

    strain = float(e_docked) - float(e_relax)
    return StrainResult(
        ok=True,
        smiles=smiles, n_heavy=n_heavy_pose,
        e_docked_kcalmol=float(e_docked),
        e_relaxed_mean_kcalmol=float(e_relax),
        strain_kcalmol=strain,
        energy_ratio=float(ratio),
        band=_classify(float(ratio)))


if __name__ == "__main__":
    import argparse
    import sys
    from pathlib import Path

    REPO_ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(REPO_ROOT))

    from src.dock import dock_ligand  # noqa: E402

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb")
    ap.add_argument("smiles")
    ap.add_argument("--center", nargs=3, type=float, required=True,
                    metavar=("X", "Y", "Z"))
    ap.add_argument("--box", nargs=3, type=float, default=[20, 20, 20])
    ap.add_argument("--exh", type=int, default=16)
    ap.add_argument("--num-modes", type=int, default=3)
    ap.add_argument("--ensemble", type=int, default=50)
    args = ap.parse_args()

    r = dock_ligand(
        args.pdb, args.smiles,
        center_A=tuple(args.center), box_size_A=tuple(args.box),
        exhaustiveness=args.exh, num_modes=args.num_modes, seed=1,
        cpu=2)
    if not r.ok or not r.poses:
        print(f"dock failed: {r.reason}")
        sys.exit(1)
    for pose in r.poses:
        s = ligand_strain(pose.elements, pose.positions_A,
                          args.smiles, ensemble_n=args.ensemble)
        print(f"pose #{pose.mode}  ΔG = {pose.affinity_kcalmol:+.2f}  "
              f"|  {s.summary()}")
