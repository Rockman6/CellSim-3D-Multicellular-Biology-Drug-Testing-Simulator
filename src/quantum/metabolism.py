"""CYP3A4 site-of-metabolism (SoM) predictor via homolytic C-H BDE.

Non-AI physics-based approach. For each C-H bond in the molecule:

    BDE(C-H) = E(R·) + E(H·) - E(R-H)

where E(R-H) is the xTB GFN2 energy of the full molecule, E(R·) is
the xTB energy of the radical left after removing that H (doublet,
multiplicity = 2), and E(H·) is the isolated hydrogen radical energy
at the same xTB level (a well-defined constant ≈ −10.7 eV at GFN2).

The C-H bond with the lowest BDE is typically the primary site of
metabolism for radical-abstraction CYP enzymes (of which CYP3A4 is
the canonical example — responsible for ~50 % of phase-I drug
metabolism).

References:
  - Rydberg & Olsen 2012 (SMARTCyp, J Chem Inf Model 52:2919) —
    original rule-based CYP3A4 SoM method using DFT-precomputed
    C-H activation energies, ~80 % top-3 accuracy.
  - Kromann 2017 (RegioSQM, Chem Sci 8:660) — semi-empirical
    (PM6) protonation-energy method; shown to agree with BDE-
    based rankings.
  - de Groot 2006 (Drug Metab Rev 38:251) — CYP3A4 mechanism
    discussion.

Limitations (honest documentation):

  - BDE alone ignores steric accessibility. A buried C-H with a
    low BDE may still be protected by the substrate's binding
    orientation. A full model would combine BDE with accessibility
    (distance to the CYP heme iron after docking into the CYP
    active site). This MVP reports BDE-ranked SoM candidates; a
    future PR adds "CYP3A4 active-site docking" to refine.
  - N-H bonds (amines, amides) are handled but often not the
    primary CYP substrate.
  - S-H, P-H, O-H bonds can also be CYP substrates in rare cases;
    we rank them too for completeness but flag them as secondary.

The output is a ranked list of candidate abstraction sites with
BDE values in kcal/mol, per-atom coordinates for the viewer,
and a "top-3" selector that biologists can check against
literature.
"""

from __future__ import annotations

import logging
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

# Hydrogen-radical E in Hartree at GFN2 — a well-defined constant
# (independent of molecule; computed once from a standalone H atom
# xtb single-point with mult=2). Value taken from a calibration run:
# xtb -> -0.393405 Hartree = -10.7047 eV. If GFN2 changes versions
# we recompute on demand (see below).
_H_RADICAL_HARTREE_DEFAULT = -0.393405

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.quantum.xtb import (  # noqa: E402
    XtbResult, _parse_xtb_out, _tool_versions, _HARTREE_TO_EV,
)

logger = logging.getLogger(__name__)


@dataclass
class SoMCandidate:
    """One candidate abstraction site."""

    # Heavy-atom (carbon / N / O / S / P) bearing the H that gets
    # abstracted.
    parent_atom_idx: int           # 0-indexed in the full molecule
    parent_element: str            # "C", "N", ...
    hydrogen_atom_idx: int         # 0-indexed in the full molecule
    position_A: list = field(default_factory=list)  # parent heavy xyz
    # Energetics
    bde_Hartree: Optional[float] = None
    bde_kcalmol: Optional[float] = None
    rank: Optional[int] = None     # 1 = primary SoM (lowest BDE)

    def summary(self) -> str:
        return (f"rank={self.rank:>2d}  {self.parent_element}"
                f"(idx={self.parent_atom_idx})  "
                f"BDE={self.bde_kcalmol:>6.1f} kcal/mol")


@dataclass
class SoMResult:
    """Envelope for a CYP3A4 SoM prediction."""

    smiles: str
    ok: bool
    reason: str = ""
    method: str = "GFN2-xTB BDE"
    enzyme: str = "CYP3A4"

    candidates: list = field(default_factory=list)  # list[SoMCandidate]
    elements: list = field(default_factory=list)
    positions_A: list = field(default_factory=list)
    e_mol_Hartree: Optional[float] = None
    e_H_radical_Hartree: Optional[float] = None

    tool_versions: dict = field(default_factory=dict)
    wall_seconds: Optional[float] = None

    def top_k(self, k: int = 3) -> list:
        """Return top-k lowest-BDE candidates."""
        return [c for c in self.candidates if c.rank is not None
                and c.rank <= k]

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] SoM {self.smiles}  {self.reason}"
        lines = [f"[OK]   SoM {self.enzyme}  {self.smiles}  "
                 f"({self.wall_seconds:.1f} s, "
                 f"{len(self.candidates)} C-H candidates)"]
        for c in self.top_k(5):
            lines.append("  " + c.summary())
        return "\n".join(lines)

    def as_dict(self) -> dict:
        d = asdict(self)
        return d


def _xtb_single_point_from_xyz_text(
    xyz_text: str, *, charge: int, multiplicity: int,
    method: str = "gfn2", timeout_s: int = 120,
    workdir: Optional[Path] = None,
) -> Optional[float]:
    """Low-level xtb call given pre-formatted XYZ text. Returns the
    total energy in Hartree, or None on failure."""
    import shutil
    import subprocess
    import tempfile

    xtb_bin = shutil.which("xtb")
    if xtb_bin is None:
        return None
    ctx = (tempfile.TemporaryDirectory(prefix="cellsim-bde-")
           if workdir is None else None)
    try:
        wd = Path(workdir or ctx.name)  # type: ignore[union-attr]
        xyz_path = wd / "mol.xyz"
        xyz_path.write_text(xyz_text)
        cmd = [xtb_bin, str(xyz_path),
               "--" + method,
               "--chrg", str(charge),
               "--uhf", str(max(0, multiplicity - 1))]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=timeout_s, cwd=str(wd))
        except subprocess.TimeoutExpired:
            return None
        if r.returncode != 0:
            logger.debug("xtb failed rc=%d: %s",
                         r.returncode, (r.stderr or r.stdout)[-200:])
            return None
        parsed = _parse_xtb_out(r.stdout)
        e = parsed.get("total_energy_Hartree")
        return e
    finally:
        if ctx is not None:
            ctx.cleanup()


def _compute_h_radical_hartree() -> float:
    """Run an xtb single-point on a standalone H atom (doublet) to
    get the current GFN2 reference for E(H·)."""
    xyz = "1\nH radical\nH     0.000000     0.000000     0.000000\n"
    e = _xtb_single_point_from_xyz_text(
        xyz, charge=0, multiplicity=2)
    return e if e is not None else _H_RADICAL_HARTREE_DEFAULT


def _build_xyz_text(title: str, elements: list, positions: list) -> str:
    lines = [f"{len(elements)}", title]
    for sym, (x, y, z) in zip(elements, positions):
        lines.append(f"{sym:<4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


def predict_cyp_som_bde(
    smiles: str,
    *,
    include_elements: tuple = ("C", "N", "O", "S", "P"),
    top_k: int = 5,
    seed: int = 1,
    timeout_s: int = 120,
    recompute_h_ref: bool = False,
) -> SoMResult:
    """Predict CYP3A4 SoM candidates via homolytic C-H BDE.

    Parameters
    ----------
    smiles : str
        Input compound.
    include_elements : tuple
        Heavy-atom elements whose X-H bonds are considered. Default
        covers CYP-relevant C/N/O/S/P-H. For CYP3A4 specifically,
        carbon C-H dominates.
    top_k : int
        Return ranks 1..top_k candidates.
    seed : int
        RDKit ETKDGv3 embed seed.
    timeout_s : int
        Per-xtb-call timeout.
    recompute_h_ref : bool
        If True, run an xtb single-point on a standalone H atom once
        at the start to get the exact GFN2 reference E(H·). If False
        (default), use the cached value from our calibration run.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    result = SoMResult(
        smiles=smiles, ok=False,
        tool_versions=_tool_versions())
    t0 = time.time()

    # Build 3D molecule with explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result.reason = f"RDKit could not parse SMILES: {smiles}"
        return result
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        result.reason = "RDKit ETKDGv3 embed failed"
        return result
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    elements: list[str] = []
    positions: list[list[float]] = []
    for i, atom in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        elements.append(atom.GetSymbol())
        positions.append([p.x, p.y, p.z])
    result.elements = elements
    result.positions_A = positions

    # E(R-H): full molecule.
    full_xyz = _build_xyz_text(f"full {smiles}", elements, positions)
    e_mol = _xtb_single_point_from_xyz_text(
        full_xyz, charge=0, multiplicity=1, timeout_s=timeout_s)
    if e_mol is None:
        result.reason = "xtb failed on full molecule"
        return result
    result.e_mol_Hartree = e_mol

    # E(H·): run once per call if requested, else use cached.
    if recompute_h_ref:
        e_H = _compute_h_radical_hartree()
    else:
        e_H = _H_RADICAL_HARTREE_DEFAULT
    result.e_H_radical_Hartree = e_H

    # Enumerate X-H bonds where X is in include_elements.
    candidates: list[SoMCandidate] = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in include_elements:
            continue
        for neighbour in atom.GetNeighbors():
            if neighbour.GetAtomicNum() != 1:
                continue
            h_idx = neighbour.GetIdx()
            parent_idx = atom.GetIdx()
            # Build radical XYZ: remove the hydrogen at h_idx.
            rad_elements = [e for i, e in enumerate(elements)
                            if i != h_idx]
            rad_positions = [p for i, p in enumerate(positions)
                             if i != h_idx]
            rad_xyz = _build_xyz_text(
                f"rad at idx={h_idx}", rad_elements, rad_positions)
            e_rad = _xtb_single_point_from_xyz_text(
                rad_xyz, charge=0, multiplicity=2, timeout_s=timeout_s)
            if e_rad is None:
                logger.debug(
                    "xtb failed on radical (H idx %d); skipping", h_idx)
                continue
            bde_h = e_rad + e_H - e_mol
            cand = SoMCandidate(
                parent_atom_idx=parent_idx,
                parent_element=atom.GetSymbol(),
                hydrogen_atom_idx=h_idx,
                position_A=list(positions[parent_idx]),
                bde_Hartree=bde_h,
                bde_kcalmol=bde_h * 627.509474,
            )
            candidates.append(cand)

    if not candidates:
        result.reason = "no X-H bonds of requested element types"
        result.wall_seconds = time.time() - t0
        return result

    # Rank ascending by BDE (lowest = primary SoM).
    candidates.sort(key=lambda c: c.bde_kcalmol)
    for rank, c in enumerate(candidates, 1):
        c.rank = rank

    result.candidates = candidates
    result.wall_seconds = time.time() - t0
    result.ok = True
    return result


def verify_som_dft(
    result: SoMResult,
    *,
    top_n: int = 3,
    functional: str = "B3LYP",
    basis: str = "def2-svp",
) -> SoMResult:
    """Re-rank the top-N xTB-BDE candidates at DFT level.

    Takes an existing SoMResult from predict_cyp_som_bde, recomputes
    BDE for the top-N candidates using PySCF DFT (paper-grade),
    and overwrites those candidates' `bde_Hartree` / `bde_kcalmol`
    fields with the DFT values. Re-ranks the candidate list by the
    new BDE. Remaining candidates keep their xTB values.

    Biologist interpretation: xTB BDE has a ~20 kcal/mol systematic
    offset vs experiment; DFT gets within ~2-5 kcal/mol. If the
    xTB top-1 matches the DFT top-1, the xTB ranking is trustworthy
    for this compound. If they disagree, DFT wins.

    Non-AI; PySCF Kohn-Sham physics end-to-end.
    """
    if not result.ok:
        return result
    import time

    # Lazy import; DFT is heavy.
    try:
        from src.quantum.dft import dft_single_point
    except ImportError as e:
        logger.warning("DFT module unavailable: %s", e)
        return result

    # Need full atom-level coords. We have them on result.positions_A.
    # But we also need them per radical (remove one H). Rebuild
    # radicals from scratch using the same pathway as the xTB code.
    top = result.candidates[:top_n]
    if not top:
        return result

    # Full-molecule DFT energy.
    t0 = time.time()
    full = dft_single_point(
        result.smiles, functional=functional, basis=basis,
        charge=0, multiplicity=1)
    if not full.ok:
        logger.warning("DFT on full molecule failed: %s", full.reason)
        return result
    e_mol_Ha = full.total_energy_Hartree

    # H-atom reference at DFT.
    h_ref = _dft_h_atom_energy(functional=functional, basis=basis)
    if h_ref is None:
        logger.warning("DFT H-atom reference failed; skipping re-rank")
        return result

    # For each candidate, compute radical energy.
    # We use the xTB-cached radical geometry (result.positions_A
    # minus the hydrogen at c.hydrogen_atom_idx) and do a DFT
    # single-point (UKS mult=2) on it.
    from pyscf import gto, dft as pyscf_dft

    updated = list(result.candidates)
    for i, cand in enumerate(top):
        h_idx = cand.hydrogen_atom_idx
        rad_elements = [e for j, e in enumerate(result.elements)
                        if j != h_idx]
        rad_positions = [p for j, p in enumerate(result.positions_A)
                         if j != h_idx]
        try:
            mol = gto.M(
                atom=[(sym, (x, y, z)) for sym, (x, y, z) in
                      zip(rad_elements, rad_positions)],
                basis=basis, charge=0, spin=1,
                max_memory=4000, verbose=0)
            mf = pyscf_dft.UKS(mol)
            mf.xc = functional.lower()
            mf.max_cycle = 60
            mf.conv_tol = 1e-8
            e_rad = mf.kernel()
        except Exception as e:
            logger.warning("DFT radical SCF failed for H%d: %s",
                           h_idx, str(e)[:80])
            continue
        if e_rad is None or not getattr(mf, "converged", False):
            continue
        bde_Ha = e_rad + h_ref - e_mol_Ha
        bde_kcal = bde_Ha * 627.509474
        # Update the candidate in place.
        rank_pos = updated.index(cand)
        updated[rank_pos] = SoMCandidate(
            parent_atom_idx=cand.parent_atom_idx,
            parent_element=cand.parent_element,
            hydrogen_atom_idx=cand.hydrogen_atom_idx,
            position_A=cand.position_A,
            bde_Hartree=float(bde_Ha),
            bde_kcalmol=float(bde_kcal),
            rank=None)  # re-ranked below

    # Re-rank by BDE ascending.
    updated.sort(key=lambda c: (c.bde_kcalmol
                                  if c.bde_kcalmol is not None
                                  else float("inf")))
    for rank, c in enumerate(updated, 1):
        c.rank = rank

    result.candidates = updated
    result.method = f"{functional}/{basis} BDE (top-{top_n}) + GFN2-xTB BDE (rest)"
    result.wall_seconds = (result.wall_seconds or 0.0) + (time.time() - t0)
    return result


def _dft_h_atom_energy(*, functional: str = "B3LYP",
                      basis: str = "def2-svp") -> Optional[float]:
    """DFT energy of an isolated H atom (doublet)."""
    try:
        from pyscf import gto, dft as pyscf_dft
    except ImportError:
        return None
    try:
        mol = gto.M(atom=[("H", (0.0, 0.0, 0.0))],
                    basis=basis, charge=0, spin=1,
                    max_memory=2000, verbose=0)
        mf = pyscf_dft.UKS(mol)
        mf.xc = functional.lower()
        mf.max_cycle = 60
        mf.conv_tol = 1e-8
        e = mf.kernel()
        if e is None or not getattr(mf, "converged", False):
            return None
        return float(e)
    except Exception as e:
        logger.debug("H-atom DFT failed: %s", e)
        return None


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles")
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--recompute-h", action="store_true",
                    help="recompute E(H·) for this run "
                         "(slower, calibrates to current xtb build)")
    ap.add_argument("--dft-verify", type=int, default=0,
                    metavar="N",
                    help="after xTB ranks top-N, re-rank those N "
                         "with DFT B3LYP/def2-SVP for paper-grade "
                         "accuracy (default 0 = skip DFT).")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    r = predict_cyp_som_bde(
        args.smiles, top_k=args.top_k,
        recompute_h_ref=args.recompute_h)

    if args.dft_verify > 0 and r.ok:
        print(f"[som] xTB ranking complete; DFT-verifying top-"
              f"{args.dft_verify} …")
        r = verify_som_dft(r, top_n=args.dft_verify)

    if args.json:
        import json
        print(json.dumps(r.as_dict(), indent=2, default=str))
    else:
        print(r.summary())
    sys.exit(0 if r.ok else 1)
