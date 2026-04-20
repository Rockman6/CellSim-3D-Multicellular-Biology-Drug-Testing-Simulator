"""PySCF DFT single-point driver — rigorous QM complement to xTB.

xTB GFN2 is fast (~0.3 s per drug-sized molecule) but semi-empirical
— its absolute energies carry a ~20 kcal/mol systematic offset vs
experiment. When a biologist needs a second opinion on a reactive-
site prediction or a BDE, rerunning the top-3 candidates at DFT
level tightens the answer.

Default: B3LYP/def2-SVP — the canonical "is this energy roughly
right?" DFT recipe. ~10–30 s for a small drug (< 40 heavy atoms)
on a CPU. Users can escalate to B3LYP/def2-TZVP or ωB97X-D for
publication-grade numbers.

Non-AI. PySCF implements the Kohn-Sham equations from published
parameterisations; no learned surrogate anywhere.

Usage:
    from src.quantum import dft_single_point
    r = dft_single_point("CC(=O)OC1=CC=CC=C1C(=O)O")
    print(r.summary())
"""

from __future__ import annotations

import logging
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

_HARTREE_TO_EV = 27.211386245988
_HARTREE_TO_KCAL = 627.509474


@dataclass
class DftResult:
    """Envelope for a PySCF DFT single-point calculation."""

    smiles: str
    ok: bool
    reason: str = ""
    method: str = "DFT"
    functional: str = "B3LYP"
    basis: str = "def2-svp"
    charge: int = 0
    multiplicity: int = 1

    total_energy_Hartree: Optional[float] = None
    total_energy_eV: Optional[float] = None
    homo_eV: Optional[float] = None
    lumo_eV: Optional[float] = None
    homo_lumo_gap_eV: Optional[float] = None
    dipole_Debye: Optional[float] = None
    n_atoms: int = 0

    elements: list = field(default_factory=list)
    positions_A: list = field(default_factory=list)

    tool_versions: dict = field(default_factory=dict)
    wall_seconds: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] DFT {self.smiles}  {self.reason}"

        def _f(x, fmt=".2f", u=""):
            return f"{x:{fmt}}{u}" if x is not None else "n/a"

        return (f"[OK]   DFT {self.functional}/{self.basis}  "
                f"{self.smiles}  atoms={self.n_atoms}  "
                f"E={_f(self.total_energy_eV, '.2f', ' eV')}  "
                f"HOMO={_f(self.homo_eV)}  LUMO={_f(self.lumo_eV)}  "
                f"gap={_f(self.homo_lumo_gap_eV, '.2f', ' eV')}  "
                f"µ={_f(self.dipole_Debye, '.2f', ' D')}  "
                f"({self.wall_seconds:.1f} s)")


def _tool_versions() -> dict:
    versions = {}
    try:
        import pyscf
        versions["pyscf"] = getattr(pyscf, "__version__", "?")
    except Exception:
        pass
    try:
        import rdkit
        versions["rdkit"] = rdkit.__version__
    except Exception:
        pass
    return versions


def _mol_xyz_3d(smiles: str, random_seed: int = 1
                 ) -> Optional[tuple]:
    """Build a 3D-embedded RDKit mol → (elements, xyz_Å)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    xyz = [[conf.GetAtomPosition(i).x,
            conf.GetAtomPosition(i).y,
            conf.GetAtomPosition(i).z]
           for i in range(mol.GetNumAtoms())]
    return elements, xyz


def dft_single_point(
    smiles: str,
    *,
    functional: str = "B3LYP",
    basis: str = "def2-svp",
    charge: int = 0,
    multiplicity: int = 1,
    random_seed: int = 1,
    max_memory_mb: int = 4000,
    cache: Optional[Any] = None,
) -> DftResult:
    """Run one PySCF DFT single-point on a SMILES.

    Returns `DftResult` with total energy, HOMO/LUMO, gap, dipole.
    Never raises on recoverable failure; `ok=False` with reason.
    """
    versions = _tool_versions()
    result = DftResult(
        smiles=smiles, ok=False,
        functional=functional, basis=basis,
        charge=charge, multiplicity=multiplicity,
        tool_versions=versions)

    # -------- cache lookup ------------------------------------
    cache_key = None
    cache_method = None
    if cache is not None:
        try:
            from src.cache.hashing import compound_hash, method_key
            cache_key = compound_hash(smiles)
            cache_method = method_key(
                "dft.single_point",
                versions.get("pyscf", "unknown"),
                {"functional": functional, "basis": basis,
                 "charge": charge, "mult": multiplicity,
                 "seed": random_seed})
            if cache_key is not None:
                hit = cache.get(cache_key, cache_method)
                if hit is not None:
                    logger.debug("dft cache HIT")
                    return DftResult(**{
                        k: v for k, v in hit.value.items()
                        if k in DftResult.__dataclass_fields__
                    })
        except Exception as e:
            logger.debug("dft cache lookup failed: %s", e)
            cache_key = None

    try:
        from pyscf import gto, dft
    except ImportError as e:
        result.reason = f"pyscf not installed: {e}"
        return result

    # Build geometry.
    geom = _mol_xyz_3d(smiles, random_seed=random_seed)
    if geom is None:
        result.reason = "RDKit embed failed for SMILES"
        return result
    elements, xyz = geom
    result.elements = elements
    result.positions_A = xyz
    result.n_atoms = len(elements)

    # Run PySCF.
    t0 = time.time()
    try:
        mol = gto.M(
            atom=[(sym, (x, y, z)) for sym, (x, y, z) in zip(elements, xyz)],
            basis=basis,
            charge=charge,
            spin=multiplicity - 1,    # 2S (unpaired electrons)
            max_memory=max_memory_mb,
            verbose=0,
        )
    except Exception as e:
        result.reason = f"PySCF geometry build failed: {str(e)[:120]}"
        result.wall_seconds = time.time() - t0
        return result

    try:
        if multiplicity == 1:
            mf = dft.RKS(mol)
        else:
            mf = dft.UKS(mol)
        mf.xc = functional.lower()
        mf.max_cycle = 60
        mf.conv_tol = 1e-8
        e_tot = mf.kernel()
    except Exception as e:
        result.reason = f"PySCF SCF failed: {str(e)[:120]}"
        result.wall_seconds = time.time() - t0
        return result

    if e_tot is None or not getattr(mf, "converged", False):
        result.reason = "SCF did not converge"
        result.wall_seconds = time.time() - t0
        return result

    result.total_energy_Hartree = float(e_tot)
    result.total_energy_eV = float(e_tot) * _HARTREE_TO_EV

    # Frontier orbitals (occupied highest / virtual lowest).
    try:
        mo_energy = mf.mo_energy
        mo_occ = mf.mo_occ
        # For restricted (RKS): 1D array. For unrestricted (UKS): tuple.
        if isinstance(mo_energy, tuple) or (hasattr(mo_energy, "ndim")
                                             and mo_energy.ndim == 2):
            # take the alpha channel for a simple HOMO/LUMO guess.
            me = mo_energy[0]
            mo = mo_occ[0]
        else:
            me = mo_energy
            mo = mo_occ
        occ_mask = mo > 0.5
        vir_mask = mo < 0.5
        if occ_mask.any() and vir_mask.any():
            homo = float(me[occ_mask].max()) * _HARTREE_TO_EV
            lumo = float(me[vir_mask].min()) * _HARTREE_TO_EV
            result.homo_eV = homo
            result.lumo_eV = lumo
            result.homo_lumo_gap_eV = lumo - homo
    except Exception as e:
        logger.debug("HOMO/LUMO extraction failed: %s", e)

    # Dipole moment (Debye).
    try:
        dip = mf.dip_moment(unit="Debye", verbose=0)
        import numpy as np
        result.dipole_Debye = float(np.linalg.norm(dip))
    except Exception as e:
        logger.debug("dipole extraction failed: %s", e)

    result.wall_seconds = time.time() - t0
    result.ok = True

    # cache put
    if cache is not None and cache_key is not None:
        try:
            cache.put(cache_key, cache_method, result.as_dict())
        except Exception as e:
            logger.debug("dft cache put failed: %s", e)

    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles")
    ap.add_argument("--functional", default="B3LYP")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    r = dft_single_point(
        args.smiles, functional=args.functional, basis=args.basis,
        charge=args.charge, multiplicity=args.mult,
        random_seed=args.seed)
    if args.json:
        import json
        print(json.dumps(r.as_dict(), indent=2, default=str))
    else:
        print(r.summary())
    sys.exit(0 if r.ok else 1)
