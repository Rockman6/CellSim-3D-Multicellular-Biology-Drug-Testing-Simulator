"""SMILES -> parameterised OpenMM System.

Pipeline:
    SMILES
      -> RDKit parse + 3D embed + conformer optimise (UFF/MMFF)
      -> OpenFF Molecule
      -> OpenFF Sage force field assign parameters
      -> AM1-BCC partial charges
      -> openmm.System (minimizable)

Every result carries a `ParametrizeResult` envelope with method,
tool versions, and the success/failure reason — no silent fallbacks.
This is the Campaign-1 UQ-envelope policy enforced from day one.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional
import hashlib
import json
import logging
import traceback

logger = logging.getLogger(__name__)


@dataclass
class ParametrizeResult:
    """Envelope around a (possibly failed) parametrisation attempt."""

    smiles: str
    ok: bool
    reason: str = ""
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    formula: Optional[str] = None
    n_atoms: Optional[int] = None
    n_conformers: Optional[int] = None
    charge_method: Optional[str] = None
    ff_version: Optional[str] = None
    tool_versions: dict = field(default_factory=dict)
    positions_nm: Optional[list] = None
    partial_charges_e: Optional[list] = None
    elements: Optional[list] = None
    bonds: Optional[list] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def hash_key(self) -> str:
        """Content hash keyed by inchi_key + method — cache-ready."""
        if not self.inchi_key:
            raise ValueError("no inchi_key — parametrisation failed")
        method = f"{self.charge_method}|{self.ff_version}"
        payload = f"{self.inchi_key}|{method}".encode()
        return hashlib.sha256(payload).hexdigest()[:16]


def _tool_versions() -> dict:
    versions = {}
    try:
        import rdkit

        versions["rdkit"] = rdkit.__version__
    except Exception:
        pass
    try:
        import openff.toolkit

        versions["openff-toolkit"] = openff.toolkit.__version__
    except Exception:
        pass
    try:
        import openmm

        versions["openmm"] = openmm.__version__
    except Exception:
        pass
    return versions


def parametrize_smiles(
    smiles: str,
    *,
    charge_method: str = "am1bcc",
    ff_name: str = "openff-2.1.0.offxml",
    random_seed: int = 1,
    cache: Optional["object"] = None,
) -> ParametrizeResult:
    """Return a ParametrizeResult for the given SMILES.

    If `cache` (an `src.cache.Cache` instance) is supplied, an
    identical-physics cache hit short-circuits the expensive
    AM1-BCC + OpenFF path (~20 s → ~ms). Only `ok=True` results
    are cached, so transient failures don't poison future runs.
    Cache key is content-addressed: `(compound_hash(smiles),
    method=parametrize|ff|charge|seed)`.

    Never raises for recoverable failures (unparsable SMILES, missing
    parameters, etc.) — the failure is captured in `result.ok = False`
    with a human-readable `reason`. Upstream callers treat the batch
    like a statistic, not an exception stream.
    """

    versions = _tool_versions()

    # -------- cache lookup ------------------------------------
    cache_key = None
    cache_method = None
    if cache is not None:
        try:
            from src.cache.hashing import compound_hash, method_key
            cache_key = compound_hash(smiles)
            cache_method = method_key(
                "parametrize", ff_name,
                {"charge": charge_method, "seed": random_seed})
            if cache_key is not None:
                hit = cache.get(cache_key, cache_method)
                if hit is not None:
                    logger.debug("cache HIT for %s / %s",
                                 cache_key, cache_method)
                    return ParametrizeResult(**{
                        k: v for k, v in hit.value.items()
                        if k in ParametrizeResult.__dataclass_fields__
                    })
        except Exception as e:
            logger.debug("cache lookup failed: %s", e)
            cache_key = None  # fall through to compute path

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:
        return ParametrizeResult(
            smiles=smiles,
            ok=False,
            reason=f"rdkit import failed: {e}",
            tool_versions=versions,
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ParametrizeResult(
            smiles=smiles,
            ok=False,
            reason="RDKit could not parse SMILES",
            tool_versions=versions,
        )

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        return ParametrizeResult(
            smiles=smiles,
            ok=False,
            reason="RDKit ETKDGv3 embed failed",
            tool_versions=versions,
        )
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e:
            return ParametrizeResult(
                smiles=smiles,
                ok=False,
                reason=f"conformer optimise failed: {e}",
                tool_versions=versions,
            )

    inchi = Chem.MolToInchi(mol)
    inchi_key = Chem.InchiToInchiKey(inchi) if inchi else None
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

    # Always capture RDKit geometry so the viewer works even when
    # downstream steps (OpenFF import, charge assignment, FF assignment)
    # later fail.
    conf = mol.GetConformer()
    rdkit_positions = [
        [conf.GetAtomPosition(i).x * 0.1,
         conf.GetAtomPosition(i).y * 0.1,
         conf.GetAtomPosition(i).z * 0.1]
        for i in range(mol.GetNumAtoms())
    ]
    rdkit_elements = [a.GetSymbol() for a in mol.GetAtoms()]
    rdkit_bonds = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
         b.GetBondTypeAsDouble())
        for b in mol.GetBonds()
    ]
    rdkit_n_atoms = mol.GetNumAtoms()
    rdkit_n_conformers = mol.GetNumConformers()

    def _rdkit_fallback(reason: str, **extra) -> ParametrizeResult:
        return ParametrizeResult(
            smiles=smiles,
            ok=False,
            reason=reason,
            inchi=inchi,
            inchi_key=inchi_key,
            formula=formula,
            n_atoms=rdkit_n_atoms,
            n_conformers=rdkit_n_conformers,
            positions_nm=rdkit_positions,
            elements=rdkit_elements,
            bonds=rdkit_bonds,
            tool_versions=versions,
            **extra,
        )

    try:
        from openff.toolkit import ForceField, Molecule
    except ImportError as e:
        return _rdkit_fallback(f"openff-toolkit not available: {e}")

    try:
        off_mol = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    except Exception as e:
        return _rdkit_fallback(f"OpenFF Molecule.from_rdkit failed: {e}")

    try:
        off_mol.assign_partial_charges(partial_charge_method=charge_method)
    except Exception as e:
        # Truncate the very long OpenFF diagnostic to the first line
        # and surface the #1 user-level cause: running outside an
        # activated conda env so AMBERHOME is unset and the
        # AmberToolsToolkitWrapper never registers.
        reason = str(e).splitlines()[0]
        import os

        if (charge_method.startswith("am1") and
                not os.environ.get("AMBERHOME")):
            reason += (
                " [likely cause: AMBERHOME is unset — run under "
                "`conda activate cellsim` or `mamba run -n cellsim`]")
        return _rdkit_fallback(
            f"{charge_method} charge assignment failed: {reason}")

    try:
        ff = ForceField(ff_name)
        system = ff.create_openmm_system(off_mol.to_topology())
    except Exception as e:
        return _rdkit_fallback(
            f"OpenFF parametrisation failed ({ff_name}): {e}",
            charge_method=charge_method)

    # Extract geometry + charges for cache + viewer.
    conf = off_mol.conformers[0]
    try:
        from openff.units import unit as offunit

        positions = conf.m_as(offunit.nanometer).tolist()
    except Exception:
        positions = [list(p) for p in conf]

    charges_e = [float(q.magnitude) if hasattr(q, "magnitude")
                 else float(q)
                 for q in off_mol.partial_charges]
    elements = [a.symbol for a in off_mol.atoms]
    bonds = [
        (b.atom1_index, b.atom2_index, float(b.bond_order))
        for b in off_mol.bonds
    ]

    # Sanity: the generated System should at least construct.
    _ = system.getNumParticles()

    result = ParametrizeResult(
        smiles=smiles,
        ok=True,
        inchi=inchi,
        inchi_key=inchi_key,
        formula=formula,
        n_atoms=off_mol.n_atoms,
        n_conformers=off_mol.n_conformers,
        charge_method=charge_method,
        ff_version=ff_name,
        tool_versions=versions,
        positions_nm=positions,
        partial_charges_e=charges_e,
        elements=elements,
        bonds=bonds,
    )

    # -------- cache put (only on ok=True) ---------------------
    if cache is not None and cache_key is not None:
        try:
            cache.put(cache_key, cache_method, result.as_dict())
        except Exception as e:
            logger.debug("cache put failed: %s", e)
    return result


def dump_json(result: ParametrizeResult) -> str:
    return json.dumps(result.as_dict(), indent=2, default=str)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Parameterise a SMILES into an OpenMM system.")
    p.add_argument("smiles", help="Input SMILES string")
    p.add_argument("--charge", default="am1bcc",
                   help="partial charge method (am1bcc, am1-mulliken, …)")
    p.add_argument("--ff", default="openff-2.1.0.offxml",
                   help="OpenFF force-field file")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--json", action="store_true",
                   help="emit the full result dict as JSON")
    args = p.parse_args()

    result = parametrize_smiles(
        args.smiles,
        charge_method=args.charge,
        ff_name=args.ff,
        random_seed=args.seed,
    )
    if args.json:
        print(dump_json(result))
    else:
        if result.ok:
            print(f"OK  {result.smiles}  formula={result.formula} "
                  f"atoms={result.n_atoms} charges={result.charge_method} "
                  f"ff={result.ff_version}  hash={result.hash_key()}")
        else:
            print(f"FAIL  {result.smiles}  reason: {result.reason}")
