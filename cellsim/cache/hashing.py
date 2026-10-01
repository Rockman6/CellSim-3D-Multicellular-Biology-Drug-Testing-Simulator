"""Stable content hashes for the physics-prior cache.

Goal: two inputs that *physically mean the same thing* must produce
the same hash, regardless of cosmetic differences (SMILES atom
ordering, PDB header metadata, …). Non-AI, deterministic,
dependency-light.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Optional


def compound_hash(smiles: str) -> Optional[str]:
    """Canonical InChI Key based hash for a SMILES.

    Uses RDKit to normalise the SMILES → canonical InChI → InChI
    Key → first 16 hex chars of SHA-256 over the key. The
    InChI Key already canonicalises stereo + tautomer information
    as far as InChI understands it, so compounds that differ only
    in SMILES atom ordering get identical hashes.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    inchi = Chem.MolToInchi(Chem.AddHs(mol))
    if not inchi:
        return None
    key = Chem.InchiToInchiKey(inchi)
    if not key:
        return None
    return hashlib.sha256(key.encode()).hexdigest()[:16]


_PDB_COORD_RE = re.compile(r"^(ATOM  |HETATM)")


def receptor_hash(pdb_path: str | Path,
                   include_hetatm: bool = False) -> Optional[str]:
    """Hash a receptor PDB by its physical atoms only.

    Strategy: extract `(residue_name, residue_id, atom_name, chain,
    x, y, z)` for every ATOM row (optionally HETATM), sort by
    (chain, residue_id, atom_name), canonicalise floats to 3
    decimal places, hash the joined record stream.

    Result is stable across header / remark / solvent differences.
    """
    p = Path(pdb_path)
    if not p.exists():
        return None
    records: list[str] = []
    for line in p.read_text().splitlines():
        if not _PDB_COORD_RE.match(line):
            continue
        if line.startswith("HETATM") and not include_hetatm:
            continue
        # skip water / common ions by convention — these can vary
        # between preps without changing the receptor identity.
        resname = line[17:20].strip()
        if resname in ("HOH", "WAT", "TIP3", "NA", "CL",
                       "NA+", "CL-", "ZN", "MG", "CA", "MN", "FE"):
            continue
        try:
            atom_name = line[12:16].strip()
            chain = line[21]
            res_id = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        records.append(
            f"{chain}|{res_id:>5d}|{resname:<4s}|{atom_name:<4s}|"
            f"{x:.3f}|{y:.3f}|{z:.3f}")
    if not records:
        return None
    records.sort()
    payload = "\n".join(records).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def method_key(name: str, version: str, extra: Optional[dict] = None) -> str:
    """Canonical method identifier. `name`+`version` plus sorted
    key=value pairs from `extra` (e.g. exhaustiveness, seed)."""
    parts = [name, version]
    if extra:
        for k in sorted(extra.keys()):
            parts.append(f"{k}={extra[k]}")
    return "|".join(parts)


def full_key(content_hash: str, method: str) -> str:
    return f"{content_hash}::{method}"
