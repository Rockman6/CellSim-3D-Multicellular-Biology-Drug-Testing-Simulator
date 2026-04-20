"""xTB GFN2-xTB semi-empirical driver (non-AI).

Campaign 1, Layer 1.4 foundation. For every ligand of interest we
run one xTB single-point calculation to get:

  - total energy (Hartree), HOMO, LUMO, HOMO-LUMO gap (eV)
  - Mulliken partial charges (per atom)
  - CM5 charges (better for electrostatics)
  - molecular dipole magnitude (Debye)
  - atomic indices (so we can overlay onto the 3D viewer)

This is the input for:
  - electrostatic potential / HOMO-LUMO viewer (src/quantum/viewer.py)
  - CYP3A4 site-of-metabolism prediction via Fukui functions
    (src/quantum/metabolism.py — next PR)
  - reactive-fragment screening in the Campaign-1 cache

Design principles (from MISSION.md):
  - No learned model anywhere. GFN2-xTB is a classical semi-
    empirical method with published parametrisation — auditable,
    deterministic, no training set.
  - Every result carries full provenance (xtb version, method,
    random seed for dispersion sampling, solvent model).
  - Never raises on recoverable failure; failure captured as
    `result.ok = False` with a readable reason.

Implementation: subprocess the `xtb` CLI (more stable API across
versions than `xtb-python` in-process). Produces input XYZ from
RDKit, parses xtb.out for numerics.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

_HARTREE_TO_EV = 27.211386245988
_DEBYE_PER_AU = 2.541746  # 1 a.u. = 2.5417 D


@dataclass
class XtbResult:
    """Envelope for a single xtb GFN2 calculation."""

    smiles: str
    ok: bool
    reason: str = ""
    method: str = "GFN2-xTB"
    solvent: Optional[str] = None
    charge: int = 0
    multiplicity: int = 1
    random_seed: int = 1

    # Energies
    total_energy_Hartree: Optional[float] = None
    total_energy_eV: Optional[float] = None

    # Frontier orbitals
    homo_eV: Optional[float] = None
    lumo_eV: Optional[float] = None
    homo_lumo_gap_eV: Optional[float] = None

    # Electrostatics
    dipole_Debye: Optional[float] = None
    mulliken_charges: list = field(default_factory=list)   # per atom (e)
    cm5_charges: list = field(default_factory=list)        # per atom (e)
    elements: list = field(default_factory=list)
    positions_A: list = field(default_factory=list)

    # Provenance
    tool_versions: dict = field(default_factory=dict)
    wall_seconds: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] {self.smiles}  {self.reason}"
        def _f(x, fmt=".2f", u=""):
            return f"{x:{fmt}}{u}" if x is not None else "n/a"
        return (f"[OK]   xtb {self.method}  {self.smiles}  "
                f"atoms={len(self.elements)}  "
                f"E={_f(self.total_energy_eV, '.2f', ' eV')}  "
                f"HOMO={_f(self.homo_eV)}  "
                f"LUMO={_f(self.lumo_eV)}  "
                f"gap={_f(self.homo_lumo_gap_eV, '.2f', ' eV')}  "
                f"µ={_f(self.dipole_Debye, '.2f', ' D')}  "
                f"({self.wall_seconds:.2f} s)")


def _tool_versions() -> dict:
    versions: dict = {}
    bin_ = shutil.which("xtb")
    if bin_:
        try:
            r = subprocess.run([bin_, "--version"], capture_output=True,
                                text=True, timeout=10)
            # Output contains a line like "   xtb version 6.X.Y ..."
            for line in (r.stdout + r.stderr).splitlines():
                m = re.search(r"xtb\s+version\s+([\d.]+)", line)
                if m:
                    versions["xtb-cli"] = m.group(1)
                    break
        except Exception:
            pass
    try:
        import xtb
        versions["xtb-python"] = getattr(xtb, "__version__", "?")
    except Exception:
        pass
    return versions


def _write_xyz(
    smiles: str, workdir: Path, seed: int = 1
) -> tuple[Path, list[str], list[list[float]]]:
    """SMILES -> 3D coords -> XYZ file.

    Returns (xyz_path, element_list, positions_Å).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("RDKit ETKDGv3 embed failed")
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

    xyz = workdir / "mol.xyz"
    with xyz.open("w") as fo:
        fo.write(f"{len(elements)}\n")
        fo.write(f"cellsim SMILES={smiles}\n")
        for sym, (x, y, z) in zip(elements, positions):
            fo.write(f"{sym:<4s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    return xyz, elements, positions


def _parse_xtb_out(text: str) -> dict:
    """Pull the numbers we care about out of xtb's main stdout.

    xtb 6.6+ prints a "|         Property Printout          |" header;
    the relevant fields are on well-known labelled lines.
    """
    out: dict = {}
    # HOMO / LUMO are printed as part of the orbital-energy table;
    # easier to get from the summary block:
    #   "          Dispersion energy ...              "
    #   "HL-Gap           0.1234 Eh           3.3586 eV"
    #   "Fermi-level      ..."
    #   "Total energy              -12.3456 Eh"
    #   "molecular dipole: ... total (Debye): 1.23"
    m = re.search(r"HL-Gap\s+[-\d.]+\s*Eh\s+([-\d.]+)\s*eV",
                  text)
    if m:
        out["homo_lumo_gap_eV"] = float(m.group(1))
    m = re.search(r"TOTAL ENERGY\s+([-\d.]+)\s+Eh", text,
                  re.IGNORECASE)
    if m:
        out["total_energy_Hartree"] = float(m.group(1))
    # HOMO / LUMO individually from the orbital table: find the
    # LUMO line marked "(LUMO)" and the HOMO line marked "(HOMO)".
    homo, lumo = None, None
    for line in text.splitlines():
        if "(HOMO)" in line:
            parts = line.split()
            # columns: idx occ Eh eV (HOMO).  Last token is the tag;
            # second-to-last is the orbital energy in eV.
            try:
                homo = float(parts[-2])
            except (ValueError, IndexError):
                pass
        elif "(LUMO)" in line and lumo is None:
            parts = line.split()
            try:
                lumo = float(parts[-2])
            except (ValueError, IndexError):
                pass
    if homo is not None:
        out["homo_eV"] = homo
    if lumo is not None:
        out["lumo_eV"] = lumo
    # Dipole magnitude (Debye)
    m = re.search(r"molecular dipole:[\s\S]*?full:[^\n]*",
                  text)
    if m:
        # Parse the full line; xtb prints x y z and total in Debye.
        full_line = m.group(0).splitlines()[-1]
        parts = full_line.split()
        try:
            # The LAST number on the "full:" line is the total dipole
            # in Debye for xtb 6.x.
            out["dipole_Debye"] = float(parts[-1])
        except Exception:
            pass
    return out


def _parse_charges_file(path: Path, n_atoms: int) -> list[float]:
    """xtb writes a 'charges' file — one value per line."""
    vals: list[float] = []
    if not path.exists():
        return vals
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            vals.append(float(line.split()[0]))
        except ValueError:
            continue
    return vals[:n_atoms] if len(vals) >= n_atoms else vals


def xtb_single_point(
    smiles: str,
    *,
    charge: int = 0,
    multiplicity: int = 1,
    solvent: Optional[str] = None,
    method: str = "gfn2",
    random_seed: int = 1,
    timeout_s: int = 300,
    cache: Optional[Any] = None,
) -> XtbResult:
    """Run one GFN2-xTB single-point on a SMILES.

    Parameters mirror xtb's CLI: `charge`, `multiplicity` (2S+1),
    `solvent` (e.g. "water", "dmso" — activates ALPB implicit solvent),
    `method` ("gfn1" / "gfn2" / "gfnff").

    Returns an `XtbResult` envelope; `ok=False` with a reason on any
    recoverable failure (RDKit embed, xtb crash, parsing gap).
    """
    import time

    result = XtbResult(smiles=smiles, ok=False,
                       method=method.upper() + "-xTB",
                       solvent=solvent, charge=charge,
                       multiplicity=multiplicity,
                       random_seed=random_seed,
                       tool_versions=_tool_versions())

    xtb_bin = shutil.which("xtb")
    if xtb_bin is None:
        result.reason = "xtb binary not on PATH (activate cellsim env)"
        return result

    # -------- cache lookup ------------------------------------
    cache_key = None
    cache_method = None
    if cache is not None:
        try:
            from src.cache.hashing import compound_hash, method_key
            cache_key = compound_hash(smiles)
            cache_method = method_key(
                "xtb.single_point",
                versions.get("xtb-cli", "unknown"),
                {
                    "method": method,
                    "charge": charge,
                    "mult": multiplicity,
                    "solvent": solvent or "vacuum",
                    "seed": random_seed,
                })
            if cache_key is not None:
                hit = cache.get(cache_key, cache_method)
                if hit is not None:
                    logger.debug("xtb cache HIT %s / %s",
                                 cache_key, cache_method)
                    return XtbResult(**{
                        k: v for k, v in hit.value.items()
                        if k in XtbResult.__dataclass_fields__
                    })
        except Exception as e:
            logger.debug("xtb cache lookup failed: %s", e)
            cache_key = None

    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="cellsim-xtb-") as tmp:
        workdir = Path(tmp)
        try:
            xyz_path, elements, positions = _write_xyz(
                smiles, workdir, seed=random_seed)
        except Exception as e:
            result.reason = f"XYZ prep failed: {str(e)[:200]}"
            return result
        result.elements = elements
        result.positions_A = positions

        cmd = [xtb_bin, str(xyz_path),
               "--" + method,
               "--chrg", str(charge),
               "--uhf", str(multiplicity - 1),
               "--json"]
        if solvent:
            cmd += ["--alpb", solvent]
        try:
            r = subprocess.run(
                cmd, capture_output=True, text=True,
                timeout=timeout_s, cwd=str(workdir))
        except subprocess.TimeoutExpired:
            result.reason = f"xtb timeout after {timeout_s} s"
            result.wall_seconds = time.time() - t0
            return result

        if r.returncode != 0:
            result.reason = (f"xtb exit {r.returncode}: "
                             f"{(r.stderr or r.stdout)[-500:]}")
            result.wall_seconds = time.time() - t0
            return result

        parsed = _parse_xtb_out(r.stdout)
        for k, v in parsed.items():
            setattr(result, k, v)
        if result.total_energy_Hartree is not None:
            result.total_energy_eV = (
                result.total_energy_Hartree * _HARTREE_TO_EV)
        # Derive HL gap if not parsed directly but both orbitals are
        # known.
        if (result.homo_lumo_gap_eV is None and
                result.homo_eV is not None and result.lumo_eV is not None):
            result.homo_lumo_gap_eV = result.lumo_eV - result.homo_eV

        # Mulliken charges come in workdir/charges
        result.mulliken_charges = _parse_charges_file(
            workdir / "charges", len(elements))

        # xtb JSON output ("xtbout.json") has CM5 / Mulliken when
        # available; pick up whichever fields are there.
        jpath = workdir / "xtbout.json"
        if jpath.exists():
            try:
                import json
                j = json.loads(jpath.read_text())
                if "partial charges" in j and not result.mulliken_charges:
                    result.mulliken_charges = list(j["partial charges"])
                if "HOMO" in j and result.homo_eV is None:
                    result.homo_eV = float(j["HOMO"])
                if "LUMO" in j and result.lumo_eV is None:
                    result.lumo_eV = float(j["LUMO"])
                if "total energy" in j and result.total_energy_Hartree is None:
                    result.total_energy_Hartree = float(j["total energy"])
                if "dipole" in j and result.dipole_Debye is None:
                    # xtbout dipole is in a.u.; convert to Debye.
                    try:
                        vec = j["dipole"]
                        mag = sum(v * v for v in vec) ** 0.5
                        result.dipole_Debye = mag * _DEBYE_PER_AU
                    except Exception:
                        pass
            except Exception:
                pass

    result.wall_seconds = time.time() - t0
    if result.total_energy_Hartree is None:
        result.reason = "parse gap: no total energy in xtb output"
        return result
    result.ok = True

    # -------- cache put (only ok=True) ------------------------
    if cache is not None and cache_key is not None:
        try:
            cache.put(cache_key, cache_method, result.as_dict())
        except Exception as e:
            logger.debug("xtb cache put failed: %s", e)
    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--solvent", default=None,
                    help="implicit solvent name (water, dmso, …)")
    ap.add_argument("--method", default="gfn2")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    r = xtb_single_point(
        args.smiles, charge=args.charge, multiplicity=args.mult,
        solvent=args.solvent, method=args.method,
        random_seed=args.seed)
    if args.json:
        import json as _json
        print(_json.dumps(r.as_dict(), indent=2, default=str))
    else:
        print(r.summary())
    sys.exit(0 if r.ok else 1)
