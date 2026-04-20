"""AutoDock Vina docking driver (non-AI).

Design principles (per MISSION.md):
  - Primary evidence is Vina's empirical scoring function; auditable,
    no black-box.
  - Every result carries full provenance: Vina version, seed,
    exhaustiveness, search box, receptor + ligand hashes.
  - Biologist-friendly units: ΔG reported in kcal/mol (default) AND
    kJ/mol (conversion factor 4.184), K_d implied via
    K_d = exp(ΔG / RT).
  - Deterministic: `--seed` is passed through; the same (receptor,
    ligand, box, seed, exhaustiveness) always yields the same poses.
  - Never raises on recoverable failure; failure captured in
    `result.ok = False` with a readable `reason`.

Pipeline:
    Layer-1.2 receptor (minimised) -> PDB temp file -> Meeko -> PDBQT
    Layer-1.1 ligand (parametrised) -> SDF temp file -> Meeko -> PDBQT
    Vina v1.2.x CLI via subprocess (stable API; Python bindings are
    optional here).
    Parse out_ligand.pdbqt to extract poses + scores.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional
import hashlib
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile

logger = logging.getLogger(__name__)

# Gas constant in kcal/(mol·K). At T = 298.15 K, RT ≈ 0.5925 kcal/mol.
_R_KCAL = 1.987204e-3
_T_K = 298.15
_KCAL_TO_KJ = 4.184


def _tool_versions() -> dict:
    versions: dict = {}
    try:
        import vina
        versions["vina-python"] = getattr(vina, "__version__", "?")
    except Exception:
        pass
    try:
        import meeko
        versions["meeko"] = getattr(meeko, "__version__", "?")
    except Exception:
        pass
    # Probe the vina binary too — the Python package version is not
    # always the same as the executable.
    vbin = shutil.which("vina")
    if vbin:
        try:
            out = subprocess.run(
                [vbin, "--version"], capture_output=True, text=True,
                timeout=10)
            versions["vina-cli"] = (out.stdout or out.stderr).strip()
        except Exception:
            pass
    return versions


@dataclass
class DockingPose:
    """One ranked docking pose."""

    mode: int                       # 1-indexed pose rank
    affinity_kcalmol: float         # Vina ΔG
    affinity_kJmol: float           # same in SI
    kd_implied_nM: float            # exp(ΔG / RT) in nanomolar
    rmsd_lb_A: float = 0.0          # Vina "rmsd lb" vs top pose
    rmsd_ub_A: float = 0.0          # Vina "rmsd ub" vs top pose
    rmsd_vs_reference_A: Optional[float] = None  # filled if caller
                                                 # supplied a reference
    # PoseBusters flags (filled by src/dock/validity.attach_posebusters).
    # The split matters for biologists:
    #   posebusters_pocket_ok — pose is in-pocket with no clashes
    #                           against protein / ions / waters.
    #                           Trustworthy for triage ranking.
    #   posebusters_geometry_ok — bond lengths, angles, chirality,
    #                           internal clash, internal energy all
    #                           pass. Required for downstream FEP.
    #                           Vina poses routinely fail this one
    #                           because Vina uses approximate sterics;
    #                           fix = run a short MD minimisation
    #                           (future Layer 1.3 PR).
    #   posebusters_ok        — every PB test passes, including the
    #                           strict rmsd_≤_2Å. Strictest gate.
    posebusters_ok: Optional[bool] = None
    posebusters_pocket_ok: Optional[bool] = None
    posebusters_geometry_ok: Optional[bool] = None
    posebusters_flags: Optional[dict] = None
    # Coordinates (Å) — per-atom positions in the same atom ordering
    # as the prep step produced.
    positions_A: list = field(default_factory=list)
    elements: list = field(default_factory=list)

    def biologist_summary(self) -> str:
        """One-line human summary in biologist units."""
        rmsd_tag = ""
        if self.rmsd_vs_reference_A is not None:
            if self.rmsd_vs_reference_A < 2.0:
                rmsd_tag = f"  ✓ crystal-RMSD {self.rmsd_vs_reference_A:.2f} Å"
            elif self.rmsd_vs_reference_A < 3.0:
                rmsd_tag = f"  ~ crystal-RMSD {self.rmsd_vs_reference_A:.2f} Å"
            else:
                rmsd_tag = f"  ✗ crystal-RMSD {self.rmsd_vs_reference_A:.2f} Å"
        pb = ""
        # Prefer the biologist-relevant pocket flag for the summary
        # line (geometry noise from Vina's approximate sterics is
        # not an interesting triage signal).
        if self.posebusters_pocket_ok is True:
            pb = "  pocket:ok"
        elif self.posebusters_pocket_ok is False:
            pb = "  pocket:fail"
        elif self.posebusters_ok is True:
            pb = "  PB:ok"
        elif self.posebusters_ok is False:
            pb = "  PB:fail"
        return (f"#{self.mode:>2d}  ΔG = {self.affinity_kcalmol:>6.2f} kcal/mol  "
                f"({self.affinity_kJmol:>6.1f} kJ/mol)  "
                f"K_d ≈ {self.kd_implied_nM:>8.2g} nM{rmsd_tag}{pb}")


@dataclass
class DockingResult:
    """Envelope for a (possibly failed) docking run."""

    receptor_source: str
    ligand_smiles: str
    ok: bool
    reason: str = ""

    # Identity
    receptor_hash: Optional[str] = None
    ligand_inchi_key: Optional[str] = None
    ligand_formula: Optional[str] = None

    # Vina inputs (full provenance)
    center_A: Optional[tuple] = None        # (x, y, z) in Å
    box_size_A: Optional[tuple] = None       # (dx, dy, dz) in Å
    exhaustiveness: Optional[int] = None
    num_modes: Optional[int] = None
    seed: Optional[int] = None

    # Results
    poses: list = field(default_factory=list)       # list[DockingPose]
    best_kcalmol: Optional[float] = None            # = poses[0].affinity_kcalmol
    best_rmsd_vs_reference_A: Optional[float] = None

    # Metadata
    tool_versions: dict = field(default_factory=dict)
    wall_seconds: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] {self.ligand_formula or self.ligand_smiles}  {self.reason}"
        lines = [
            f"[OK]   dock {self.ligand_formula or self.ligand_smiles}  "
            f"into {Path(self.receptor_source).name}  "
            f"(exh={self.exhaustiveness}, modes={self.num_modes}, "
            f"seed={self.seed}, {self.wall_seconds:.1f} s)",
            f"  best ΔG = {self.best_kcalmol:.2f} kcal/mol  "
            f"({self.best_kcalmol * _KCAL_TO_KJ:.1f} kJ/mol)",
        ]
        for p in self.poses[:5]:
            lines.append(f"  " + p.biologist_summary())
        if len(self.poses) > 5:
            lines.append(f"  ... ({len(self.poses) - 5} more)")
        return "\n".join(lines)

    def cache_key(self) -> str:
        """Cache key per `src/cache/` schema: keyed on receptor +
        ligand + search-box + seed + exhaustiveness."""
        if not (self.receptor_hash and self.ligand_inchi_key):
            raise ValueError("cannot cache-key a failed docking result")
        parts = [self.receptor_hash, self.ligand_inchi_key,
                 str(self.center_A), str(self.box_size_A),
                 str(self.exhaustiveness), str(self.num_modes),
                 str(self.seed)]
        return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def _kd_from_dG(dG_kcalmol: float) -> float:
    """Return K_d in nanomolar, implied from ΔG in kcal/mol at 298 K.

    K_d = exp(ΔG / RT). Biologists usually think in nM, so convert
    from molar to nanomolar (× 1e9).
    """
    return math.exp(dG_kcalmol / (_R_KCAL * _T_K)) * 1e9


def _receptor_hash(pdb_text: str) -> str:
    return hashlib.sha256(pdb_text.encode()).hexdigest()[:16]


def _prep_receptor_pdbqt(pdb_path: Path, workdir: Path) -> Path:
    """Convert a minimised PDB to PDBQT via Meeko's CLI.

    The incoming PDB is expected to be pre-minimised (Layer 1.2
    pipeline); we strip solvent / ions / heteroatoms and run
    `mk_prepare_receptor.py` to produce the AutoDock atom types Vina
    expects.
    """
    out = workdir / "receptor.pdbqt"
    script = shutil.which("mk_prepare_receptor.py")
    if script is None:
        raise RuntimeError(
            "mk_prepare_receptor.py not found — ensure Meeko is installed "
            "in the active conda env")

    # Strip waters / ions / TER into a cleaned PDB first.
    cleaned = workdir / "receptor_clean.pdb"
    with pdb_path.open() as fi, cleaned.open("w") as fo:
        for line in fi:
            if line.startswith(("HEADER", "TITLE", "REMARK")):
                fo.write(line)
            elif line.startswith("ATOM"):
                resname = line[17:20].strip()
                if resname in ("HOH", "WAT", "TIP3", "NA", "CL",
                               "NA+", "CL-", "ZN", "MG", "CA"):
                    continue
                fo.write(line)
            elif line.startswith(("TER", "END")):
                fo.write(line)

    cmd = [script,
           "--read_pdb", str(cleaned),
           "-p", str(out),
           # Robust defaults for the zoo of X-ray PDBs biologists
           # throw at the tool: take the primary altloc of any
           # residue with alternates, and drop residues whose
           # chemistry doesn't match an AMBER template (e.g.
           # selenomethionine, phosphorylated residues, metals
           # pre-stripped above).
           "--default_altloc", "A",
           "-a"]
    try:
        r = subprocess.run(
            cmd, capture_output=True, text=True, timeout=180)
    except Exception as e:
        raise RuntimeError(f"mk_prepare_receptor.py crashed: {e}")
    if r.returncode != 0 or not out.exists():
        raise RuntimeError(
            f"mk_prepare_receptor.py failed (rc={r.returncode}): "
            f"{(r.stderr or r.stdout)[:500]}")
    return out


def _prep_ligand_pdbqt(smiles: str, workdir: Path, seed: int = 1) -> Path:
    """SMILES → 3D conformer → PDBQT via Meeko (0.7+ API)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

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

    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        raise RuntimeError("Meeko produced no ligand setups")
    # Meeko 0.7 returns (pdbqt_string, is_ok, error_message).
    pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        raise RuntimeError(f"Meeko PDBQT write failed: {err}")
    out = workdir / "ligand.pdbqt"
    out.write_text(pdbqt_str)
    return out


def _parse_vina_output(pdbqt_out: Path) -> list[DockingPose]:
    """Parse Vina's output PDBQT into DockingPose objects.

    Vina writes MODEL/ENDMDL blocks; the `REMARK VINA RESULT:` line
    of each model carries `affinity rmsd_lb rmsd_ub`.
    """
    poses: list[DockingPose] = []
    current_mode = 0
    current_score = None
    current_lb = 0.0
    current_ub = 0.0
    positions: list = []
    elements: list = []

    for line in pdbqt_out.read_text().splitlines():
        if line.startswith("MODEL"):
            current_mode = int(line.split()[-1])
            positions = []
            elements = []
            current_score = None
        elif line.startswith("REMARK VINA RESULT:"):
            parts = line.split()
            # REMARK VINA RESULT: <affinity> <rmsd_lb> <rmsd_ub>
            try:
                current_score = float(parts[3])
                current_lb = float(parts[4])
                current_ub = float(parts[5])
            except (IndexError, ValueError):
                current_score = None
        elif line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                # Element extraction from PDBQT. The canonical PDB
                # element column (76-78) is only 2 chars wide in
                # AutoDock-Vina output and gets chopped for two-
                # letter elements like Cl/Br. The atom-NAME column
                # [12:16] is 4 chars wide and reliably carries the
                # full element symbol (e.g. " Cl ", " C  ", " N  ")
                # — Meeko writes the periodic-table symbol here
                # rather than a raw atom name.
                elem = line[12:16].strip()
                if not elem:
                    elem = line[76:78].strip()
                # Strip trailing digits (some AutoDock variants
                # write " C1 " / " N2 " for numbered atoms).
                elem = "".join(c for c in elem if c.isalpha())
            except ValueError:
                continue
            positions.append([x, y, z])
            elements.append(elem.title())
        elif line.startswith("ENDMDL"):
            if current_score is not None and positions:
                poses.append(DockingPose(
                    mode=current_mode,
                    affinity_kcalmol=current_score,
                    affinity_kJmol=current_score * _KCAL_TO_KJ,
                    kd_implied_nM=_kd_from_dG(current_score),
                    rmsd_lb_A=current_lb,
                    rmsd_ub_A=current_ub,
                    positions_A=positions,
                    elements=elements,
                ))
    return poses


def dock_ligand(
    receptor_pdb: str | Path,
    ligand_smiles: str,
    *,
    center_A: tuple[float, float, float],
    box_size_A: tuple[float, float, float] = (22.0, 22.0, 22.0),
    exhaustiveness: int = 16,
    num_modes: int = 9,
    seed: int = 1,
    cpu: int = 0,
    timeout_s: int = 900,
    cache: Optional[Any] = None,
) -> DockingResult:
    """Dock a ligand SMILES into a prepared receptor PDB.

    `receptor_pdb` should be a pre-minimised PDB (e.g. from
    `src/md/protein.load_protein_pdb(...).minimised_positions_nm`
    dumped via OpenMM's PDBFile.writeFile).

    If `cache` is provided (an `src.cache.Cache` instance), an
    identical (ligand, receptor, search box, exhaustiveness, seed,
    num_modes) tuple is short-circuited from the store instead of
    re-running Vina (~5 s → ms). Only `ok=True` results are cached.

    Returns a `DockingResult` with up to `num_modes` poses sorted by
    Vina score ascending (best first).
    """
    import time

    versions = _tool_versions()
    receptor_pdb = Path(receptor_pdb)
    result = DockingResult(
        receptor_source=str(receptor_pdb),
        ligand_smiles=ligand_smiles,
        ok=False,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
        center_A=tuple(center_A),
        box_size_A=tuple(box_size_A),
        tool_versions=versions,
    )

    if not receptor_pdb.exists():
        result.reason = f"receptor not found: {receptor_pdb}"
        return result

    try:
        result.receptor_hash = _receptor_hash(receptor_pdb.read_text())
    except Exception as e:
        result.reason = f"receptor read failed: {e}"
        return result

    # Resolve InChI key for cache + provenance.
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(ligand_smiles)
        if mol is not None:
            result.ligand_inchi_key = Chem.InchiToInchiKey(
                Chem.MolToInchi(Chem.AddHs(mol)))
            result.ligand_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        pass

    # -------- cache lookup ------------------------------------
    # Key the cache entry on (ligand + receptor + search box + Vina
    # knobs + seed). Compound + receptor hashes are already
    # content-addressed so the key is stable across renames.
    cache_key = None
    cache_method = None
    if cache is not None:
        try:
            from src.cache.hashing import compound_hash, method_key
            lig_h = compound_hash(ligand_smiles)
            rec_h = result.receptor_hash
            if lig_h is not None and rec_h is not None:
                cache_key = f"{lig_h}+{rec_h}"
                # Truncate center/box floats to 3 dp so tiny
                # floating-point noise doesn't miss the cache.
                cache_method = method_key(
                    "vina.dock",
                    (versions.get("vina-cli")
                     or versions.get("vina-python") or "unknown"),
                    {
                        "cx": round(center_A[0], 3),
                        "cy": round(center_A[1], 3),
                        "cz": round(center_A[2], 3),
                        "bx": round(box_size_A[0], 3),
                        "by": round(box_size_A[1], 3),
                        "bz": round(box_size_A[2], 3),
                        "exh": exhaustiveness,
                        "modes": num_modes,
                        "seed": seed,
                    })
                hit = cache.get(cache_key, cache_method)
                if hit is not None:
                    logger.debug("vina cache HIT %s / %s",
                                 cache_key, cache_method)
                    try:
                        return _inflate_docking_result(hit.value)
                    except Exception as e:
                        logger.debug("vina cache inflate failed: %s", e)
                        # Fall through and recompute on corrupt rows.
        except Exception as e:
            logger.debug("vina cache lookup failed: %s", e)
            cache_key = None

    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="cellsim-dock-") as tmp:
        workdir = Path(tmp)

        try:
            rec_pdbqt = _prep_receptor_pdbqt(receptor_pdb, workdir)
        except Exception as e:
            result.reason = f"receptor prep failed: {e}"
            return result

        try:
            lig_pdbqt = _prep_ligand_pdbqt(ligand_smiles, workdir, seed=seed)
        except Exception as e:
            result.reason = f"ligand prep failed: {e}"
            return result

        out_pdbqt = workdir / "out.pdbqt"
        vina_bin = shutil.which("vina")
        if vina_bin is None:
            result.reason = "vina executable not found in PATH"
            return result

        cmd = [
            vina_bin,
            "--receptor", str(rec_pdbqt),
            "--ligand", str(lig_pdbqt),
            "--center_x", f"{center_A[0]:.3f}",
            "--center_y", f"{center_A[1]:.3f}",
            "--center_z", f"{center_A[2]:.3f}",
            "--size_x", f"{box_size_A[0]:.3f}",
            "--size_y", f"{box_size_A[1]:.3f}",
            "--size_z", f"{box_size_A[2]:.3f}",
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--seed", str(seed),
            "--out", str(out_pdbqt),
        ]
        if cpu > 0:
            cmd += ["--cpu", str(cpu)]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout_s)
        except subprocess.TimeoutExpired:
            result.reason = f"vina timeout after {timeout_s} s"
            result.wall_seconds = time.time() - t0
            return result
        except Exception as e:
            result.reason = f"vina invocation failed: {e}"
            result.wall_seconds = time.time() - t0
            return result

        if proc.returncode != 0:
            result.reason = (
                f"vina exited {proc.returncode}: "
                f"{(proc.stderr or proc.stdout)[:400]}")
            result.wall_seconds = time.time() - t0
            return result

        try:
            result.poses = _parse_vina_output(out_pdbqt)
        except Exception as e:
            result.reason = f"output parse failed: {e}"
            result.wall_seconds = time.time() - t0
            return result

    result.wall_seconds = time.time() - t0

    if not result.poses:
        result.reason = "no poses produced"
        return result

    result.best_kcalmol = result.poses[0].affinity_kcalmol
    result.ok = True

    # -------- cache put (only ok=True) ------------------------
    if cache is not None and cache_key is not None:
        try:
            cache.put(cache_key, cache_method, result.as_dict())
        except Exception as e:
            logger.debug("vina cache put failed: %s", e)
    return result


def _inflate_docking_result(value: dict) -> DockingResult:
    """Reconstruct a DockingResult (including DockingPose list) from
    a cached JSON dict."""
    pose_dicts = value.pop("poses", None) or []
    poses = []
    for pd in pose_dicts:
        pf = {k: v for k, v in pd.items()
              if k in DockingPose.__dataclass_fields__}
        poses.append(DockingPose(**pf))
    fields = {k: v for k, v in value.items()
              if k in DockingResult.__dataclass_fields__}
    # Center / box can come back as lists; normalise to tuples.
    if isinstance(fields.get("center_A"), list):
        fields["center_A"] = tuple(fields["center_A"])
    if isinstance(fields.get("box_size_A"), list):
        fields["box_size_A"] = tuple(fields["box_size_A"])
    result = DockingResult(**fields)
    result.poses = poses
    return result


def attach_reference_rmsd(
    result: DockingResult, reference_positions_A: list
) -> DockingResult:
    """Fill in `rmsd_vs_reference_A` for each pose using an RDKit
    no-alignment all-atom RMSD against the supplied reference.

    `reference_positions_A` must have the same element order and
    length as each pose's `positions_A`. If it doesn't, the RMSD is
    skipped for that pose (reason logged).
    """
    import numpy as np

    ref = np.array(reference_positions_A, dtype=float)
    for p in result.poses:
        pose_arr = np.array(p.positions_A, dtype=float)
        if pose_arr.shape != ref.shape:
            logger.debug(
                "ref shape %s != pose shape %s; skipping rmsd",
                ref.shape, pose_arr.shape)
            continue
        c_ref = ref.mean(axis=0)
        c_pose = pose_arr.mean(axis=0)
        d = (pose_arr - c_pose) - (ref - c_ref)
        p.rmsd_vs_reference_A = float(
            np.sqrt((d * d).sum(axis=1).mean()))
    if result.poses and result.poses[0].rmsd_vs_reference_A is not None:
        result.best_rmsd_vs_reference_A = \
            result.poses[0].rmsd_vs_reference_A
    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", required=True,
                    help="minimised receptor PDB")
    ap.add_argument("--ligand-smiles", required=True)
    ap.add_argument("--center", required=True,
                    help="comma-separated x,y,z in Å "
                         "(centre of search box)")
    ap.add_argument("--box", default="22,22,22",
                    help="comma-separated dx,dy,dz in Å (default 22,22,22)")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=0)
    args = ap.parse_args()

    center = tuple(float(x) for x in args.center.split(","))
    box = tuple(float(x) for x in args.box.split(","))

    r = dock_ligand(
        args.receptor, args.ligand_smiles,
        center_A=center, box_size_A=box,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes, seed=args.seed, cpu=args.cpu)
    print(r.summary())
    sys.exit(0 if r.ok else 1)
