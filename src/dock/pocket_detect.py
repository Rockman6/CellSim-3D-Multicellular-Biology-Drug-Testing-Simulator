"""Auto-detect druggable binding pockets on a protein PDB (non-AI).

Wrapper around **fpocket** (Le Guilloux et al 2009, BMC
Bioinformatics 10:168), the canonical non-ML geometric pocket
detector. fpocket builds a Voronoi diagram of the protein surface
and clusters alpha-spheres into candidate pockets, then ranks them
by a trained logistic "drug score" — which despite its name is a
fixed published model (not updated against new data), deterministic
given a PDB input.

From a biologist's perspective this solves the "what's the
--center argument for my own protein?" problem that would otherwise
require a cocrystal. You hand CellSim a PDB; it hands you back a
list of ranked pocket centroids with drug scores, volumes, and
recommended Vina search boxes.

Output: list of `PocketCandidate` records sorted by drug score
(descending = most druggable first).

Usage:
    from src.dock import detect_pockets
    pockets = detect_pockets("my_receptor.pdb")
    for p in pockets[:3]:
        print(p.summary())

    # drop the top pocket straight into dock_ligand:
    r = dock_ligand(
        "my_receptor.pdb", "CC(=O)OC1=CC=CC=C1C(=O)O",
        center_A=pockets[0].center_A,
        box_size_A=pockets[0].suggested_box_A,
        ...)
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
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class PocketCandidate:
    """One fpocket-detected binding-site candidate."""

    rank: int                      # 1 = most druggable
    drug_score: Optional[float] = None
    pocket_score: Optional[float] = None
    volume_A3: Optional[float] = None
    n_alpha_spheres: Optional[int] = None
    hydrophobicity_score: Optional[float] = None
    polarity_score: Optional[int] = None

    # Geometry
    center_A: Optional[tuple] = None           # (x, y, z)
    suggested_box_A: Optional[tuple] = None    # (dx, dy, dz)
    bounding_box_min_A: Optional[tuple] = None
    bounding_box_max_A: Optional[tuple] = None

    # Residues lining the pocket (for display / papers).
    residues: list = field(default_factory=list)   # [{chain, id, name}]
    n_atoms: Optional[int] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        box = (f"{self.suggested_box_A[0]:.1f}×"
               f"{self.suggested_box_A[1]:.1f}×"
               f"{self.suggested_box_A[2]:.1f}"
               if self.suggested_box_A else "?")
        return (f"pocket#{self.rank}  drug={self.drug_score:.3f}  "
                f"volume={self.volume_A3:.0f} Å³  "
                f"center=({self.center_A[0]:.1f},"
                f"{self.center_A[1]:.1f},"
                f"{self.center_A[2]:.1f})  "
                f"box={box}  residues={len(self.residues)}")


_HEADER_RE = re.compile(
    r"HEADER\s+\d+\s*-\s*(.+?)\s*:\s*(.+)$")


def _parse_pocket_atm(path: Path) -> PocketCandidate:
    """Parse fpocket's pocketN_atm.pdb into a PocketCandidate."""
    import numpy as np

    text = path.read_text()
    fields: dict = {}
    positions: list = []
    residues_seen: dict = {}

    for line in text.splitlines():
        m = _HEADER_RE.match(line)
        if m:
            key = m.group(1).strip().lower()
            val = m.group(2).strip()
            fields[key] = val
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            positions.append([x, y, z])
            chain = line[21]
            try:
                rid = int(line[22:26])
            except ValueError:
                continue
            rname = line[17:20].strip()
            key = (chain, rid)
            if key not in residues_seen:
                residues_seen[key] = rname

    pc = PocketCandidate(rank=0)

    # Numeric fields
    def _get_float(k: str) -> Optional[float]:
        v = fields.get(k)
        if v is None:
            return None
        try:
            return float(v.split()[0])
        except (ValueError, IndexError):
            return None

    pc.drug_score = _get_float("drug score")
    pc.pocket_score = _get_float("pocket score")
    pc.volume_A3 = _get_float("pocket volume (monte carlo)")
    n_alpha = _get_float("number of alpha spheres")
    pc.n_alpha_spheres = int(n_alpha) if n_alpha is not None else None
    pc.hydrophobicity_score = _get_float("hydrophobicity score")
    pol = _get_float("polarity score")
    pc.polarity_score = int(pol) if pol is not None else None

    if positions:
        arr = np.array(positions)
        ctr = arr.mean(axis=0)
        mn, mx = arr.min(axis=0), arr.max(axis=0)
        # Suggested Vina search box: pocket extent + 6 Å padding on
        # each side, clamped to a minimum of 18 Å (Vina struggles
        # with boxes smaller than ~15 Å for any reasonable ligand)
        # and maximum 30 Å (larger boxes explode the search space).
        extent = (mx - mn)
        box = [min(30.0, max(18.0, float(e + 12.0))) for e in extent]
        pc.center_A = tuple(float(c) for c in ctr)
        pc.suggested_box_A = tuple(box)
        pc.bounding_box_min_A = tuple(float(v) for v in mn)
        pc.bounding_box_max_A = tuple(float(v) for v in mx)
        pc.n_atoms = len(positions)

    pc.residues = [
        dict(chain=k[0], id=k[1], name=v)
        for k, v in sorted(residues_seen.items(),
                           key=lambda kv: (kv[0][0], kv[0][1]))
    ]
    return pc


def detect_pockets(
    pdb_path: str | Path,
    *,
    keep_fpocket_output: bool = False,
    top_k: int = 10,
) -> list[PocketCandidate]:
    """Run fpocket on a PDB and return the top-k candidates ranked
    by fpocket's drug score (descending).

    Returns an empty list if fpocket is unavailable or produces no
    pockets; never raises for recoverable failures.
    """
    pdb_path = Path(pdb_path)
    fp_bin = shutil.which("fpocket")
    if fp_bin is None:
        logger.warning("fpocket not on PATH (activate cellsim conda env)")
        return []
    if not pdb_path.exists():
        logger.warning("PDB not found: %s", pdb_path)
        return []

    # fpocket writes to a sibling directory `<pdb_stem>_out/` next to
    # the input PDB — copy into a temp dir so we don't pollute the
    # repo's benchmarks/ tree.
    with tempfile.TemporaryDirectory(prefix="cellsim-fpocket-") as tmp:
        tmp = Path(tmp)
        local_pdb = tmp / pdb_path.name
        local_pdb.write_bytes(pdb_path.read_bytes())
        cmd = [fp_bin, "-f", str(local_pdb)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=300)
        except subprocess.TimeoutExpired:
            logger.warning("fpocket timed out on %s", pdb_path)
            return []
        if r.returncode != 0:
            logger.warning("fpocket exit %d: %s",
                           r.returncode, (r.stderr or r.stdout)[:300])
            return []

        out_dir = tmp / f"{local_pdb.stem}_out" / "pockets"
        if not out_dir.exists():
            logger.warning("fpocket produced no pockets dir")
            return []

        candidates: list[PocketCandidate] = []
        for pdb_file in sorted(out_dir.glob("pocket*_atm.pdb")):
            try:
                pc = _parse_pocket_atm(pdb_file)
                candidates.append(pc)
            except Exception as e:
                logger.debug("failed to parse %s: %s", pdb_file, e)

        # Sort by drug score desc (most druggable first).
        candidates.sort(
            key=lambda p: (-(p.drug_score or -1.0), -(p.volume_A3 or 0.0)))
        for i, pc in enumerate(candidates[:top_k], 1):
            pc.rank = i
        if keep_fpocket_output:
            dst = Path(str(pdb_path) + ".fpocket_out")
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(tmp / f"{local_pdb.stem}_out", dst)
            logger.info("fpocket output preserved at %s", dst)
    return candidates[:top_k]


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--keep", action="store_true",
                    help="preserve fpocket's output dir next to "
                         "the input PDB for inspection")
    args = ap.parse_args()

    pockets = detect_pockets(
        args.pdb, top_k=args.top, keep_fpocket_output=args.keep)
    if not pockets:
        print("no pockets detected")
        sys.exit(1)
    print(f"[fpocket] {len(pockets)} top pockets on "
          f"{Path(args.pdb).name}:")
    for p in pockets:
        print(f"  {p.summary()}")
