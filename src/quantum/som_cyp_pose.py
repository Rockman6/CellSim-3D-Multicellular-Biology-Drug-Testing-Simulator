"""CYP3A4 SoM prediction with heme-accessibility filtering.

The honest limitation of BDE-only SoM prediction (documented in
BENCHMARKS.md): BDE ranking picks the most reactive C-H, but CYP
metabolism also requires the H to be within reach of the activating
Fe=O ferryl species (~6 Å from Fe). BDE ≠ kinetic preference.

This module closes that gap by:
  1. Docking the compound into 1TQN (CYP3A4 crystal with heme).
  2. Running the xTB BDE SoM predictor on the compound (independent
     of the docking).
  3. For the docked top pose, mapping SoM-candidate heavy atoms
     to pose atoms via SMILES heavy-atom order (Meeko preserves
     this through Vina).
  4. Computing distance from each SoM parent atom to the heme Fe.
  5. Filtering candidates: only atoms within `max_fe_distance_A`
     are 'accessible' — everything else is struck from the ranking.
  6. Top-1 = lowest BDE among the accessible set.

The resulting primary SoM should match the clinically reported
metabolism site on compounds that have both a low-BDE C-H and a
low-Fe-distance C-H — which covers most CYP3A4 substrates.

Non-AI. All physics: xTB BDE + Vina pose + Euclidean distance +
SMARTS-validation of the ranked outcome.
"""

from __future__ import annotations

import logging
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.cyp_inhibition import _FE_COORD_1TQN, _CYP3A4_PDB, _CYP3A4_BOX  # noqa: E402
from src.dock import dock_ligand  # noqa: E402
from src.quantum.metabolism import predict_cyp_som_bde, SoMResult  # noqa: E402

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
    if e in ("H", "HD", "HS"):
        return "H"
    return e.title()


def _parent_pose_idx(cand, som, pose) -> Optional[int]:
    """Map a SoMCandidate.parent_atom_idx to the matching heavy-atom
    index in the docked pose. Returns None on mismatch."""
    n_heavy_before = sum(
        1 for i, e in enumerate(som.elements)
        if i < cand.parent_atom_idx and e.upper() != "H")
    pose_heavy = [(i, _normalise_elem(e))
                  for i, e in enumerate(pose.elements)
                  if _normalise_elem(e) != "H"]
    if n_heavy_before >= len(pose_heavy):
        return None
    pose_full_idx, pose_elem = pose_heavy[n_heavy_before]
    if pose_elem.title() != cand.parent_element.title():
        return None
    return pose_full_idx


@dataclass
class CypPoseSoMResult:
    smiles: str
    ok: bool
    reason: str = ""

    n_candidates: int = 0
    fe_coord: tuple = _FE_COORD_1TQN
    max_fe_distance_A: float = 6.0

    # Per-candidate bundled info; ordered by accessible-BDE rank.
    accessible_rank: list = field(default_factory=list)
    # Each entry: dict(rank_post_filter, parent_atom_idx,
    #                  parent_element, bde_kcalmol, fe_distance_A,
    #                  accessible: bool)

    dock_dG_kcalmol: Optional[float] = None
    strain_band: Optional[str] = None
    strain_kcalmol: Optional[float] = None
    strain_ratio: Optional[float] = None
    wall_seconds: Optional[float] = None

    def top1_parent_idx(self) -> Optional[int]:
        for c in self.accessible_rank:
            if c.get("accessible"):
                return c["parent_atom_idx"]
        return None

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] CYP pose-SoM {self.smiles}  {self.reason}"
        strain_tag = (f"  strain={self.strain_band}"
                      if self.strain_band else "")
        lines = [
            f"[OK]   CYP3A4 pose-SoM  {self.smiles}  "
            f"n_candidates={self.n_candidates}  "
            f"dock_ΔG={self.dock_dG_kcalmol:+.2f}{strain_tag}  "
            f"max_Fe={self.max_fe_distance_A:.1f} Å  "
            f"wall {self.wall_seconds:.1f} s",
        ]
        for c in self.accessible_rank[:5]:
            tag = "✓" if c["accessible"] else "—"
            lines.append(
                f"  {tag} {c['parent_element']}(idx={c['parent_atom_idx']}): "
                f"BDE={c['bde_kcalmol']:.1f}  "
                f"Fe-dist={c['fe_distance_A']:.2f} Å")
        return "\n".join(lines)

    def as_dict(self) -> dict:
        return asdict(self)


def predict_cyp_som_with_heme_access(
    smiles: str,
    *,
    max_fe_distance_A: float = 6.0,
    exhaustiveness: int = 32,
    num_modes: int = 9,
    seed: int = 1,
    cpu: int = 2,
    cache: Optional["object"] = None,
    ensemble_pose_search: bool = True,
) -> CypPoseSoMResult:
    """Run xTB SoM + CYP3A4 docking; filter SoM by Fe accessibility.

    If `ensemble_pose_search=True` (default), iterate through ALL
    Vina poses and pick the one that places the lowest-BDE
    candidate closest to Fe (within `max_fe_distance_A`). This
    matches the biological fact that CYP3A4's active site is
    flexible and the catalytically-productive pose is not always
    Vina's best-scored pose — energy and catalytic competence
    are distinct properties.
    """
    import math

    result = CypPoseSoMResult(
        smiles=smiles, ok=False,
        max_fe_distance_A=max_fe_distance_A)
    t0 = time.time()

    # Step 1: xTB BDE ranking (independent of receptor).
    som = predict_cyp_som_bde(smiles)
    if not som.ok or not som.candidates:
        result.reason = f"SoM predictor failed: {som.reason}"
        return result
    result.n_candidates = len(som.candidates)

    # Step 2: dock into CYP3A4.
    dock_r = dock_ligand(
        _CYP3A4_PDB, smiles,
        center_A=_FE_COORD_1TQN, box_size_A=_CYP3A4_BOX,
        exhaustiveness=exhaustiveness, num_modes=num_modes,
        seed=seed, cpu=cpu, cache=cache)
    if not dock_r.ok or not dock_r.poses:
        result.reason = f"CYP3A4 dock failed: {dock_r.reason}"
        result.wall_seconds = time.time() - t0
        return result
    result.dock_dG_kcalmol = float(dock_r.poses[0].affinity_kcalmol)

    # Step 2a: ensemble search — pick the Vina pose that places the
    # lowest-BDE candidate within reach of Fe. Fall back to pose #1
    # if no pose reaches the heme for any candidate.
    top_pose = dock_r.poses[0]
    if ensemble_pose_search:
        best_score = None
        best_pose = top_pose
        for cand_pose in dock_r.poses:
            # For this pose, compute (min BDE among accessible
            # candidates, min Fe distance across all candidates).
            # Lower best_bde + any accessibility = better.
            min_accessible_bde = None
            for cand in som.candidates:
                pose_idx = _parent_pose_idx(cand, som, cand_pose)
                if pose_idx is None:
                    continue
                xyz = cand_pose.positions_A[pose_idx]
                d = math.sqrt(sum((xyz[k] - _FE_COORD_1TQN[k]) ** 2
                                    for k in range(3)))
                if d <= max_fe_distance_A:
                    if (min_accessible_bde is None
                            or cand.bde_kcalmol < min_accessible_bde):
                        min_accessible_bde = cand.bde_kcalmol
            if min_accessible_bde is not None:
                if best_score is None or min_accessible_bde < best_score:
                    best_score = min_accessible_bde
                    best_pose = cand_pose
        top_pose = best_pose

    # Step 3: map SoM candidate parent atoms → pose atoms.
    # The BDE predictor indexes atoms in its internal RDKit mol
    # (which preserves SMILES order + adds Hs). The docking pose
    # has heavy atoms in SMILES order followed by Hs. Both paths
    # use Meeko/RDKit which preserves the heavy-atom ordering, so
    # mapping is: som.positions_A[parent_atom_idx] is the embedded
    # geometry BUT the DOCKED geometry is on top_pose.positions_A.
    #
    # To find the corresponding heavy-atom index in the pose, we
    # need: (a) how many heavy atoms precede parent_atom_idx in
    # som.elements, and (b) find that heavy-atom position in the
    # pose (which lists heavy atoms first via Meeko's convention).
    #
    # Walk the SoM `elements` array; every time we hit the parent
    # idx, count how many heavy atoms came before and use that
    # as the pose index.
    entries = []
    for cand in som.candidates:
        pose_full_idx = _parent_pose_idx(cand, som, top_pose)
        if pose_full_idx is None:
            continue
        xyz = top_pose.positions_A[pose_full_idx]
        dx, dy, dz = (xyz[i] - _FE_COORD_1TQN[i] for i in range(3))
        fe_dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        entries.append(dict(
            parent_atom_idx=cand.parent_atom_idx,
            parent_element=cand.parent_element,
            bde_kcalmol=cand.bde_kcalmol,
            fe_distance_A=fe_dist,
            accessible=(fe_dist <= max_fe_distance_A),
        ))

    # Rank entries by BDE ascending WITHIN the accessible set first;
    # accessible candidates float to the top. Inaccessible entries
    # keep their BDE ordering but below any accessible ones.
    entries.sort(key=lambda c: (
        not c["accessible"], c["bde_kcalmol"]))
    for rank, e in enumerate(entries, 1):
        e["rank_post_filter"] = rank
    result.accessible_rank = entries

    if not any(c["accessible"] for c in entries):
        result.reason = (
            f"no SoM candidate within {max_fe_distance_A} Å of Fe; "
            "likely a docking pose that doesn't place any C-H in "
            "reach — try --max-fe > 6 or refine pose")
        result.wall_seconds = time.time() - t0
        return result

    # Pose-trust: UFF-ensemble strain on the top catalytically-
    # productive pose. If strain=reject, the SoM prediction is
    # resting on a non-physical docking pose and should not be
    # trusted clinically — flagged but not suppressed, so the
    # user can see the full ranking and make the call.
    try:
        from src.dock.strain import ligand_strain
        s = ligand_strain(
            top_pose.elements, top_pose.positions_A,
            smiles, ensemble_n=20)
        if s.ok:
            result.strain_band = s.band
            result.strain_kcalmol = round(s.strain_kcalmol, 2)
            result.strain_ratio = round(s.energy_ratio, 3)
    except Exception as e:
        logger.debug("strain on CYP pose-SoM failed: %s", e)

    result.ok = True
    result.wall_seconds = time.time() - t0
    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles")
    ap.add_argument("--max-fe", type=float, default=6.0,
                    help="max Fe-atom distance (Å) to count as "
                         "accessible (default 6.0)")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--cache", default=None)
    args = ap.parse_args()

    cache = None
    if args.cache:
        from src.cache import Cache
        cache = Cache(args.cache)

    r = predict_cyp_som_with_heme_access(
        args.smiles,
        max_fe_distance_A=args.max_fe,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes, seed=args.seed, cpu=args.cpu,
        cache=cache)
    print(r.summary())
    sys.exit(0 if r.ok else 1)
