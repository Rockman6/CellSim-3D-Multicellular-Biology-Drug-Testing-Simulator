"""CYP3A4 inhibition prediction via docking into the heme active site.

CYP3A4 metabolises ~50 % of marketed drugs; inhibitors cause
drug-drug interactions (DDI) via competitive blockade of substrate
oxidation. Every pharma pipeline screens for CYP3A4 inhibition
early (ICH E14, FDA DDI guidance).

CellSim's non-AI approach: dock the compound into a published
CYP3A4 crystal structure (1TQN, 2.05 Å) with the search box
centred on the heme iron. A compound that docks tight (< -7
kcal/mol) and places any atom within ~7 Å of Fe is a likely
competitive inhibitor.

Honest caveats (documented):
  - Vina's empirical scoring doesn't model metal coordination,
    so imidazole / pyridine coordinating inhibitors (type II
    binders: ketoconazole, clotrimazole, posaconazole) are
    under-ranked.
  - True hit requires IC50 measurement in microsomes or
    recombinant CYP3A4 — this is a triage filter, not a
    replacement.

Usage:
    from src.dock import cyp3a4_inhibition
    r = cyp3a4_inhibition("CC(=O)OC1=CC=CC=C1C(=O)O")
    print(r.summary())
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import dock_ligand  # noqa: E402

logger = logging.getLogger(__name__)

# Heme-Fe coordinate in 1TQN (chain A, residue 508). Taken from
# HETATM row; see benchmarks/dock/1tqn.pdb line starting
# 'HETATM 3810 FE   HEM A 508'.
_FE_COORD_1TQN = (-15.846, -23.032, -11.293)
_CYP3A4_PDB = REPO_ROOT / "benchmarks" / "dock" / "1tqn.pdb"
# Active-site box centred on Fe + offset upward along +y (the
# known substrate-access direction in CYP3A4). 28-Å cube covers
# the wide active cavity.
_CYP3A4_BOX = (28.0, 28.0, 28.0)


@dataclass
class CypInhibitionResult:
    """Envelope for a CYP3A4 inhibition screen."""

    smiles: str
    ok: bool
    reason: str = ""
    ligand_formula: Optional[str] = None

    dG_kcalmol: Optional[float] = None
    Kd_implied_nM: Optional[float] = None
    min_distance_to_Fe_A: Optional[float] = None

    inhibitor_risk: Optional[str] = None     # "low" | "medium" | "high"
    classification_reason: str = ""
    wall_seconds: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] CYP3A4 inhibition {self.smiles}  {self.reason}"
        return (f"[OK]   CYP3A4 inhibition  {self.ligand_formula or self.smiles}"
                f"  ΔG = {self.dG_kcalmol:+.2f} kcal/mol  "
                f"min(Fe-atom) = {self.min_distance_to_Fe_A:.2f} Å  "
                f"risk = {self.inhibitor_risk}  "
                f"({self.classification_reason})  "
                f"wall {self.wall_seconds:.1f} s")


def _classify(dG_kcalmol: float, fe_dist_A: float
               ) -> tuple[str, str]:
    """Classify DDI / inhibition risk."""
    if dG_kcalmol < -8.0 and fe_dist_A < 7.0:
        return "high", (f"ΔG {dG_kcalmol:+.2f} < -8 and Fe-atom "
                        f"{fe_dist_A:.1f} Å < 7 Å")
    if dG_kcalmol < -7.0 or fe_dist_A < 5.0:
        return "medium", (f"ΔG {dG_kcalmol:+.2f} and Fe-atom "
                          f"{fe_dist_A:.1f} Å above either "
                          "medium threshold")
    return "low", (f"ΔG {dG_kcalmol:+.2f} > -7 and Fe-atom "
                   f"{fe_dist_A:.1f} Å > 5 Å — below thresholds")


def cyp3a4_inhibition(
    smiles: str,
    *,
    exhaustiveness: int = 32,
    num_modes: int = 3,
    seed: int = 1,
    cpu: int = 2,
    cache: Optional["object"] = None,
) -> CypInhibitionResult:
    """Predict CYP3A4 inhibition / DDI risk for a compound.

    Docks into 1TQN (CYP3A4 crystal with heme, 2.05 Å) with a
    28-Å box on the heme Fe, then classifies risk from
    (ΔG, min-atom-to-Fe-distance).
    """
    import math
    import time

    result = CypInhibitionResult(smiles=smiles, ok=False)
    if not _CYP3A4_PDB.exists():
        result.reason = f"CYP3A4 PDB missing: {_CYP3A4_PDB}"
        return result

    t0 = time.time()
    r = dock_ligand(
        _CYP3A4_PDB, smiles,
        center_A=_FE_COORD_1TQN,
        box_size_A=_CYP3A4_BOX,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes, seed=seed, cpu=cpu,
        cache=cache)
    if not r.ok or not r.poses:
        result.reason = f"CYP3A4 dock failed: {r.reason}"
        result.wall_seconds = time.time() - t0
        return result

    top = r.poses[0]
    result.ligand_formula = r.ligand_formula
    result.dG_kcalmol = float(top.affinity_kcalmol)
    try:
        import math
        result.Kd_implied_nM = math.exp(
            top.affinity_kcalmol / (1.987204e-3 * 298.15)) * 1e9
    except Exception:
        pass

    # Minimum distance from any ligand heavy atom to the Fe.
    def _sq(a, b):
        return sum((ai - bi) ** 2 for ai, bi in zip(a, b))
    heavy = [(e, xyz) for e, xyz in zip(top.elements, top.positions_A)
             if e.upper() != "H"]
    min_d2 = min(_sq(xyz, _FE_COORD_1TQN) for _, xyz in heavy)
    result.min_distance_to_Fe_A = math.sqrt(min_d2)

    risk, reason = _classify(
        result.dG_kcalmol, result.min_distance_to_Fe_A)
    result.inhibitor_risk = risk
    result.classification_reason = reason

    result.ok = True
    result.wall_seconds = time.time() - t0
    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--num-modes", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--cache", default=None,
                    help="SQLite cache path for Vina reuse")
    args = ap.parse_args()

    cache = None
    if args.cache:
        from src.cache import Cache
        cache = Cache(args.cache)

    r = cyp3a4_inhibition(
        args.smiles,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes, seed=args.seed, cpu=args.cpu,
        cache=cache)
    print(r.summary())
    sys.exit(0 if r.ok else 1)
