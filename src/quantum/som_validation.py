"""Layer 1.4 CYP3A4 SoM literature-validation harness.

Given a YAML bundle of (SMILES, known_som_smarts) triples, runs
predict_cyp_som_bde (and optionally DFT rescoring) on each, then
checks whether the top-1 predicted parent atom index falls inside
the SMARTS-matched set for the literature site. Reports per-
compound correct / incorrect and aggregate accuracy.

Usage:
    from src.quantum.som_validation import run_som_validation
    r = run_som_validation(
        "benchmarks/quantum/cyp3a4_som_validation.yaml",
        dft_verify=3)
    print(r.summary())

CLI:
    cellsim som-validate benchmarks/quantum/cyp3a4_som_validation.yaml
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

from src.quantum.metabolism import (  # noqa: E402
    predict_cyp_som_bde, verify_som_dft,
)

logger = logging.getLogger(__name__)


@dataclass
class SoMValidationPoint:
    name: str
    smiles: str
    known_som_smarts: str
    known_site_text: str
    ok: bool = False
    reason: str = ""
    pred_top1_idx: Optional[int] = None
    pred_top1_element: Optional[str] = None
    pred_top1_bde: Optional[float] = None
    correct: Optional[bool] = None
    method_used: Optional[str] = None


@dataclass
class SoMValidationResult:
    yaml_path: str
    entries: list = field(default_factory=list)
    dft_verified: int = 0
    wall_seconds: Optional[float] = None

    def correct_count(self) -> int:
        return sum(1 for p in self.entries if p.correct is True)

    def incorrect_count(self) -> int:
        return sum(1 for p in self.entries if p.correct is False)

    def ok_count(self) -> int:
        return sum(1 for p in self.entries if p.ok)

    def summary(self) -> str:
        n = len(self.entries)
        ok = self.ok_count()
        correct = self.correct_count()
        lines = [
            f"[SoM validation]  {self.yaml_path}",
            f"  n = {n}  predictor ran ok: {ok}/{n}  "
            f"top-1 matches literature: {correct}/{ok}  "
            f"(DFT rescored: {self.dft_verified})",
            f"  wall = {(self.wall_seconds or 0.0):.1f} s",
        ]
        for p in self.entries:
            if not p.ok:
                tag = "!"
                note = p.reason
            elif p.correct:
                tag = "✓"
                note = (f"predicted top-1 at idx={p.pred_top1_idx} "
                        f"({p.pred_top1_element}) matches "
                        f"'{p.known_site_text}'")
            else:
                tag = "✗"
                note = (f"predicted top-1 at idx={p.pred_top1_idx} "
                        f"({p.pred_top1_element}) NOT in "
                        f"'{p.known_site_text}'")
            lines.append(f"  {tag} {p.name:<14s}  {note}")
        return "\n".join(lines)

    def as_dict(self) -> dict:
        return asdict(self)


def run_som_validation(
    yaml_path: str | Path,
    *,
    dft_verify: int = 0,
    heme_accessibility: bool = True,
    max_fe_distance_A: float = 10.0,
    cache_path: Optional[str] = None,
) -> SoMValidationResult:
    """Run the SoM literature-validation harness.

    `dft_verify`: if > 0, apply verify_som_dft(top_n=dft_verify) on
    each prediction before checking against the SMARTS. Much slower
    (minutes per compound) but more accurate.

    `heme_accessibility`: dock each compound into CYP3A4 (1TQN) and
    re-rank SoM candidates by (accessible first, then BDE ascending).
    Matches the biological mechanism: only C-H bonds within
    `max_fe_distance_A` of the heme iron can be oxidised by the ferryl
    species.

    DEFAULT ON (changed 2026-07). Bare BDE ranking asks only "which
    C-H is weakest?", never "can the enzyme actually reach it?", and
    that is not how CYP3A4 works — a buried weak C-H is not a site of
    metabolism. Measured on the bundled literature set, bare BDE scores
    1/3 and the heme re-rank 2/3, so the accurate path is now the
    default and the fast path is opt-out (`--no-heme-access`). Cost:
    one docking run per compound.
    """
    import yaml as pyyaml
    from rdkit import Chem
    from rdkit.Chem import AllChem

    data = pyyaml.safe_load(Path(yaml_path).read_text())
    entries_cfg = data.get("entries", []) or []

    result = SoMValidationResult(
        yaml_path=str(yaml_path),
        dft_verified=dft_verify,
    )
    t0 = time.time()

    # Open optional cache (shared across compounds).
    _cache = None
    if cache_path:
        try:
            from src.cache import Cache
            _cache = Cache(cache_path)
        except Exception as e:
            logger.debug("cache open failed: %s", e)

    for cfg in entries_cfg:
        p = SoMValidationPoint(
            name=cfg["name"],
            smiles=cfg["smiles"],
            known_som_smarts=cfg["known_som_smarts"],
            known_site_text=cfg.get("known_site_text", ""),
        )
        print(f"[som-validate] {p.name}", flush=True)

        # Run predictor.
        try:
            som = predict_cyp_som_bde(p.smiles)
        except Exception as e:
            p.reason = f"predictor crashed: {e}"
            result.entries.append(p)
            continue
        if not som.ok or not som.candidates:
            p.reason = f"predictor failed: {som.reason}"
            result.entries.append(p)
            continue

        # Optionally DFT-rescore.
        if dft_verify > 0:
            try:
                som = verify_som_dft(som, top_n=dft_verify)
            except Exception as e:
                print(f"  DFT rescoring skipped: {e}")

        # Heme-accessibility re-rank (CYP3A4 pose docking).
        pose_top1_idx = None
        pose_top1_elem = None
        pose_top1_bde = None
        if heme_accessibility:
            try:
                from src.quantum.som_cyp_pose import \
                    predict_cyp_som_with_heme_access
                cyp = predict_cyp_som_with_heme_access(
                    p.smiles,
                    max_fe_distance_A=max_fe_distance_A,
                    exhaustiveness=16, num_modes=3,
                    cache=_cache)
                if cyp.ok and cyp.accessible_rank:
                    top_acc = cyp.accessible_rank[0]
                    pose_top1_idx = top_acc["parent_atom_idx"]
                    pose_top1_elem = top_acc["parent_element"]
                    pose_top1_bde = top_acc["bde_kcalmol"]
                    p.method_used = (som.method +
                                     f" + heme_access(<{max_fe_distance_A} Å)")
            except Exception as e:
                print(f"  heme-accessibility rescoring failed: {e}")

        p.ok = True
        if p.method_used is None:
            p.method_used = som.method

        # Top-1 candidate: use heme-filtered rank if available.
        if pose_top1_idx is not None:
            p.pred_top1_idx = pose_top1_idx
            p.pred_top1_element = pose_top1_elem
            p.pred_top1_bde = pose_top1_bde
        else:
            top = som.candidates[0]
            p.pred_top1_idx = top.parent_atom_idx
            p.pred_top1_element = top.parent_element
            p.pred_top1_bde = top.bde_kcalmol

        # Check correctness: rebuild the RDKit mol from the SAME
        # pathway the predictor used (so atom indices match), then
        # run SMARTS match on heavy-atom skeleton.
        mol = Chem.AddHs(Chem.MolFromSmiles(p.smiles))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        patt = Chem.MolFromSmarts(p.known_som_smarts)
        if patt is None:
            p.reason = f"bad SMARTS: {p.known_som_smarts}"
            result.entries.append(p)
            continue

        # Collect all atom indices in ANY match of the SMARTS pattern.
        matches = mol.GetSubstructMatches(patt)
        valid_idxs: set[int] = set()
        for m in matches:
            valid_idxs.update(m)
        p.correct = (p.pred_top1_idx in valid_idxs)
        result.entries.append(p)

        tag = "✓" if p.correct else "✗"
        print(f"  {tag} top-1 idx={p.pred_top1_idx} "
              f"({p.pred_top1_element}); SMARTS matches "
              f"{len(matches)} positions "
              f"→ {'in match' if p.correct else 'NOT in match'}")

    result.wall_seconds = time.time() - t0
    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("yaml", help="validation-set YAML")
    ap.add_argument("--dft-verify", type=int, default=0, metavar="N",
                    help="after xTB top-N ranking, DFT-rescore those "
                         "N for every compound (slow but paper-grade)")
    ap.add_argument("--no-heme-access", action="store_true",
                    help="DISABLE the CYP3A4 (1TQN) heme-accessibility "
                         "re-rank. Faster, but markedly less accurate: "
                         "bare BDE ranking ignores whether a C-H can "
                         "physically reach the ferryl oxygen, and on "
                         "the bundled literature set it scores 1/3 vs "
                         "2/3 with the re-rank. Only use for a quick "
                         "look, never for a reported number.")
    ap.add_argument("--max-fe", type=float, default=10.0,
                    help="Fe-distance cutoff when --heme-access "
                         "(default 10.0 Å)")
    ap.add_argument("--cache", default=None,
                    help="SQLite cache for repeated runs")
    args = ap.parse_args()

    r = run_som_validation(
        args.yaml,
        dft_verify=args.dft_verify,
        heme_accessibility=not args.no_heme_access,
        max_fe_distance_A=args.max_fe,
        cache_path=args.cache)
    print()
    print(r.summary())
    print()
    ok = r.ok_count()
    correct = r.correct_count()
    pct = (correct / ok * 100) if ok else 0.0
    print(f"  AGGREGATE: {correct}/{ok} predictions match "
          f"literature ({pct:.0f} %)")
    sys.exit(0 if ok == len(r.entries) and correct == ok else 1)
