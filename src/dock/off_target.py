"""Off-target / selectivity screen — one compound vs N receptors.

Drop-in wet-lab toxicity triage: given a lead compound and a set
of potential off-target receptors (kinases, GPCRs, CYPs, ion
channels), dock it into each and report a ranked per-receptor ΔG.
fpocket auto-detects the binding site on every receptor so the
biologist doesn't need to know the pocket coordinates up front.

Output:
    target                ΔG(kcal/mol)   K_d       fpocket_drug
    target_A (my target)        -9.5     1.2 nM    0.81
    kinase_B (off-target)       -7.1     8.4 µM    0.74
    gpcr_X                      -6.0    100 µM    0.52
    ...

The canonical selectivity story for a lead compound: its ΔG drop
from the intended target to the strongest off-target is the
in-silico proxy for selectivity (biologists want >= 3 kcal/mol =
~200× K_d selectivity).
"""

from __future__ import annotations

import logging
import math
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import detect_pockets, dock_ligand  # noqa: E402

logger = logging.getLogger(__name__)


@dataclass
class OffTargetEntry:
    """One receptor's result in the off-target panel."""

    name: str
    receptor_pdb: str
    ok: bool
    reason: str = ""
    dG_kcalmol: Optional[float] = None
    dG_kJmol: Optional[float] = None
    kd_implied_nM: Optional[float] = None
    pocket_drug_score: Optional[float] = None
    pocket_center_A: Optional[tuple] = None
    pocket_box_A: Optional[tuple] = None
    wall_seconds: Optional[float] = None
    pb_pocket_ok: Optional[bool] = None

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class OffTargetResult:
    """Envelope for a full off-target panel run."""

    ligand_smiles: str
    ligand_formula: Optional[str] = None
    entries: list = field(default_factory=list)  # list[OffTargetEntry]
    wall_seconds: Optional[float] = None

    def sorted_by_affinity(self) -> list:
        """Entries sorted by ΔG ascending (strongest binder first)."""
        ok_first = [e for e in self.entries if e.ok and e.dG_kcalmol is not None]
        ok_first.sort(key=lambda e: e.dG_kcalmol)
        failed = [e for e in self.entries if e not in ok_first]
        return ok_first + failed

    def selectivity_kcalmol(self) -> Optional[float]:
        """ΔG(best) − ΔG(second-best). Positive = the best target
        is preferentially bound vs the next best off-target."""
        ranked = self.sorted_by_affinity()
        ok = [e for e in ranked if e.ok and e.dG_kcalmol is not None]
        if len(ok) < 2:
            return None
        return ok[1].dG_kcalmol - ok[0].dG_kcalmol

    def summary(self) -> str:
        ranked = self.sorted_by_affinity()
        lines = [f"[OK]  off-target  {self.ligand_formula or self.ligand_smiles}  "
                 f"n_receptors={len(self.entries)}  "
                 f"({self.wall_seconds:.1f} s)"]
        for i, e in enumerate(ranked, 1):
            if e.ok:
                kd = e.kd_implied_nM or float("nan")
                kd_s = (f"{kd:.1f} nM" if kd < 1e3 else
                        f"{kd/1e3:.1f} µM" if kd < 1e6 else
                        f"{kd/1e6:.1f} mM")
                pocket = (f"pb:ok" if e.pb_pocket_ok else
                          f"pb:? " if e.pb_pocket_ok is None else
                          f"pb:no")
                drug = (f"drug={e.pocket_drug_score:.2f}"
                        if e.pocket_drug_score is not None else
                        "drug=?  ")
                lines.append(
                    f"  #{i:>2d}  {e.name:<30s}  "
                    f"ΔG = {e.dG_kcalmol:+.2f}  K_d ≈ {kd_s:<10s}  "
                    f"{drug}  {pocket}")
            else:
                lines.append(
                    f"  -   {e.name:<30s}  FAIL  {e.reason[:60]}")
        sel = self.selectivity_kcalmol()
        if sel is not None:
            lines.append(
                f"  selectivity (best vs 2nd-best): "
                f"ΔΔG = {sel:+.2f} kcal/mol  "
                f"({'good' if sel > 3.0 else 'weak'})")
        return "\n".join(lines)

    def as_dict(self) -> dict:
        return asdict(self)


def _kd_from_dG_nM(dG_kcalmol: float, T_K: float = 298.15) -> float:
    R_kcal = 1.987204e-3
    return math.exp(dG_kcalmol / (R_kcal * T_K)) * 1e9


def off_target_screen(
    ligand_smiles: str,
    receptors: list,  # list of (name, pdb_path) tuples OR list of PDB paths
    *,
    exhaustiveness: int = 32,
    num_modes: int = 3,
    seed: int = 1,
    cpu_per_job: int = 2,
    timeout_s: int = 600,
    cache: Optional["object"] = None,
    auto_pocket_rank: int = 1,
) -> OffTargetResult:
    """Dock one ligand into N receptors; report ranked ΔG per receptor.

    `receptors` is either `[(name, pdb_path), …]` or `[pdb_path,
    …]` (names inferred from path stems).

    Each receptor's binding-site is auto-detected via fpocket
    (top-ranked pocket by default). If fpocket can't find a pocket
    that entry fails with a readable reason.
    """
    # Normalise input.
    norm: list[tuple[str, str]] = []
    for item in receptors:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            name, pdb = item
        else:
            pdb = str(item)
            name = Path(pdb).stem
        norm.append((name, pdb))

    result = OffTargetResult(
        ligand_smiles=ligand_smiles, entries=[])

    # Resolve formula for the result envelope.
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(ligand_smiles)
        if mol is not None:
            result.ligand_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        pass

    t0_all = time.time()
    for name, pdb in norm:
        t0 = time.time()
        entry = OffTargetEntry(name=name, receptor_pdb=pdb, ok=False)

        # Auto pocket.
        pockets = detect_pockets(pdb, top_k=auto_pocket_rank + 1)
        if not pockets or len(pockets) < auto_pocket_rank:
            entry.reason = (f"fpocket found {len(pockets)} pockets; "
                            f"cannot satisfy rank={auto_pocket_rank}")
            entry.wall_seconds = time.time() - t0
            result.entries.append(entry)
            print(f"  [off-target] {name:<30s}  FAIL: {entry.reason}",
                  flush=True)
            continue
        pkt = pockets[auto_pocket_rank - 1]
        entry.pocket_drug_score = pkt.drug_score
        entry.pocket_center_A = pkt.center_A
        entry.pocket_box_A = pkt.suggested_box_A

        # Dock.
        r = dock_ligand(
            pdb, ligand_smiles,
            center_A=pkt.center_A, box_size_A=pkt.suggested_box_A,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes, seed=seed,
            cpu=cpu_per_job, timeout_s=timeout_s, cache=cache)
        if not r.ok or not r.poses:
            entry.reason = r.reason or "no poses"
            entry.wall_seconds = time.time() - t0
            result.entries.append(entry)
            print(f"  [off-target] {name:<30s}  FAIL: {entry.reason}",
                  flush=True)
            continue
        top = r.poses[0]
        entry.ok = True
        entry.dG_kcalmol = float(top.affinity_kcalmol)
        entry.dG_kJmol = float(top.affinity_kJmol)
        entry.kd_implied_nM = float(_kd_from_dG_nM(top.affinity_kcalmol))
        entry.pb_pocket_ok = top.posebusters_pocket_ok
        entry.wall_seconds = time.time() - t0
        result.entries.append(entry)
        print(f"  [off-target] {name:<30s}  "
              f"ΔG = {entry.dG_kcalmol:+.2f}  "
              f"drug_score = {pkt.drug_score:.2f}  "
              f"({entry.wall_seconds:.1f} s)",
              flush=True)

    result.wall_seconds = time.time() - t0_all
    return result


def write_csv(result: OffTargetResult, path: str | Path) -> None:
    import csv
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["rank", "name", "receptor_pdb", "ok", "reason",
            "dG_kcalmol", "dG_kJmol", "Kd_nM",
            "pocket_drug_score", "pb_pocket_ok",
            "pocket_center_x", "pocket_center_y", "pocket_center_z",
            "wall_s"]
    with path.open("w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(cols)
        for rank, e in enumerate(result.sorted_by_affinity(), 1):
            row = [
                rank if e.ok else "",
                e.name, e.receptor_pdb, e.ok, e.reason,
                (f"{e.dG_kcalmol:.3f}" if e.dG_kcalmol is not None
                 else ""),
                (f"{e.dG_kJmol:.3f}" if e.dG_kJmol is not None
                 else ""),
                (f"{e.kd_implied_nM:.3g}"
                 if e.kd_implied_nM is not None else ""),
                (f"{e.pocket_drug_score:.3f}"
                 if e.pocket_drug_score is not None else ""),
                e.pb_pocket_ok if e.pb_pocket_ok is not None else "",
                (f"{e.pocket_center_A[0]:.2f}"
                 if e.pocket_center_A else ""),
                (f"{e.pocket_center_A[1]:.2f}"
                 if e.pocket_center_A else ""),
                (f"{e.pocket_center_A[2]:.2f}"
                 if e.pocket_center_A else ""),
                (f"{e.wall_seconds:.1f}"
                 if e.wall_seconds is not None else ""),
            ]
            w.writerow(row)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ligand-smiles", required=True)
    ap.add_argument("--receptors", required=True,
                    help="comma-separated PDB paths OR a file "
                         "containing 'name<TAB>pdb' per line")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu-per-job", type=int, default=2)
    ap.add_argument("--cache", default=None,
                    help="SQLite cache path")
    ap.add_argument("--out-csv", default=None)
    ap.add_argument("--save-plot", default=None,
                    help="save bar-chart PNG of ΔG per receptor")
    ap.add_argument("--intended-target", default=None,
                    help="highlight this receptor name in the plot")
    args = ap.parse_args()

    receptors: list = []
    p = Path(args.receptors)
    if p.exists() and p.is_file() and p.suffix in (".tsv", ".txt"):
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split(None, 1)
            if len(parts) == 2:
                receptors.append((parts[0], parts[1]))
            else:
                receptors.append(parts[0])
    else:
        for item in args.receptors.split(","):
            item = item.strip()
            if not item:
                continue
            if "=" in item:
                # name=path syntax
                name, path = item.split("=", 1)
                receptors.append((name.strip(), path.strip()))
            else:
                receptors.append(item)

    cache = None
    if args.cache:
        from src.cache import Cache
        cache = Cache(args.cache)

    r = off_target_screen(
        args.ligand_smiles, receptors,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes, seed=args.seed,
        cpu_per_job=args.cpu_per_job, cache=cache)
    print()
    print(r.summary())
    if args.out_csv:
        write_csv(r, args.out_csv)
        print(f"\n[off-target] wrote {args.out_csv}")
    if args.save_plot:
        from src.dock.off_target_viewer import render_off_target_result
        render_off_target_result(
            r, intended_target=args.intended_target,
            save=args.save_plot, show=False)
    sys.exit(0)
