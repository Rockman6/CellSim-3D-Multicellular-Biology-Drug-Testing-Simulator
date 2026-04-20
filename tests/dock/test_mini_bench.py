"""Layer 1.3 mini-docking benchmark — aggregate pose-recovery %.

Runs re-docking on every entry in benchmarks/dock/mini_bench.yaml
and reports the aggregate pass-rate using the canonical gate:

    top-1 RMSD <= 2.5 Å  AND  best-of-top-3 < 2.0 Å

This is the industry-standard "how good is your docker?" number.
Vina papers report ~ 75-85 % on the Astex Diverse Set (85 cocrystals)
at exhaustiveness=32. Our mini set is curated for difficulty
diversity (small/medium/large ligands, four unrelated protein
families), so per-entry success matters more than a precise %.

Gate: >= 3/4 entries pass. Lower than Astex's 75 % because our set
deliberately includes a hard case (indinavir with ~ 35 heavy atoms
and many rotatable bonds; Vina's published top-3 success rate on
systems this size is ~ 60 %).

Run:
    conda activate cellsim
    python tests/dock/test_mini_bench.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import (  # noqa: E402
    attach_crystal_rmsd, dock_ligand,
)

try:
    import yaml  # PyYAML comes with environment.yml
except ImportError:  # degrade gracefully
    yaml = None

BENCH = REPO_ROOT / "benchmarks" / "dock" / "mini_bench.yaml"

TOP1_GATE_A = 2.5
TOP3_GATE_A = 2.0
PASS_FRACTION_GATE = 0.66   # >= 2/3 on the 3-cocrystal set


def _load_entries() -> list[dict]:
    if yaml is None:
        # Minimal fallback parser — pyyaml is in environment.yml so
        # only trips in an under-provisioned env.
        raise RuntimeError(
            "PyYAML not installed — activate cellsim conda env")
    data = yaml.safe_load(BENCH.read_text())
    return data.get("entries", [])


def _run_one(entry: dict, *, exhaustiveness: int = 64,
             num_modes: int = 9, seed: int = 1,
             cpu: int = 2) -> dict:
    pdb = REPO_ROOT / "benchmarks" / "dock" / entry["pdb"]
    smi = entry["ligand_smiles"]
    resname = str(entry["ligand_resname"])
    center = tuple(entry["crystal_center_A"])
    box = tuple(entry["search_box_A"])

    t0 = time.time()
    r = dock_ligand(pdb, smi, center_A=center, box_size_A=box,
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes, seed=seed, cpu=cpu)
    if r.ok:
        try:
            r = attach_crystal_rmsd(
                r, crystal_pdb=pdb, ligand_resname=resname)
        except Exception as e:
            print(f"  WARN: crystal RMSD failed: {e}")
    wall = time.time() - t0

    rec = dict(name=entry["name"], pdb=entry["pdb"], ok=r.ok,
               reason=r.reason if not r.ok else "",
               wall_s=round(wall, 1))
    if not r.ok:
        return rec
    top = r.poses[0]
    top3 = [p for p in r.poses[:3]
            if p.rmsd_vs_reference_A is not None]
    top3_best = (min((p.rmsd_vs_reference_A for p in top3))
                 if top3 else None)
    rec.update(
        dG_kcalmol=round(top.affinity_kcalmol, 2),
        top1_rmsd_A=(round(top.rmsd_vs_reference_A, 2)
                     if top.rmsd_vs_reference_A is not None else None),
        top3_best_rmsd_A=(round(top3_best, 2)
                          if top3_best is not None else None),
    )
    top1 = rec["top1_rmsd_A"]
    top3b = rec["top3_best_rmsd_A"]
    rec["pass_gate"] = (
        top1 is not None and top1 <= TOP1_GATE_A
        and top3b is not None and top3b < TOP3_GATE_A)
    return rec


def test_mini_bench():
    entries = _load_entries()
    assert entries, "no entries in mini_bench.yaml"

    t0 = time.time()
    rows: list[dict] = []
    for entry in entries:
        print(f"[bench] {entry['name']:<32s} ({entry['pdb']})",
              flush=True)
        rec = _run_one(entry)
        rows.append(rec)
        if rec["ok"]:
            tag = "✓ PASS" if rec.get("pass_gate") else "✗ FAIL"
            print(f"  {tag}  ΔG={rec['dG_kcalmol']:+.2f} kcal/mol  "
                  f"top1={rec.get('top1_rmsd_A')} Å  "
                  f"top3b={rec.get('top3_best_rmsd_A')} Å  "
                  f"({rec['wall_s']:.1f} s)")
        else:
            print(f"  ERROR: {rec['reason'][:100]}")

    passed = sum(1 for r in rows if r.get("pass_gate"))
    n = len(rows)
    frac = passed / n
    wall = time.time() - t0

    print()
    print(f"[bench] {passed}/{n} = {frac * 100:.0f} %  (gate "
          f"{PASS_FRACTION_GATE * 100:.0f} %)   wall={wall:.1f} s")
    for r in rows:
        tag = ("✓" if r.get("pass_gate")
               else "✗" if r["ok"]
               else "!")
        top1 = f"{r.get('top1_rmsd_A', '-')} Å" if r["ok"] else "-"
        top3b = f"{r.get('top3_best_rmsd_A', '-')} Å" if r["ok"] else "-"
        print(f"  {tag} {r['name']:<32s} top1 {top1:<8s} "
              f"top3_best {top3b}")

    assert frac >= PASS_FRACTION_GATE, (
        f"aggregate pose recovery {frac*100:.0f}% < "
        f"{PASS_FRACTION_GATE*100:.0f}% gate "
        f"({passed}/{n})")


if __name__ == "__main__":
    try:
        test_mini_bench()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
