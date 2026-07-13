#!/usr/bin/env python
"""Layer 1.7 — blind pose-recovery + PoseBusters validation harness.

Campaign-1 exit criteria #1 (PDBBind blind pose recovery >= 75% within
2 Å) and #3 (PoseBusters physical validity >= 95%) require a real,
external, blind benchmark — not our own 3-cocrystal toy set. This
harness closes that gap.

For every entry `{pdb, resname?}` in the manifest it:
  1. fetches the crystal structure from RCSB (cached),
  2. picks the drug ligand (the given resname, else the largest
     non-buffer HETATM group) and fetches its AUTHORITATIVE SMILES from
     RCSB's chemical-component API — never from memory, so the
     benchmark is not circular on our own inputs,
  3. derives the search box from the crystal ligand centroid + extent,
  4. re-docks the SMILES with AutoDock Vina (the docker never sees the
     crystal pose),
  5. scores the top pose's heavy-atom RMSD against the crystal pose and
     runs PoseBusters.

Aggregate gates:
  - pose recovery: fraction with top-1 RMSD <= 2.0 Å  (criterion #1)
  - PoseBusters:   fraction with posebusters_ok         (criterion #3)

Entries that error (bad resname, unparseable ligand, fetch failure) are
reported and EXCLUDED from the gate denominators — never silently
counted as pass or fail.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

STRUCT_DIR = REPO_ROOT / "benchmarks" / "pdbbind" / "structures"

# HETATM residue names that are NOT the drug: water, ions, common
# cryo/crystallisation additives and buffers. Used to auto-pick the
# real ligand when the manifest doesn't name one.
_NON_DRUG_HET = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4",
    "1PE", "2PE", "ACT", "ACY", "FMT", "DMS", "MPD", "TRS", "EPE",
    "MES", "BME", "IMD", "TLA", "CIT", "FLC", "BOG", "NAG", "MAN",
    "BMA", "IOD", "NO3", "NH4", "AZI", "CO3", "BCT", "OXL",
    # monatomic ions
    "NA", "K", "CL", "BR", "MG", "CA", "MN", "ZN", "FE", "FE2",
    "CU", "NI", "CO", "CD", "HG", "CS", "RB", "SR", "BA", "LI",
}


def _fetch(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "cellsim-blind-bench"})
    with urllib.request.urlopen(req, timeout=timeout) as f:
        return f.read()


def fetch_pdb(pdb_id: str) -> Path:
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    out = STRUCT_DIR / f"{pdb_id.lower()}.pdb"
    if not out.exists() or out.stat().st_size == 0:
        data = _fetch(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb")
        out.write_bytes(data)
    return out


def fetch_ligand_smiles(resname: str) -> str | None:
    """Authoritative isomeric SMILES for a PDB chemical component."""
    try:
        d = json.loads(_fetch(
            f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname.upper()}"))
    except Exception:
        return None
    desc = d.get("rcsb_chem_comp_descriptor", {}) or {}
    return desc.get("SMILES_stereo") or desc.get("SMILES")


def pick_ligand_resname(pdb_path: Path) -> str | None:
    """Largest non-buffer HETATM group by heavy-atom count."""
    counts: dict[str, int] = {}
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        rn = line[17:20].strip()
        if rn in _NON_DRUG_HET:
            continue
        elem = (line[76:78].strip() or line[12:14].strip()).title()
        if elem == "H":
            continue
        counts[rn] = counts.get(rn, 0) + 1
    if not counts:
        return None
    # Prefer a drug-sized ligand (>= 6 heavy atoms); ignore tiny groups.
    drug = {k: v for k, v in counts.items() if v >= 6}
    pool = drug or counts
    return max(pool, key=pool.get)


def _center_box(coords: list[dict], pad: float = 8.0,
                lo: float = 18.0, hi: float = 30.0) -> tuple:
    xs = [c["x"] for c in coords]
    ys = [c["y"] for c in coords]
    zs = [c["z"] for c in coords]
    center = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2,
              (min(zs) + max(zs)) / 2)
    span = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    side = min(hi, max(lo, span + pad))
    return center, (side, side, side)


def run_one(entry: dict, *, exhaustiveness: int, num_modes: int,
            seed: int, cpu: int, refine_poses: bool = False) -> dict:
    from src.dock import (
        dock_ligand, attach_crystal_rmsd, attach_posebusters,
        extract_hetatm_ligand, refine_pose_openff)

    pdb_id = entry["pdb"].replace(".pdb", "")
    rec: dict = {"pdb": pdb_id, "name": entry.get("name", pdb_id)}
    t0 = time.time()
    try:
        pdb_path = fetch_pdb(pdb_id)
        resname = entry.get("resname") or pick_ligand_resname(pdb_path)
        if not resname:
            rec.update(error="no drug ligand found in structure")
            return rec
        rec["resname"] = resname
        smiles = fetch_ligand_smiles(resname)
        if not smiles:
            rec.update(error=f"no SMILES for ligand {resname}")
            return rec
        rec["smiles"] = smiles
        # ONE ligand copy for the box centre. Crystals with the ligand
        # in several chains (e.g. imatinib in 2hyy, 4 copies) would
        # otherwise centre the box on the centroid of ALL copies —
        # between the binding sites — and dock into empty space.
        coords = extract_hetatm_ligand(pdb_path, resname, first_copy_only=True)
        if len(coords) < 4:
            rec.update(error=f"ligand {resname} has < 4 heavy atoms in file")
            return rec
        center, box = _center_box(coords)
        r = dock_ligand(pdb_path, smiles, center_A=center, box_size_A=box,
                        exhaustiveness=exhaustiveness, num_modes=num_modes,
                        seed=seed, cpu=cpu)
        if not r.ok:
            rec.update(error=f"dock failed: {r.reason}")
            return rec
        # Optional OpenMM geometry refinement. NORMALLY UNNECESSARY: the
        # pose is reconstructed via Meeko's reverse conversion (correct
        # atom mapping + true Vina geometry), which passes PoseBusters
        # directly. Refinement is off by default because it is slow
        # (AM1-BCC per pose) and its restraints move the pose enough to
        # hurt the crystal RMSD. Kept opt-in for pre-FEP relaxation use.
        if refine_poses:
            try:
                r.poses = [refine_pose_openff(p, r.ligand_smiles)
                           for p in r.poses]
            except Exception as e:  # noqa: BLE001
                rec["refine_warn"] = f"{type(e).__name__}: {str(e)[:80]}"
        r = attach_crystal_rmsd(r, crystal_pdb=pdb_path, ligand_resname=resname)
        # PoseBusters for PHYSICAL validity only (criterion #3): pass the
        # receptor (clash / in-pocket / geometry tests) but NOT the
        # crystal — including the crystal folds the RMSD<=2 test into
        # posebusters_ok, which would conflate physical validity with
        # pose-recovery accuracy (criterion #1, measured separately).
        r = attach_posebusters(r, receptor_pdb=pdb_path)
        top = r.poses[0]
        top3 = [p.rmsd_vs_reference_A for p in r.poses[:3]
                if p.rmsd_vs_reference_A is not None]
        rec.update(
            n_lig_heavy=len(coords),
            dG_kcalmol=round(top.affinity_kcalmol, 2),
            top1_rmsd_A=(round(top.rmsd_vs_reference_A, 2)
                         if top.rmsd_vs_reference_A is not None else None),
            top3_best_rmsd_A=(round(min(top3), 2) if top3 else None),
            posebusters_ok=top.posebusters_ok,
        )
    except Exception as e:  # noqa: BLE001 — one bad entry must not kill the run
        rec.update(error=f"{type(e).__name__}: {str(e)[:160]}")
    finally:
        rec["wall_s"] = round(time.time() - t0, 1)
    return rec


def summarise(rows: list[dict]) -> dict:
    scored = [r for r in rows if r.get("top1_rmsd_A") is not None]
    errored = [r for r in rows if r.get("error")]
    n = len(scored)
    recov = [r for r in scored if r["top1_rmsd_A"] <= 2.0]
    recov3 = [r for r in scored
              if r.get("top3_best_rmsd_A") is not None
              and r["top3_best_rmsd_A"] <= 2.0]
    pb_scored = [r for r in scored if r.get("posebusters_ok") is not None]
    pb_ok = [r for r in pb_scored if r["posebusters_ok"]]
    return {
        "n_entries": len(rows),
        "n_scored": n,
        "n_errored": len(errored),
        # top-1 = Vina's own #1-ranked pose (scoring + sampling);
        # top-3 = best of the top 3 (sampling power, ranking-independent).
        "pose_recovery_top1_frac": (len(recov) / n) if n else None,
        "pose_recovery_top3_frac": (len(recov3) / n) if n else None,
        "pose_recovery_gate": 0.75,
        "posebusters_frac": (len(pb_ok) / len(pb_scored)) if pb_scored else None,
        "posebusters_gate": 0.95,
        "errored_pdbs": [r["pdb"] for r in errored],
    }


def main(argv=None) -> int:
    import yaml
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest", help="YAML with entries: [{pdb, resname?, name?}]")
    ap.add_argument("--out-csv", default=None)
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None,
                    help="only run the first N entries (smoke)")
    ap.add_argument("--refine", action="store_true",
                    help="opt-in OpenMM pose geometry refinement "
                         "(slow, normally unnecessary — Meeko "
                         "reconstruction already yields valid geometry)")
    args = ap.parse_args(argv)

    entries = yaml.safe_load(Path(args.manifest).read_text()).get("entries", [])
    if args.limit:
        entries = entries[:args.limit]
    print(f"[blind-dock-bench] {args.manifest}  ({len(entries)} entries)")
    rows = []
    for i, e in enumerate(entries, 1):
        r = run_one(e, exhaustiveness=args.exhaustiveness,
                    num_modes=args.num_modes, seed=args.seed, cpu=args.cpu,
                    refine_poses=args.refine)
        rows.append(r)
        if r.get("error"):
            print(f"  [{i}/{len(entries)}] {r['pdb']:6s} ERROR: {r['error']}"
                  f"  ({r['wall_s']}s)", flush=True)
        else:
            pb = {True: "PB✓", False: "PB✗", None: "PB?"}[r.get("posebusters_ok")]
            print(f"  [{i}/{len(entries)}] {r['pdb']:6s} {r.get('resname',''):4s} "
                  f"top1={r['top1_rmsd_A']}Å top3={r['top3_best_rmsd_A']}Å "
                  f"dG={r['dG_kcalmol']} {pb}  ({r['wall_s']}s)", flush=True)

    s = summarise(rows)
    print()
    print(f"[blind-dock-bench] scored {s['n_scored']}/{s['n_entries']} "
          f"({s['n_errored']} errored: {s['errored_pdbs']})")
    if s["pose_recovery_top1_frac"] is not None:
        pr = s["pose_recovery_top1_frac"]
        pr3 = s["pose_recovery_top3_frac"]
        print(f"  pose recovery top-1 <= 2.0 Å: {pr:.0%}  "
              f"[gate >= 75%] {'PASS' if pr >= 0.75 else 'FAIL'}")
        print(f"  pose recovery top-3 <= 2.0 Å: {pr3:.0%}  "
              f"(sampling power; ranking-independent)")
    if s["posebusters_frac"] is not None:
        pb = s["posebusters_frac"]
        print(f"  PoseBusters validity:           {pb:.0%}  "
              f"[gate >= 95%] {'PASS' if pb >= 0.95 else 'FAIL'}")

    if args.out_csv:
        cols = ["pdb", "name", "resname", "smiles", "n_lig_heavy",
                "dG_kcalmol", "top1_rmsd_A", "top3_best_rmsd_A",
                "posebusters_ok", "wall_s", "error"]
        outp = Path(args.out_csv)
        outp.parent.mkdir(parents=True, exist_ok=True)
        with outp.open("w", newline="", encoding="utf-8") as fo:
            w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        (outp.with_suffix(".summary.json")).write_text(json.dumps(s, indent=2))
        print(f"  wrote {outp} + {outp.with_suffix('.summary.json').name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
