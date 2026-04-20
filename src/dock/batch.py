"""Multi-ligand batch screen — the wet-lab-replacement use case.

Given:
  - one receptor PDB + search box
  - a list of N ligand SMILES (with names)

Produce a ranked CSV that a biologist can hand to wet-lab:

    rank  name  smiles  formula  dG_kcalmol  Kd_nM_or_µM  pocket_ok  crystal_rmsd_A  notes
    1     compound_A  ...       -9.8  12.3 nM  True  -  top
    2     compound_B  ...       -8.6  120 nM  True  -
    ...

Plus a side PNG/HTML rank dashboard showing per-compound ΔG +
uncertainty bars, sorted by affinity.

This is what Campaign 1 Layer 1.3 actually delivers to wet-lab labs:
"given a shortlist of 100 compounds, tell me which 10 are most
worth synthesising." That's the 10–100× triage reduction claim in
MISSION.md, made concrete.

Cache-aware: if a (ligand, receptor, box, seed, exhaustiveness) tuple
has already been docked, reuse the cached pose instead of re-running
Vina. This is what makes the screen fast on iterative campaigns.

Multiprocessing: each compound docked in parallel workers. Vina
itself uses multiple CPU threads per call, but one-compound-per-
worker plus --cpu=1 per Vina is the standard efficient pattern for
docking N compounds.
"""

from __future__ import annotations

import argparse
import csv
import logging
import multiprocessing as mp
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.chem import compute_admet  # noqa: E402
from src.dock import (  # noqa: E402
    DockingResult,
    attach_crystal_rmsd,
    attach_posebusters,
    dock_ligand,
)

# src.uq.monte_carlo_dock imports src.dock internally, so it's
# imported lazily inside the worker to avoid a circular-import
# failure at module load time.

logger = logging.getLogger(__name__)


@dataclass
class BatchConfig:
    """Configuration bundle passed to each worker."""

    receptor_pdb: str
    center_A: tuple
    box_size_A: tuple
    exhaustiveness: int = 16          # fast screen default
    num_modes: int = 5
    seed: int = 1
    cpu_per_job: int = 1              # Vina threads per compound
    timeout_s: int = 600
    crystal_pdb: Optional[str] = None
    crystal_resname: Optional[str] = None
    run_posebusters: bool = True
    run_admet: bool = True            # Lipinski / QED / logP / logS
    mc_samples: int = 0               # 0 = single-seed; N>=2 = MC over
                                      # seeds [seed, seed+N) → ΔG ± CI


def _kd_label(kd_nM: float) -> str:
    if kd_nM < 1.0:
        return f"{kd_nM * 1000:.1f} pM"
    if kd_nM < 1e3:
        return f"{kd_nM:.1f} nM"
    if kd_nM < 1e6:
        return f"{kd_nM / 1e3:.1f} µM"
    return f"{kd_nM / 1e6:.1f} mM"


def _worker(task: tuple) -> dict:
    """Dock one compound and return a flat dict record."""
    name, smiles, cfg_tuple = task
    cfg = BatchConfig(**cfg_tuple)

    t0 = time.time()
    # ADMET is cheap (sub-ms); compute regardless of docking outcome
    # so a compound that fails to dock still has drug-likeness flags
    # in the output (biologist triage signal even on dock-fails).
    admet = None
    if cfg.run_admet:
        try:
            admet = compute_admet(smiles)
        except Exception as e:
            logger.debug("ADMET failed on %s: %s", name, e)

    # If MC sampling requested, run multiple seeds up front to
    # get the ΔG distribution. The pose geometry is taken from the
    # *best* seed (the one that found the lowest ΔG), which is what
    # biologists intuitively want to see.
    mc_stats = None
    dock_seed = cfg.seed
    if cfg.mc_samples >= 2:
        try:
            from src.uq import monte_carlo_dock  # lazy: see module top
            mc_stats = monte_carlo_dock(
                receptor_pdb=cfg.receptor_pdb,
                ligand_smiles=smiles,
                center_A=cfg.center_A, box_size_A=cfg.box_size_A,
                n_samples=cfg.mc_samples, base_seed=cfg.seed,
                exhaustiveness=cfg.exhaustiveness,
                num_modes=cfg.num_modes,
                cpu_per_job=cfg.cpu_per_job, workers=1,
                timeout_s=cfg.timeout_s)
            if mc_stats.ok and mc_stats.best_seed is not None:
                dock_seed = mc_stats.best_seed
        except Exception as e:
            logger.debug("MC failed on %s: %s", name, e)

    try:
        r = dock_ligand(
            cfg.receptor_pdb, smiles,
            center_A=cfg.center_A, box_size_A=cfg.box_size_A,
            exhaustiveness=cfg.exhaustiveness,
            num_modes=cfg.num_modes, seed=dock_seed,
            cpu=cfg.cpu_per_job, timeout_s=cfg.timeout_s)
        if r.ok and cfg.crystal_pdb and cfg.crystal_resname:
            try:
                r = attach_crystal_rmsd(
                    r, crystal_pdb=cfg.crystal_pdb,
                    ligand_resname=cfg.crystal_resname)
            except Exception as e:
                logger.debug("crystal RMSD failed for %s: %s", name, e)
        if r.ok and cfg.run_posebusters:
            try:
                r = attach_posebusters(
                    r, receptor_pdb=cfg.receptor_pdb,
                    crystal_pdb=cfg.crystal_pdb,
                    crystal_resname=cfg.crystal_resname)
            except Exception as e:
                logger.debug("PoseBusters failed for %s: %s", name, e)
    except Exception as e:
        logger.exception("docking crashed on %s", name)
        return dict(name=name, smiles=smiles, ok=False,
                    reason=f"crashed: {e}",
                    wall_s=time.time() - t0)

    wall = time.time() - t0

    def _admet_fields() -> dict:
        if admet is None or not admet.ok:
            return {}
        return dict(
            MW=round(admet.MW, 2),
            logP=round(admet.logP, 3),
            TPSA=round(admet.tpsa, 1),
            HBA=admet.hba, HBD=admet.hbd, rotb=admet.rotb,
            ro5_pass=admet.ro5_pass,
            ro5_violations=admet.ro5_violations,
            QED=round(admet.qed, 3) if admet.qed is not None else None,
            logS=round(admet.logS_ESOL, 2),
            solubility=admet.solubility_class,
        )

    def _mc_fields() -> dict:
        if mc_stats is None or not mc_stats.ok:
            return {}
        return dict(
            dG_mean_kcalmol=round(mc_stats.dG_mean_kcalmol, 3),
            dG_std_kcalmol=round(mc_stats.dG_std_kcalmol, 3),
            dG_ci95_lo=round(mc_stats.dG_ci95_lo_kcalmol, 3),
            dG_ci95_hi=round(mc_stats.dG_ci95_hi_kcalmol, 3),
            mc_n_ok=mc_stats.n_ok,
            mc_n_samples=mc_stats.n_samples,
        )

    if not r.ok:
        return dict(name=name, smiles=smiles, ok=False,
                    reason=r.reason, wall_s=wall,
                    **_admet_fields(), **_mc_fields())

    top = r.poses[0]
    record = dict(
        name=name, smiles=smiles, ok=True, reason="",
        formula=r.ligand_formula or "",
        inchi_key=r.ligand_inchi_key or "",
        dG_kcalmol=round(top.affinity_kcalmol, 3),
        dG_kJmol=round(top.affinity_kJmol, 3),
        Kd_nM=round(top.kd_implied_nM, 3),
        Kd_human=_kd_label(top.kd_implied_nM),
        crystal_rmsd_A=(round(top.rmsd_vs_reference_A, 3)
                         if top.rmsd_vs_reference_A is not None
                         else None),
        pocket_ok=top.posebusters_pocket_ok,
        geometry_ok=top.posebusters_geometry_ok,
        pb_all_ok=top.posebusters_ok,
        n_poses=len(r.poses),
        seed=r.seed, exhaustiveness=r.exhaustiveness,
        wall_s=round(wall, 1),
        **_admet_fields(),
        **_mc_fields(),
    )
    return record


def load_smi(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split(None, 1)
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else smi
        rows.append((smi, name))
    return rows


def run_batch(
    smi_file: Path, cfg: BatchConfig, *,
    workers: int = 1, out_csv: Optional[Path] = None,
) -> list[dict]:
    """Dock every SMILES in `smi_file` against `cfg.receptor_pdb`.

    Returns records sorted by ΔG ascending (best affinity first).
    If `out_csv` is provided, also writes the sorted CSV.
    """
    entries = load_smi(smi_file)
    n = len(entries)
    print(f"[batch] {n} compounds -> {Path(cfg.receptor_pdb).name}  "
          f"workers={workers} exh={cfg.exhaustiveness}", flush=True)

    cfg_tuple = cfg.__dict__
    tasks = [(name, smi, cfg_tuple) for smi, name in entries]

    records: list[dict] = []
    t0 = time.time()
    if workers <= 1:
        for i, task in enumerate(tasks, 1):
            rec = _worker(task)
            records.append(rec)
            _progress_line(i, n, rec)
    else:
        with mp.Pool(workers) as pool:
            for i, rec in enumerate(
                    pool.imap_unordered(_worker, tasks, chunksize=1), 1):
                records.append(rec)
                _progress_line(i, n, rec)
    dt = time.time() - t0
    print(f"[batch] wall={dt:.1f} s  "
          f"({n/dt:.2f} compounds/s aggregate)", flush=True)

    # Sort: successful first, ranked by ΔG ascending.
    ok = sorted([r for r in records if r.get("ok")],
                key=lambda r: r["dG_kcalmol"])
    fail = [r for r in records if not r.get("ok")]
    sorted_records = ok + fail
    # Assign rank within successful set.
    for rank, r in enumerate(ok, 1):
        r["rank"] = rank
    for r in fail:
        r["rank"] = None

    if out_csv:
        _write_csv(sorted_records, out_csv)
        print(f"[batch] wrote {out_csv}", flush=True)

    return sorted_records


def _progress_line(i: int, n: int, rec: dict) -> None:
    if rec.get("ok"):
        tag = f"ΔG={rec['dG_kcalmol']:>6.2f} Kd={rec['Kd_human']:<10s}"
        if rec.get("dG_mean_kcalmol") is not None:
            tag += (f" MC={rec['dG_mean_kcalmol']:+.2f}"
                    f"±{rec['dG_std_kcalmol']:.2f}")
        pb = ("pocket:ok" if rec.get("pocket_ok") else
              "pocket:fail" if rec.get("pocket_ok") is False else
              "pocket:n/a")
    else:
        tag = "FAIL"
        pb = (rec.get("reason") or "?")[:60]
    print(f"  [{i:>3d}/{n}]  {rec['name']:<25s}  {tag}  {pb}",
          flush=True)


def _write_csv(records: list[dict], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["rank", "name", "smiles", "formula", "inchi_key",
            "dG_kcalmol", "dG_kJmol", "Kd_nM", "Kd_human",
            "dG_mean_kcalmol", "dG_std_kcalmol",
            "dG_ci95_lo", "dG_ci95_hi",
            "mc_n_ok", "mc_n_samples",
            "crystal_rmsd_A", "pocket_ok", "geometry_ok", "pb_all_ok",
            "MW", "logP", "TPSA", "HBA", "HBD", "rotb",
            "ro5_pass", "ro5_violations", "QED", "logS", "solubility",
            "n_poses", "seed", "exhaustiveness", "wall_s",
            "ok", "reason"]
    with out.open("w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for rec in records:
            w.writerow(rec)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smi", type=Path, required=True,
                    help="input SMILES file (tab or whitespace; "
                         "'smiles [name]' per line)")
    ap.add_argument("--receptor", required=True,
                    help="receptor PDB")
    ap.add_argument("--center", required=True,
                    help="comma-separated x,y,z Å")
    ap.add_argument("--box", default="22,22,22",
                    help="comma-separated dx,dy,dz Å")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--num-modes", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--cpu-per-job", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--crystal-pdb", default=None)
    ap.add_argument("--crystal-resname", default=None)
    ap.add_argument("--no-posebusters", action="store_true")
    ap.add_argument("--no-admet", action="store_true",
                    help="skip ADMET descriptors (Lipinski / QED / "
                         "logP / logS) in the output")
    ap.add_argument("--mc", type=int, default=0, metavar="N",
                    help="run N Monte-Carlo seeds per compound "
                         "(adds ΔG mean ± SD + 95%% CI per row; "
                         "0 = disable, >=2 = enable)")
    ap.add_argument("--out-csv", type=Path, default=None,
                    help="output CSV path (default: print to stdout)")
    args = ap.parse_args(argv)

    cfg = BatchConfig(
        receptor_pdb=args.receptor,
        center_A=tuple(float(x) for x in args.center.split(",")),
        box_size_A=tuple(float(x) for x in args.box.split(",")),
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        seed=args.seed, cpu_per_job=args.cpu_per_job,
        timeout_s=args.timeout,
        crystal_pdb=args.crystal_pdb,
        crystal_resname=args.crystal_resname,
        run_posebusters=not args.no_posebusters,
        run_admet=not args.no_admet,
        mc_samples=args.mc,
    )
    records = run_batch(
        args.smi, cfg, workers=args.workers, out_csv=args.out_csv)

    # Biologist-friendly summary table on stdout.
    print()
    print("RANK  NAME                       ΔG(kcal)   K_d        "
          "POCKET  RMSD      Ro5  QED   logS   class")
    print("-" * 104)
    for r in records[:20]:
        if r.get("ok"):
            rmsd = (f"{r['crystal_rmsd_A']:.2f} Å"
                    if r.get("crystal_rmsd_A") is not None else "-")
            pocket = ("✓" if r.get("pocket_ok")
                      else "✗" if r.get("pocket_ok") is False
                      else "?")
            ro5 = ("✓" if r.get("ro5_pass") else
                   f"✗×{r.get('ro5_violations') or 0}")
            qed = (f"{r['QED']:.2f}" if r.get("QED") is not None
                   else "-")
            logs = (f"{r['logS']:+.2f}" if r.get("logS") is not None
                    else "-")
            sol = (r.get("solubility") or "-")[:15]
            print(f"{r['rank']:>4d}  {r['name']:<25s}  "
                  f"{r['dG_kcalmol']:>8.2f}   {r['Kd_human']:<9s}  "
                  f" {pocket}      {rmsd:<8s}  {ro5:<4s} "
                  f"{qed:<5s} {logs:<6s} {sol}")
        else:
            print(f" -    {r['name']:<25s}  FAIL  "
                  f"{(r.get('reason') or '')[:50]}")
    if len(records) > 20:
        print(f"... ({len(records) - 20} more; see CSV)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
