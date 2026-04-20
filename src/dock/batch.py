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
    refine_poses: bool = False        # OpenMM ligand-only minimise
                                      # → fixes bond / angle / clash
                                      # geometry so PB geometry_ok
                                      # flips True on most poses
    cache_path: Optional[str] = None  # if set, content-addressed
                                      # Vina cache at that SQLite path;
                                      # second run of an identical
                                      # (ligand, receptor, box, exh,
                                      # seed) returns instantly
    export_poses_dir: Optional[str] = None  # if set, per-compound
                                            # SDF + PDB files of all
                                            # poses get written here
    run_strain: bool = True           # post-dock UFF-ensemble strain
                                      # diagnostic on top-pose via
                                      # src.dock.strain. Surfaces
                                      # Vina artefacts (high-score but
                                      # conformationally non-physical
                                      # poses). ~0.5-2 s per compound
                                      # at ensemble_n=20.
    strain_ensemble_n: int = 20       # conformers for UFF ensemble;
                                      # PoseBusters default is 100
                                      # (slower, tighter ratio CI)
    strain_gate: bool = True          # if top-1 pose is strain-
                                      # rejected, fall through to
                                      # the next good/acceptable
                                      # pose. Keeps the reported ΔG
                                      # physically trustworthy.


def _triage_call(record: dict) -> tuple[str, str]:
    """Synthesise a per-compound triage verdict for a wet-lab user.

    Combines the ΔG tier with the pose-trust signals (strain band,
    PoseBusters pocket flag) and the critical ADMET safety flags
    (mutagenicity and hERG) into one actionable call:

        follow_up     ΔG strong, pose trustworthy, ADMET clean
        review        one signal is ambiguous — needs human judgement
        deprioritise  ΔG weak or pose questionable, but not disqualified
        drop          ΔG too weak OR pose non-physical OR known-bad ADMET

    Returns `(call, reason)`. The reason is a short human-readable
    string biologists can paste into lab notes.
    """
    dG = record.get("dG_kcalmol")
    if dG is None:
        return "drop", "no ΔG"

    strain = record.get("strain_band")
    pocket_ok = record.get("pocket_ok")
    mut = record.get("mutagenic_risk")
    herg = record.get("herg_risk")
    ro5 = record.get("ro5_violations")

    # Drop criteria — any one is disqualifying.
    if strain == "reject":
        return "drop", "non-physical pose (strain)"
    if mut == "high":
        return "drop", "high mutagenicity (Ames SMARTS hit)"
    if dG > -6.0:
        return "drop", "ΔG > -6 (too weak to triage)"
    if ro5 is not None and ro5 >= 3:
        return "drop", f"{ro5} Ro5 violations"

    # Review flags (ambiguity).
    reasons: list[str] = []
    if strain == "suspicious":
        reasons.append("suspicious pose strain")
    if pocket_ok is False:
        reasons.append("PoseBusters pocket-fail")
    if herg == "high":
        reasons.append("hERG alert")
    if mut == "medium":
        reasons.append("medium mutagenicity")
    if reasons:
        # Weak ΔG + any review flag → deprioritise; strong ΔG with
        # a flag is worth a biologist's time → review.
        if dG > -7.5:
            return "deprioritise", "; ".join(reasons)
        return "review", "; ".join(reasons)

    # Deprioritise vs follow_up purely on ΔG tier.
    if dG > -7.5:
        return "deprioritise", f"ΔG {dG:+.2f} borderline"
    return "follow_up", f"ΔG {dG:+.2f}, pose trustworthy"


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

    # Open the cache once per worker (SQLite is cheap and one
    # connection per worker is the canonical pattern).
    _cache_obj = None
    if cfg.cache_path:
        try:
            from src.cache import Cache
            _cache_obj = Cache(cfg.cache_path)
        except Exception as e:
            logger.debug("cache open failed: %s", e)

    try:
        r = dock_ligand(
            cfg.receptor_pdb, smiles,
            center_A=cfg.center_A, box_size_A=cfg.box_size_A,
            exhaustiveness=cfg.exhaustiveness,
            num_modes=cfg.num_modes, seed=dock_seed,
            cpu=cfg.cpu_per_job, timeout_s=cfg.timeout_s,
            cache=_cache_obj)
        if r.ok and cfg.crystal_pdb and cfg.crystal_resname:
            try:
                r = attach_crystal_rmsd(
                    r, crystal_pdb=cfg.crystal_pdb,
                    ligand_resname=cfg.crystal_resname)
            except Exception as e:
                logger.debug("crystal RMSD failed for %s: %s", name, e)
        # Optional OpenMM refinement step (fixes Vina's approximate
        # geometry — runs BEFORE PoseBusters so the refined geometry
        # is what gets graded.
        if r.ok and cfg.refine_poses:
            try:
                from src.dock import refine_pose_openff
                r.poses = [refine_pose_openff(p, r.ligand_smiles)
                           for p in r.poses]
            except Exception as e:
                logger.debug("refine failed for %s: %s", name, e)

        if r.ok and cfg.run_posebusters:
            try:
                r = attach_posebusters(
                    r, receptor_pdb=cfg.receptor_pdb,
                    crystal_pdb=cfg.crystal_pdb,
                    crystal_resname=cfg.crystal_resname)
            except Exception as e:
                logger.debug("PoseBusters failed for %s: %s", name, e)

        # Optional per-compound pose export (SDF + PDB).
        if r.ok and cfg.export_poses_dir:
            try:
                from src.dock import export_poses_sdf, export_poses_pdb
                export_dir = Path(cfg.export_poses_dir)
                export_dir.mkdir(parents=True, exist_ok=True)
                safe_name = "".join(c if c.isalnum() or c in "-_."
                                     else "_" for c in name)
                sdf_path = export_dir / f"{safe_name}.sdf"
                pdb_path = export_dir / f"{safe_name}.pdb"
                export_poses_sdf(r, sdf_path)
                export_poses_pdb(r, pdb_path)
            except Exception as e:
                logger.debug("pose export failed for %s: %s", name, e)
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
            bbb_permeable=admet.bbb_permeable,
            herg_risk=admet.herg_risk,
            herg_alerts=(",".join(admet.herg_alerts)
                         if admet.herg_alerts else ""),
            mutagenic_risk=admet.mutagenic_risk,
            mutagenic_alerts=(",".join(admet.mutagenic_alerts)
                              if admet.mutagenic_alerts else ""),
        )

    strain_result = None
    strain_promoted_from: Optional[int] = None  # 1-indexed pose rank
    if cfg.run_strain and r.ok and r.poses:
        try:
            from src.dock.strain import ligand_strain
            # Strain-aware top-pose selection. Score every pose;
            # if the best-ΔG pose is strain-rejected, promote the
            # best-ΔG pose that passes strain (band good /
            # acceptable). Tag the record with which rank was
            # promoted so biologists see that it isn't the raw
            # Vina #1. If NO pose passes, keep pose #1 (current
            # behaviour — honesty: Vina couldn't find any
            # physically-reasonable pose for this compound).
            pose_strains = []
            for i, pose in enumerate(r.poses):
                s = ligand_strain(
                    pose.elements, pose.positions_A,
                    r.ligand_smiles,
                    ensemble_n=cfg.strain_ensemble_n)
                pose_strains.append((i, pose, s))

            top_i, top_pose, top_strain = pose_strains[0]
            if (cfg.strain_gate and top_strain.ok
                    and top_strain.band == "reject"):
                for i, pose, s in pose_strains:
                    if s.ok and s.band in ("good", "acceptable"):
                        top_i, top_pose, top_strain = i, pose, s
                        strain_promoted_from = i + 1
                        # Rewrite r.poses so downstream consumers
                        # (CSV dG_kcalmol, strain_* fields, SDF
                        # export) see the promoted pose as rank 1.
                        r.poses[0] = pose
                        break
            strain_result = top_strain
        except Exception as e:
            logger.debug("strain failed for %s: %s", name, e)

    def _strain_fields() -> dict:
        if strain_result is None or not strain_result.ok:
            return {}
        out = dict(
            strain_band=strain_result.band,
            strain_kcalmol=round(
                strain_result.strain_kcalmol, 2),
            strain_ratio=round(strain_result.energy_ratio, 3),
        )
        if strain_promoted_from is not None:
            out["strain_promoted_from_rank"] = strain_promoted_from
        return out

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
        **_strain_fields(),
    )
    call, why = _triage_call(record)
    record["triage"] = call
    record["triage_reason"] = why
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
        # Side-effect: also render the triage-breakdown PNG next to
        # the CSV so biologists get the dashboard without running a
        # separate subcommand. Silently skip if matplotlib or the
        # viewer module isn't importable.
        try:
            from src.dock.triage_viewer import render_triage_dashboard
            png = Path(out_csv).with_suffix(".triage.png")
            render_triage_dashboard(
                [{k: (str(v) if v is not None else "")
                  for k, v in r.items()} for r in sorted_records],
                png,
                title=f"{Path(cfg.receptor_pdb).stem} — "
                      f"{sum(1 for r in sorted_records if r.get('ok'))}"
                      " compounds")
            print(f"[batch] wrote {png}", flush=True)
        except Exception as e:
            logger.debug("triage PNG render skipped: %s", e)

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
        if rec.get("strain_band"):
            tag += f"  strain:{rec['strain_band']}"
        if rec.get("strain_promoted_from_rank"):
            tag += (f"  (promoted pose #"
                    f"{rec['strain_promoted_from_rank']})")
    else:
        tag = "FAIL"
        pb = (rec.get("reason") or "?")[:60]
    print(f"  [{i:>3d}/{n}]  {rec['name']:<25s}  {tag}  {pb}",
          flush=True)


def _write_csv(records: list[dict], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["rank", "name", "triage", "triage_reason",
            "smiles", "formula", "inchi_key",
            "dG_kcalmol", "dG_kJmol", "Kd_nM", "Kd_human",
            "dG_mean_kcalmol", "dG_std_kcalmol",
            "dG_ci95_lo", "dG_ci95_hi",
            "mc_n_ok", "mc_n_samples",
            "crystal_rmsd_A", "pocket_ok", "geometry_ok", "pb_all_ok",
            "strain_band", "strain_kcalmol", "strain_ratio",
            "strain_promoted_from_rank",
            "MW", "logP", "TPSA", "HBA", "HBD", "rotb",
            "ro5_pass", "ro5_violations", "QED", "logS", "solubility",
            "bbb_permeable", "herg_risk", "herg_alerts",
            "mutagenic_risk", "mutagenic_alerts",
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
    ap.add_argument("--center", default=None,
                    help="comma-separated x,y,z Å. Omit together "
                         "with --box to let fpocket auto-detect the "
                         "top druggable pocket (biologist-friendly "
                         "default for unknown targets).")
    ap.add_argument("--box", default=None,
                    help="comma-separated dx,dy,dz Å. If omitted and "
                         "--center is omitted, uses fpocket's "
                         "suggested box.")
    ap.add_argument("--auto-pocket-rank", type=int, default=1,
                    help="which fpocket-detected pocket to target "
                         "when --center is not supplied (1 = best "
                         "drug score; default 1)")
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
    ap.add_argument("--profile-top-k", type=int, default=0,
                    metavar="K",
                    help="after ranking, run src.chem.profile on "
                         "the top K successful compounds and save "
                         "one PNG per compound next to --out-csv. "
                         "0 = skip (default)")
    ap.add_argument("--refine-poses", action="store_true",
                    help="post-dock OpenMM ligand-only minimisation "
                         "to fix Vina's approximate bond / angle / "
                         "clash geometry. Flips PoseBusters "
                         "geometry_ok True on most poses; required "
                         "before downstream FEP. Adds ~1-2 s per "
                         "pose.")
    ap.add_argument("--cache", type=str, default=None, metavar="PATH",
                    help="SQLite cache path (content-addressed). "
                         "Same (ligand, receptor, box, exh, seed) "
                         "returns from cache on repeat runs "
                         "(~1000x speedup). Use one cache per project.")
    ap.add_argument("--export-poses-dir", type=str, default=None,
                    metavar="DIR",
                    help="per-compound SDF + PDB files of all poses "
                         "written here. Open in PyMOL / ChimeraX / "
                         "any SDF viewer.")
    ap.add_argument("--no-strain", action="store_true",
                    help="skip the UFF-ensemble ligand-strain check "
                         "(default: on). Strain flags Vina poses "
                         "that score well but are conformationally "
                         "non-physical; adds ~1 s per compound.")
    ap.add_argument("--no-strain-gate", action="store_true",
                    help="disable the strain-aware top-pose "
                         "promotion (default: on). When on, a "
                         "strain-rejected top pose is replaced by "
                         "the best-ΔG pose whose strain is good/"
                         "acceptable; the reported ΔG stays "
                         "physically trustworthy.")
    ap.add_argument("--strain-ensemble", type=int, default=20,
                    metavar="N",
                    help="number of conformers for the UFF-ensemble "
                         "strain baseline (default 20; PoseBusters "
                         "default 100 is tighter but slower).")
    ap.add_argument("--out-csv", type=Path, default=None,
                    help="output CSV path (default: print to stdout)")
    args = ap.parse_args(argv)

    # Resolve search box: either CLI args or fpocket auto-detect.
    if args.center is None:
        from src.dock import detect_pockets
        print(f"[batch] --center not given; running fpocket on "
              f"{args.receptor} to auto-detect binding pocket "
              f"(rank={args.auto_pocket_rank})", flush=True)
        pockets = detect_pockets(args.receptor,
                                  top_k=max(5, args.auto_pocket_rank))
        if not pockets or len(pockets) < args.auto_pocket_rank:
            print(f"[batch] fpocket found {len(pockets)} pockets; "
                  f"cannot satisfy --auto-pocket-rank="
                  f"{args.auto_pocket_rank}. Aborting.",
                  file=sys.stderr)
            return 2
        pkt = pockets[args.auto_pocket_rank - 1]
        auto_center = pkt.center_A
        auto_box = pkt.suggested_box_A
        print(f"[batch] auto-pocket #{pkt.rank}: "
              f"drug={pkt.drug_score:.3f}  "
              f"volume={pkt.volume_A3:.0f} Å³  "
              f"center={auto_center}  box={auto_box}", flush=True)
    else:
        auto_center = tuple(float(x) for x in args.center.split(","))
        auto_box = (tuple(float(x) for x in args.box.split(","))
                    if args.box else (22.0, 22.0, 22.0))

    cfg = BatchConfig(
        receptor_pdb=args.receptor,
        center_A=auto_center,
        box_size_A=auto_box,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        seed=args.seed, cpu_per_job=args.cpu_per_job,
        timeout_s=args.timeout,
        crystal_pdb=args.crystal_pdb,
        crystal_resname=args.crystal_resname,
        run_posebusters=not args.no_posebusters,
        run_admet=not args.no_admet,
        mc_samples=args.mc,
        refine_poses=args.refine_poses,
        cache_path=args.cache,
        export_poses_dir=args.export_poses_dir,
        run_strain=not args.no_strain,
        strain_ensemble_n=args.strain_ensemble,
        strain_gate=not args.no_strain_gate,
    )
    records = run_batch(
        args.smi, cfg, workers=args.workers, out_csv=args.out_csv)

    # Optionally produce a per-compound drug profile dashboard for
    # the top-K successful hits (ties chem + quantum + ADMET + dock
    # into a one-page biologist summary).
    if args.profile_top_k > 0:
        from src.chem.profile import render_profile
        top_ok = [r for r in records if r.get("ok")][:args.profile_top_k]
        if args.out_csv is not None:
            out_dir = Path(args.out_csv).resolve().parent
        else:
            out_dir = Path(".").resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n[batch] generating profile PNGs for top "
              f"{len(top_ok)} hits in {out_dir}", flush=True)
        for rec in top_ok:
            png = out_dir / f"profile_{rec['rank']:02d}_{rec['name']}.png"
            try:
                render_profile(
                    rec["name"], rec["smiles"],
                    save=str(png), show=False,
                    seed=args.seed)
            except SystemExit as e:
                print(f"  profile {rec['name']} skipped: {e}",
                      flush=True)
            except Exception as e:
                print(f"  profile {rec['name']} failed: {e}",
                      flush=True)

    # Biologist-friendly summary table on stdout.
    print()
    print("RANK  NAME                       TRIAGE        "
          "ΔG(kcal)   K_d        "
          "POCKET  STRAIN       Ro5  QED   logS")
    print("-" * 112)
    for r in records[:20]:
        if r.get("ok"):
            pocket = ("✓" if r.get("pocket_ok")
                      else "✗" if r.get("pocket_ok") is False
                      else "?")
            ro5 = ("✓" if r.get("ro5_pass") else
                   f"✗×{r.get('ro5_violations') or 0}")
            qed = (f"{r['QED']:.2f}" if r.get("QED") is not None
                   else "-")
            logs = (f"{r['logS']:+.2f}" if r.get("logS") is not None
                    else "-")
            triage = r.get("triage") or "-"
            strain = r.get("strain_band") or "-"
            print(f"{r['rank']:>4d}  {r['name']:<25s}  "
                  f"{triage:<12s}  "
                  f"{r['dG_kcalmol']:>8.2f}   {r['Kd_human']:<9s}  "
                  f" {pocket}      {strain:<11s} "
                  f"{ro5:<4s} {qed:<5s} {logs:<6s}")
        else:
            print(f" -    {r['name']:<25s}  FAIL  "
                  f"{(r.get('reason') or '')[:50]}")
    if len(records) > 20:
        print(f"... ({len(records) - 20} more; see CSV)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
