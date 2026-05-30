"""Milestone A accuracy gate: hydration ΔG vs FreeSolv.

Runs compute_hydration_dg on every compound in a FreeSolv YAML,
compares to published experimental values, and reports MAE,
Pearson r, Spearman ρ. Gates on MAE ≤ 1.5 kcal/mol.

At production parameters (≥ 50 ps per window × 11 windows × 2
legs per compound) a 12-compound run takes ~hours on a GPU;
this is deliberately a workflow_dispatch / manual target, not
a PR gate.

Usage:
    cellsim fep-freesolv benchmarks/fep/freesolv_12.yaml \\
        --production-steps 25000 --n-windows 11 \\
        --out-csv /tmp/freesolv/results.csv

Set --production-steps 200 for a quick end-to-end sanity pass
(~20 s per compound × 12 = ~5 min, MAE will be large —
undersampled — but ranking is informative).
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class FreeSolvPoint:
    name: str
    smiles: str
    dG_expt_kcalmol: float
    dG_pred_kcalmol: Optional[float] = None
    uncertainty_kcalmol: Optional[float] = None
    residual_kcalmol: Optional[float] = None
    wall_seconds: Optional[float] = None
    ok: bool = False
    reason: str = ""
    # GHMC acceptance summary per compound — minimum of (vacuum
    # leg, solvent leg) for mean and per-window min. Consumed by
    # `cellsim fep-report` as a per-compound gate.
    ghmc_accept_mean: Optional[float] = None
    ghmc_accept_min: Optional[float] = None


@dataclass
class FreeSolvResult:
    yaml_path: str
    entries: list = field(default_factory=list)
    n_ok: int = 0
    n_total: int = 0
    mae_kcalmol: Optional[float] = None
    rmse_kcalmol: Optional[float] = None
    pearson_r: Optional[float] = None
    spearman_rho: Optional[float] = None
    wall_seconds: Optional[float] = None

    def summary(self) -> str:
        lines = [
            f"[FreeSolv validation]  {self.yaml_path}",
            f"  n = {self.n_total}  ok = {self.n_ok}  "
            f"wall = {(self.wall_seconds or 0):.1f} s",
        ]
        if self.mae_kcalmol is not None:
            gate = "✓ PASS" if self.mae_kcalmol <= 1.5 else "✗ FAIL"
            lines.append(
                f"  MAE = {self.mae_kcalmol:.2f} kcal/mol  "
                f"(Milestone A gate ≤ 1.5 → {gate})")
        if self.rmse_kcalmol is not None:
            lines.append(
                f"  RMSE = {self.rmse_kcalmol:.2f} kcal/mol")
        if self.pearson_r is not None:
            lines.append(
                f"  Pearson r = {self.pearson_r:+.3f}")
        if self.spearman_rho is not None:
            lines.append(
                f"  Spearman ρ = {self.spearman_rho:+.3f}")
        for p in self.entries:
            if not p.ok:
                lines.append(f"  ! {p.name:<18s} FAIL: {p.reason}")
            else:
                lines.append(
                    f"    {p.name:<18s}  expt={p.dG_expt_kcalmol:+.2f}  "
                    f"pred={p.dG_pred_kcalmol:+.2f}  "
                    f"resid={p.residual_kcalmol:+.2f}  "
                    f"({p.wall_seconds:.1f} s)")
        return "\n".join(lines)

    def as_dict(self) -> dict:
        return asdict(self)


_CSV_COLS = [
    "name", "smiles", "dG_expt_kcalmol",
    "dG_pred_kcalmol", "uncertainty_kcalmol",
    "residual_kcalmol", "ghmc_accept_mean",
    "ghmc_accept_min", "wall_seconds", "ok", "reason",
]


def _write_csv(csv_path: Path, entries: list) -> None:
    """Atomic CSV write: write to a .tmp sibling, then rename. Avoids
    leaving a half-written file if the process is killed mid-write."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=_CSV_COLS,
                           extrasaction="ignore")
        w.writeheader()
        for p in entries:
            w.writerow({k: getattr(p, k, "") for k in _CSV_COLS})
    tmp.replace(csv_path)


def _load_resume_rows(csv_path: Path) -> tuple[list, set]:
    """Load previously-completed FreeSolvPoint rows from a partial
    CSV. Returns (rows, names_completed). Used by --resume so a
    crashed multi-day CPU run picks up where it left off without
    re-doing compounds with non-empty dG_pred_kcalmol."""
    rows: list = []
    names: set = set()
    if not csv_path.exists():
        return rows, names
    with csv_path.open("r", encoding="utf-8-sig", newline="") as fi:
        reader = csv.DictReader(fi)
        for r in reader:
            pred = (r.get("dG_pred_kcalmol") or "").strip()
            if not pred or pred == "None":
                # Compound was attempted but didn't complete; re-run
                # it fresh by not adding to names_completed.
                continue
            def _f(key):
                v = (r.get(key) or "").strip()
                if not v or v == "None":
                    return None
                try:
                    return float(v)
                except ValueError:
                    return None
            p = FreeSolvPoint(
                name=r.get("name", "").strip(),
                smiles=r.get("smiles", "").strip(),
                dG_expt_kcalmol=float(
                    r.get("dG_expt_kcalmol") or "nan"),
                dG_pred_kcalmol=_f("dG_pred_kcalmol"),
                uncertainty_kcalmol=_f("uncertainty_kcalmol"),
                residual_kcalmol=_f("residual_kcalmol"),
                wall_seconds=_f("wall_seconds"),
                ok=((r.get("ok", "") or "").strip().lower()
                    in ("true", "1", "yes")),
                reason=r.get("reason", "").strip(),
                ghmc_accept_mean=_f("ghmc_accept_mean"),
                ghmc_accept_min=_f("ghmc_accept_min"),
            )
            rows.append(p)
            names.add(p.name)
    return rows, names


def run_freesolv_validation(
    yaml_path: str | Path,
    *,
    n_windows: int = 11,
    n_production_steps: int = 25000,
    n_equilibration_steps: int = 2500,
    sample_stride: int = 250,
    seed: int = 1,
    out_csv: Optional[Path] = None,
    resume: bool = False,
) -> FreeSolvResult:
    """Run the FreeSolv subset and compute MAE + correlations.

    If `out_csv` is provided, the CSV is written **after every
    compound completes** (not just at the end), so a 30+ hour CPU
    run that dies mid-way leaves all completed compounds intact for
    the next `--resume` invocation. The legacy main() also writes a
    final summary CSV when out_csv is given.

    If `resume=True` and `out_csv` exists, compounds with a non-
    empty `dG_pred_kcalmol` in that file are kept and skipped in
    the bench loop. Mid-compound crashes (no pred written) re-run
    cleanly because only completed rows are loaded.
    """
    import yaml as pyyaml
    import numpy as np

    from src.fep import compute_hydration_dg

    data = pyyaml.safe_load(Path(yaml_path).read_text())
    result = FreeSolvResult(yaml_path=str(yaml_path))
    t0 = time.time()

    completed_names: set = set()
    if resume and out_csv:
        prior_rows, completed_names = _load_resume_rows(out_csv)
        result.entries.extend(prior_rows)
        if completed_names:
            print(f"[freesolv] --resume: {len(completed_names)} "
                  f"compound(s) already in {out_csv}, skipping: "
                  f"{sorted(completed_names)}", flush=True)

    for entry in data.get("entries", []):
        if entry["name"] in completed_names:
            continue
        p = FreeSolvPoint(
            name=entry["name"],
            smiles=entry["smiles"],
            dG_expt_kcalmol=float(entry["dG_hydration_kcalmol"]))
        print(f"[freesolv] {p.name}  {p.smiles}", flush=True)
        ts = time.time()
        r = compute_hydration_dg(
            p.smiles,
            n_windows=n_windows,
            n_equilibration_steps=n_equilibration_steps,
            n_production_steps=n_production_steps,
            sample_stride=sample_stride, seed=seed)
        p.wall_seconds = time.time() - ts
        if r.ok:
            p.dG_pred_kcalmol = r.dG_hydration_kcalmol
            p.uncertainty_kcalmol = r.uncertainty_kcalmol
            p.residual_kcalmol = (
                r.dG_hydration_kcalmol - p.dG_expt_kcalmol)
            # Minimum-of-both-legs is the honest per-compound number:
            # if either leg had poor acceptance the overall ΔG is
            # compromised. Empty lists → None (legacy / non-GHMC).
            all_accepts = (
                list(r.ghmc_acceptance_vac)
                + list(r.ghmc_acceptance_solv))
            if all_accepts:
                p.ghmc_accept_mean = sum(all_accepts) / len(all_accepts)
                p.ghmc_accept_min = min(all_accepts)
            p.ok = True
            accept_tag = ""
            if p.ghmc_accept_mean is not None:
                accept_tag = (
                    f"  GHMC mean={p.ghmc_accept_mean:.0%} "
                    f"min={p.ghmc_accept_min:.0%}")
            print(f"  pred = {p.dG_pred_kcalmol:+.2f}  "
                  f"expt = {p.dG_expt_kcalmol:+.2f}  "
                  f"resid = {p.residual_kcalmol:+.2f}"
                  f"{accept_tag}",
                  flush=True)
        else:
            p.reason = r.reason
            print(f"  FAIL: {r.reason}", flush=True)
        result.entries.append(p)
        # Atomic incremental write — if a multi-day CPU run dies
        # mid-loop (lid close / kernel panic / power), every
        # completed compound survives in the CSV for --resume.
        if out_csv:
            try:
                _write_csv(out_csv, result.entries)
            except OSError as e:
                logger.warning(
                    "incremental CSV write failed: %s "
                    "(continuing run)", e)

    result.n_total = len(result.entries)
    ok = [p for p in result.entries if p.ok]
    result.n_ok = len(ok)
    if len(ok) >= 2:
        expts = np.array([p.dG_expt_kcalmol for p in ok])
        preds = np.array([p.dG_pred_kcalmol for p in ok])
        resids = preds - expts
        result.mae_kcalmol = float(np.mean(np.abs(resids)))
        result.rmse_kcalmol = float(np.sqrt(np.mean(resids ** 2)))
        if len(ok) >= 3 and expts.std() > 0 and preds.std() > 0:
            result.pearson_r = float(
                np.corrcoef(preds, expts)[0, 1])
            rx = np.argsort(np.argsort(preds))
            ry = np.argsort(np.argsort(expts))
            if rx.std() > 0 and ry.std() > 0:
                result.spearman_rho = float(
                    np.corrcoef(rx, ry)[0, 1])

    result.wall_seconds = time.time() - t0
    return result


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("yaml", type=Path,
                    help="FreeSolv YAML (e.g. benchmarks/fep/"
                         "freesolv_12.yaml)")
    ap.add_argument("--n-windows", type=int, default=11)
    ap.add_argument("--production-steps", type=int, default=25000,
                    help="Langevin steps per λ window (2 fs step; "
                         "25000 = 50 ps production at default "
                         "timestep)")
    ap.add_argument("--equilibration-steps", type=int, default=2500)
    ap.add_argument("--sample-stride", type=int, default=250)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-csv", type=Path, default=None,
                    help="write per-compound results as CSV "
                         "(written incrementally after each "
                         "compound, so a crashed run preserves "
                         "completed compounds for --resume)")
    ap.add_argument("--resume", action="store_true",
                    help="if --out-csv exists, skip compounds "
                         "with a non-empty dG_pred_kcalmol already "
                         "in it. For restarting a crashed multi-day "
                         "CPU run without losing completed work.")
    args = ap.parse_args(argv)

    r = run_freesolv_validation(
        args.yaml,
        n_windows=args.n_windows,
        n_production_steps=args.production_steps,
        n_equilibration_steps=args.equilibration_steps,
        sample_stride=args.sample_stride, seed=args.seed,
        out_csv=args.out_csv, resume=args.resume)
    print()
    print(r.summary())
    if args.out_csv:
        # Final write — same atomic helper. Idempotent with the
        # per-compound incremental writes above; this ensures the
        # last row reaches the file even if the loop body's writer
        # didn't run for the last compound (e.g. all skipped via
        # --resume).
        _write_csv(args.out_csv, r.entries)
        print(f"\nwrote {args.out_csv}")

    if r.mae_kcalmol is not None and r.mae_kcalmol <= 1.5:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
