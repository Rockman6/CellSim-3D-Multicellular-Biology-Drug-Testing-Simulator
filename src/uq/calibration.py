"""Layer 1.7 calibration harness — CellSim ΔG vs experimental ΔG.

Given a YAML of (SMILES, receptor, experimental ΔG) triples, docks
each compound, computes predicted ΔG, and reports:

    - predicted vs experimental scatter (stored for a viewer)
    - Pearson r, Spearman ρ, MAE, RMSE
    - split-conformal quantile q (from src.uq.conformal) so future
      predictions can be wrapped with calibrated CI bounds

This is the honest "does CellSim correlate with reality?" number
every non-AI docking paper is expected to cite. CellSim's mini-set
is small (4-10 entries per target); users bring their own larger
calibration sets for publication work.

Non-AI: pure physics docking + deterministic statistics.
"""

from __future__ import annotations

import logging
import math
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import dock_ligand  # noqa: E402
from src.uq.conformal import ConformalBounds  # noqa: E402

logger = logging.getLogger(__name__)


@dataclass
class CalibrationPoint:
    name: str
    smiles: str
    dG_expt_kcalmol: float
    dG_pred_kcalmol: Optional[float] = None
    residual_kcalmol: Optional[float] = None


@dataclass
class CalibrationResult:
    """Envelope for a calibration run."""

    receptor_pdb: str
    ok: bool
    reason: str = ""

    points: list = field(default_factory=list)  # list[CalibrationPoint]

    # Summary statistics (None if too few ok points)
    n_points: int = 0
    n_ok: int = 0
    pearson_r: Optional[float] = None
    spearman_rho: Optional[float] = None
    mae_kcalmol: Optional[float] = None
    rmse_kcalmol: Optional[float] = None
    conformal_q95_kcalmol: Optional[float] = None

    wall_seconds: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] calibration {self.receptor_pdb}  {self.reason}"
        lines = [
            f"[OK]   calibration  {Path(self.receptor_pdb).name}  "
            f"n_ok={self.n_ok}/{self.n_points}  "
            f"(wall {self.wall_seconds:.1f} s)"
        ]
        if self.pearson_r is not None:
            lines.append(
                f"  Pearson r      = {self.pearson_r:+.3f}")
        if self.spearman_rho is not None:
            lines.append(
                f"  Spearman ρ     = {self.spearman_rho:+.3f}")
        if self.mae_kcalmol is not None:
            lines.append(
                f"  MAE            = {self.mae_kcalmol:.2f} kcal/mol")
        if self.rmse_kcalmol is not None:
            lines.append(
                f"  RMSE           = {self.rmse_kcalmol:.2f} kcal/mol")
        if self.conformal_q95_kcalmol is not None:
            lines.append(
                f"  Conformal q95  = ±{self.conformal_q95_kcalmol:.2f} "
                "kcal/mol (95 % CI on future preds)")
        lines.append("  NOTE: small calibration sets give wide CI on "
                     "Pearson / Spearman; ~ 20+ points for "
                     "publication-grade claims.")
        return "\n".join(lines)


def _pearson(xs: list, ys: list) -> Optional[float]:
    import numpy as np
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if len(xs) < 2 or xs.std() == 0 or ys.std() == 0:
        return None
    return float(np.corrcoef(xs, ys)[0, 1])


def _spearman(xs: list, ys: list) -> Optional[float]:
    import numpy as np
    if len(xs) < 2:
        return None
    # Rank + Pearson on ranks.
    def _rank(v):
        a = np.argsort(np.argsort(v))
        return a.astype(float)
    return _pearson(_rank(np.array(xs, dtype=float)).tolist(),
                    _rank(np.array(ys, dtype=float)).tolist())


def run_calibration(
    yaml_path: str | Path,
    *,
    exhaustiveness: int = 32,
    num_modes: int = 9,
    seed: int = 1,
    cpu: int = 2,
    cache: Optional["object"] = None,
) -> CalibrationResult:
    """Dock every entry in `yaml_path`, compute correlation + conformal."""
    try:
        import yaml
    except ImportError:
        return CalibrationResult(
            receptor_pdb="?", ok=False,
            reason="PyYAML not installed (activate cellsim env)")

    import numpy as np
    import time

    data = yaml.safe_load(Path(yaml_path).read_text())
    rec = data.get("receptor", {})
    rec_pdb = str(REPO_ROOT / rec["pdb"])
    center = tuple(rec["crystal_center_A"])
    box = tuple(rec["search_box_A"])

    points: list[CalibrationPoint] = []
    t0 = time.time()

    for entry in data.get("entries", []):
        p = CalibrationPoint(
            name=entry["name"],
            smiles=entry["smiles"],
            dG_expt_kcalmol=float(entry["dG_expt_kcalmol"]),
        )
        print(f"[calib] docking {p.name}", flush=True)
        r = dock_ligand(
            rec_pdb, p.smiles,
            center_A=center, box_size_A=box,
            exhaustiveness=exhaustiveness, num_modes=num_modes,
            seed=seed, cpu=cpu, cache=cache)
        if r.ok and r.poses:
            p.dG_pred_kcalmol = float(r.poses[0].affinity_kcalmol)
            p.residual_kcalmol = (p.dG_pred_kcalmol
                                  - p.dG_expt_kcalmol)
            print(f"  pred = {p.dG_pred_kcalmol:+.2f}  "
                  f"expt = {p.dG_expt_kcalmol:+.2f}  "
                  f"resid = {p.residual_kcalmol:+.2f}")
        else:
            print(f"  FAIL: {r.reason}")
        points.append(p)

    ok_points = [p for p in points
                  if p.dG_pred_kcalmol is not None]
    result = CalibrationResult(
        receptor_pdb=rec_pdb, ok=True,
        points=points, n_points=len(points), n_ok=len(ok_points),
        wall_seconds=time.time() - t0)

    if len(ok_points) < 2:
        result.reason = (
            f"only {len(ok_points)}/{len(points)} compounds docked; "
            "cannot compute correlation")
        result.ok = False
        return result

    preds = [p.dG_pred_kcalmol for p in ok_points]
    truths = [p.dG_expt_kcalmol for p in ok_points]
    result.pearson_r = _pearson(preds, truths)
    result.spearman_rho = _spearman(preds, truths)
    resid = [abs(p.residual_kcalmol) for p in ok_points]
    result.mae_kcalmol = float(sum(resid) / len(resid))
    result.rmse_kcalmol = float(
        (sum(r * r for r in resid) / len(resid)) ** 0.5)

    # Split-conformal quantile on the calibration residuals.
    if len(ok_points) >= 2:
        cb = ConformalBounds(alpha=0.05)
        cb.calibrate(preds, truths)
        result.conformal_q95_kcalmol = cb.q

    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("yaml", help="calibration set YAML")
    ap.add_argument("--exh", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--cache", default=None,
                    help="SQLite cache path for Vina reuse")
    args = ap.parse_args()

    cache = None
    if args.cache:
        from src.cache import Cache
        cache = Cache(args.cache)

    r = run_calibration(
        args.yaml, exhaustiveness=args.exh, num_modes=args.num_modes,
        seed=args.seed, cpu=args.cpu, cache=cache)
    print()
    print(r.summary())
    sys.exit(0 if r.ok else 1)
