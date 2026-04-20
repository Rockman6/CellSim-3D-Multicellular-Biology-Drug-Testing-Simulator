"""Monte-Carlo over Vina docking seeds → ΔG ± CI.

AutoDock Vina's stochastic Monte-Carlo search gives slightly
different results on different random seeds. The field convention
is to accept a single-seed result and hope it's close to the
population mean, but this loses honest uncertainty information.

This module does N independent Vina runs with distinct seeds and
reports the sampling distribution of the top-pose ΔG: mean, SD,
and a 95 % CI from the sample percentiles. Biologists now see
"ΔG = -7.4 ± 0.3 kcal/mol" instead of a naked -7.4 that hides
its own variability.

This is the Campaign-1 Layer-1.6 primary UQ path per MISSION.md:
mechanistic / statistical sampling, no neural ensembles, no
learned bounds. MAPIE conformal prediction is a separate follow-up
layer (post-hoc non-parametric wrapper).

Usage:
    from src.uq import monte_carlo_dock
    r = monte_carlo_dock(
        receptor_pdb='benchmarks/dock/1stp.pdb',
        ligand_smiles='OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12',
        center_A=(11.12, 1.68, -10.75), box_size_A=(20, 20, 20),
        n_samples=8, base_seed=1,
        exhaustiveness=16, num_modes=3, workers=4)
    print(r.summary())

Honest disclosure: Vina's stochastic run-to-run jitter is a LOWER
bound on the real uncertainty. True uncertainty also includes the
choice of force field, charges, protonation, pose refinement, etc.
Sobol sensitivity over those knobs is a separate follow-up PR.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import dock_ligand  # noqa: E402

logger = logging.getLogger(__name__)


@dataclass
class DockingMCResult:
    """Monte-Carlo envelope over N seeds."""

    ligand_smiles: str
    receptor_pdb: str
    ok: bool
    reason: str = ""

    n_samples: int = 0
    n_ok: int = 0
    base_params: dict = field(default_factory=dict)
    seeds: list = field(default_factory=list)

    # Per-sample top-pose ΔG values (only for ok samples).
    dG_samples_kcalmol: list = field(default_factory=list)

    # Aggregate statistics (only when n_ok >= 2).
    dG_mean_kcalmol: Optional[float] = None
    dG_std_kcalmol: Optional[float] = None
    dG_ci95_lo_kcalmol: Optional[float] = None
    dG_ci95_hi_kcalmol: Optional[float] = None

    # Honest best ΔG seen across samples and its seed (for the
    # biologist who wants the single-best pose for inspection).
    best_dG_kcalmol: Optional[float] = None
    best_seed: Optional[int] = None

    wall_seconds: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return (f"[FAIL] MC dock {self.ligand_smiles}  "
                    f"{self.reason}")
        if self.dG_mean_kcalmol is None:
            return (f"[OK]   MC dock {self.ligand_smiles}  "
                    f"n_ok={self.n_ok}/{self.n_samples}  (too few "
                    "for statistics)")
        return (f"[OK]   MC dock {self.ligand_smiles}  "
                f"n_ok={self.n_ok}/{self.n_samples}  "
                f"ΔG = {self.dG_mean_kcalmol:+.2f} "
                f"± {self.dG_std_kcalmol:.2f} kcal/mol  "
                f"95 % CI [{self.dG_ci95_lo_kcalmol:+.2f}, "
                f"{self.dG_ci95_hi_kcalmol:+.2f}]  "
                f"best={self.best_dG_kcalmol:+.2f}@seed{self.best_seed}  "
                f"({self.wall_seconds:.1f} s)")


def _one_seed(task: tuple) -> tuple:
    """Worker: one Vina run at a specific seed."""
    (receptor_pdb, smiles, center_A, box_size_A,
     exhaustiveness, num_modes, cpu, timeout_s, seed) = task
    try:
        r = dock_ligand(
            receptor_pdb, smiles,
            center_A=center_A, box_size_A=box_size_A,
            exhaustiveness=exhaustiveness, num_modes=num_modes,
            seed=seed, cpu=cpu, timeout_s=timeout_s)
    except Exception as e:
        return seed, None, f"crash: {e}"
    if not r.ok or not r.poses:
        return seed, None, r.reason or "no poses"
    return seed, float(r.poses[0].affinity_kcalmol), ""


def monte_carlo_dock(
    *,
    receptor_pdb: str,
    ligand_smiles: str,
    center_A: tuple,
    box_size_A: tuple = (22.0, 22.0, 22.0),
    n_samples: int = 8,
    base_seed: int = 1,
    exhaustiveness: int = 16,
    num_modes: int = 3,
    cpu_per_job: int = 1,
    workers: int = 1,
    timeout_s: int = 600,
) -> DockingMCResult:
    """Run N Vina dockings at seeds [base_seed, base_seed + N) and
    report the top-pose ΔG sampling distribution.

    Returns a DockingMCResult with `ok=True` iff at least 2 of N
    runs produced a valid top pose (needed for a SD). Seeds are
    consecutive integers starting from `base_seed` so the MC is
    deterministic given (base_seed, n_samples).
    """
    import numpy as np

    seeds = list(range(int(base_seed), int(base_seed) + int(n_samples)))
    tasks = [
        (receptor_pdb, ligand_smiles, tuple(center_A), tuple(box_size_A),
         exhaustiveness, num_modes, cpu_per_job, timeout_s, s)
        for s in seeds
    ]

    result = DockingMCResult(
        ligand_smiles=ligand_smiles,
        receptor_pdb=receptor_pdb,
        ok=False, n_samples=n_samples,
        base_params={
            "center_A": tuple(center_A),
            "box_size_A": tuple(box_size_A),
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "cpu_per_job": cpu_per_job,
            "timeout_s": timeout_s,
            "base_seed": base_seed,
        },
        seeds=seeds,
    )

    t0 = time.time()
    if workers <= 1:
        outcomes = [_one_seed(t) for t in tasks]
    else:
        with mp.Pool(workers) as pool:
            outcomes = list(pool.imap_unordered(_one_seed, tasks))

    dG_values: list[float] = []
    best_dG = None
    best_seed = None
    for seed, dG, err in outcomes:
        if dG is None:
            logger.debug("seed %d failed: %s", seed, err)
            continue
        dG_values.append(float(dG))
        if best_dG is None or dG < best_dG:
            best_dG = float(dG)
            best_seed = int(seed)

    result.n_ok = len(dG_values)
    result.dG_samples_kcalmol = dG_values
    result.best_dG_kcalmol = best_dG
    result.best_seed = best_seed
    result.wall_seconds = time.time() - t0

    if len(dG_values) < 2:
        result.reason = (f"only {len(dG_values)}/{n_samples} seeds "
                         "produced a valid pose")
        return result

    arr = np.array(dG_values)
    result.dG_mean_kcalmol = float(arr.mean())
    result.dG_std_kcalmol = float(arr.std(ddof=1))
    # 95 % CI via the 2.5-97.5 percentile (no parametric assumption).
    lo, hi = np.percentile(arr, [2.5, 97.5])
    result.dG_ci95_lo_kcalmol = float(lo)
    result.dG_ci95_hi_kcalmol = float(hi)
    result.ok = True
    return result


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligand-smiles", required=True)
    ap.add_argument("--center", required=True,
                    help="comma-separated x,y,z in Å")
    ap.add_argument("--box", default="22,22,22")
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--base-seed", type=int, default=1)
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--num-modes", type=int, default=3)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--cpu-per-job", type=int, default=1)
    args = ap.parse_args()

    center = tuple(float(x) for x in args.center.split(","))
    box = tuple(float(x) for x in args.box.split(","))

    r = monte_carlo_dock(
        receptor_pdb=args.receptor,
        ligand_smiles=args.ligand_smiles,
        center_A=center, box_size_A=box,
        n_samples=args.n, base_seed=args.base_seed,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        workers=args.workers, cpu_per_job=args.cpu_per_job)
    print(r.summary())
    if r.ok:
        print(f"  samples: {r.dG_samples_kcalmol}")
    sys.exit(0 if r.ok else 1)
