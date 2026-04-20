"""Sobol global-sensitivity analysis over docking inputs.

Extends Layer 1.6 beyond seed-Monte-Carlo to answer the biologist-
UX question "which knob actually moves the ΔG prediction?". Uses
SALib's Saltelli sampling + Sobol analysis — non-AI, published
methodology (Saltelli 2010 Comput Phys Commun 181:259).

Default inputs swept:
  - exhaustiveness    (8 → 64, integer)
  - box_scale_factor  (0.8 → 1.2, multiplicative)
  - center_jitter_A   (0 → 2 Å, isotropic Gaussian shift)

For each Saltelli sample, Vina is run and the top-pose ΔG recorded.
Sobol.analyze then decomposes the ΔG variance into:
  - S1 (first-order):  fraction of variance explained by each input
                       alone
  - ST (total-effect): includes that input plus all its interactions

Biologist interpretation:
  - If S1[exhaustiveness] is large, users should max out exh for
    production runs.
  - If S1[center_jitter] is large, the docking is sensitive to
    search-box placement and users should be careful to centre on
    the crystal pose / fpocket output.

Compute budget: Saltelli sampling needs N × (2D + 2) runs for D
inputs. With D=3 and N=8, that's 64 Vina runs per compound. Feasible
for one-off analysis of a top hit, not for the CI PR gate.
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

logger = logging.getLogger(__name__)


@dataclass
class SobolResult:
    """Envelope for a Sobol sensitivity run."""

    receptor_pdb: str
    ligand_smiles: str
    ok: bool
    reason: str = ""

    input_names: list = field(default_factory=list)
    input_bounds: list = field(default_factory=list)
    n_base: int = 0
    n_total_runs: int = 0
    n_ok_runs: int = 0

    # Per-input: S1, S1_conf, ST, ST_conf (95 % CI half-width from
    # Saltelli analysis).
    s1: dict = field(default_factory=dict)
    s1_conf: dict = field(default_factory=dict)
    st: dict = field(default_factory=dict)
    st_conf: dict = field(default_factory=dict)

    dG_samples: list = field(default_factory=list)
    wall_seconds: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] Sobol  {self.reason}"
        lines = [
            f"[OK]   Sobol dock  N_base={self.n_base}  "
            f"runs={self.n_ok_runs}/{self.n_total_runs}  "
            f"({self.wall_seconds:.1f} s)"
        ]
        for name in self.input_names:
            s1 = self.s1.get(name, float("nan"))
            s1c = self.s1_conf.get(name, float("nan"))
            st = self.st.get(name, float("nan"))
            stc = self.st_conf.get(name, float("nan"))
            lines.append(
                f"  {name:<18s}  S1 = {s1:+.3f} ± {s1c:.3f}   "
                f"ST = {st:+.3f} ± {stc:.3f}")
        return "\n".join(lines)


def _run_sample(args: tuple) -> Optional[float]:
    """Worker: one Vina run with Saltelli-sampled parameters.

    Args packed into a tuple because multiprocessing workers need
    picklable inputs and the SALib numpy array gets unpacked into
    python-native types here.
    """
    (receptor, smiles, center_base, box_base, seed, timeout_s,
     exh, box_scale, center_jitter, jitter_dir) = args
    import numpy as np
    dx = center_jitter * jitter_dir[0]
    dy = center_jitter * jitter_dir[1]
    dz = center_jitter * jitter_dir[2]
    center = (center_base[0] + dx, center_base[1] + dy,
              center_base[2] + dz)
    box = tuple(b * box_scale for b in box_base)
    try:
        r = dock_ligand(
            receptor, smiles,
            center_A=center, box_size_A=box,
            exhaustiveness=int(round(exh)), num_modes=3,
            seed=seed, cpu=1, timeout_s=timeout_s)
    except Exception:
        return None
    if not r.ok or not r.poses:
        return None
    return float(r.poses[0].affinity_kcalmol)


def sobol_dock(
    *,
    receptor_pdb: str,
    ligand_smiles: str,
    center_A: tuple,
    box_size_A: tuple,
    n_base: int = 8,
    seed: int = 1,
    workers: int = 1,
    timeout_s: int = 600,
    include_inputs: tuple = ("exhaustiveness", "box_scale",
                              "center_jitter_A"),
) -> SobolResult:
    """Sobol sensitivity over docking inputs.

    `n_base` controls the Saltelli base-sample count; total runs are
    n_base * (2D + 2) where D = len(include_inputs). Results only
    become stable for n_base ≥ 32; below that the CI on Sobol
    indices is wide. The default 8 is a fast demo — honest signal
    for "is exhaustiveness dominant?" but not publication-quality.
    """
    import multiprocessing as mp
    import time

    import numpy as np

    try:
        from SALib.sample import saltelli
        from SALib.analyze import sobol as sobol_analyze
    except ImportError as e:
        return SobolResult(
            receptor_pdb=receptor_pdb, ligand_smiles=ligand_smiles,
            ok=False, reason=f"SALib not installed: {e}")

    names = list(include_inputs)
    bounds_map = {
        "exhaustiveness":  [8.0, 64.0],
        "box_scale":       [0.8, 1.2],
        "center_jitter_A": [0.0, 2.0],
    }
    unknown = [n for n in names if n not in bounds_map]
    if unknown:
        return SobolResult(
            receptor_pdb=receptor_pdb, ligand_smiles=ligand_smiles,
            ok=False, reason=f"unknown input(s): {unknown}")

    problem = {
        "num_vars": len(names),
        "names": names,
        "bounds": [bounds_map[n] for n in names],
    }
    X = saltelli.sample(problem, n_base,
                        calc_second_order=False)
    n_total = X.shape[0]
    logger.info("Sobol will run %d Vina evaluations", n_total)

    # Fixed isotropic jitter direction per run; biologist expects
    # reproducibility so tie it to the run index.
    rng = np.random.default_rng(seed)
    jitter_dirs = rng.normal(size=(n_total, 3))
    jitter_dirs /= np.linalg.norm(jitter_dirs, axis=1, keepdims=True) + 1e-12

    def _unpack(i: int) -> tuple:
        row = X[i]
        values = dict(zip(names, row))
        exh = values.get("exhaustiveness", 32.0)
        box_scale = values.get("box_scale", 1.0)
        center_jitter = values.get("center_jitter_A", 0.0)
        return (receptor_pdb, ligand_smiles,
                tuple(center_A), tuple(box_size_A),
                seed + i, timeout_s,
                exh, box_scale, center_jitter,
                tuple(jitter_dirs[i]))

    tasks = [_unpack(i) for i in range(n_total)]

    t0 = time.time()
    if workers <= 1:
        raw = [_run_sample(t) for t in tasks]
    else:
        with mp.Pool(workers) as pool:
            raw = list(pool.imap(_run_sample, tasks))

    # SALib needs a fully-populated numeric Y; impute missing with
    # the mean of observed values (the conservative default).
    observed = [v for v in raw if v is not None]
    if len(observed) < 0.5 * n_total:
        return SobolResult(
            receptor_pdb=receptor_pdb, ligand_smiles=ligand_smiles,
            ok=False, reason=f"only {len(observed)}/{n_total} "
                              "Vina runs succeeded; aborting analysis",
            n_base=n_base, n_total_runs=n_total,
            n_ok_runs=len(observed),
            dG_samples=[v for v in raw if v is not None],
            wall_seconds=time.time() - t0)
    mean_dG = float(np.mean(observed))
    Y = np.array([v if v is not None else mean_dG for v in raw])

    try:
        Si = sobol_analyze.analyze(
            problem, Y, calc_second_order=False,
            print_to_console=False)
    except Exception as e:
        return SobolResult(
            receptor_pdb=receptor_pdb, ligand_smiles=ligand_smiles,
            ok=False, reason=f"SALib analyze failed: {e}",
            n_base=n_base, n_total_runs=n_total,
            n_ok_runs=len(observed),
            dG_samples=observed, wall_seconds=time.time() - t0)

    return SobolResult(
        receptor_pdb=receptor_pdb, ligand_smiles=ligand_smiles,
        ok=True,
        input_names=names,
        input_bounds=[bounds_map[n] for n in names],
        n_base=n_base, n_total_runs=n_total,
        n_ok_runs=len(observed),
        s1={n: float(v) for n, v in zip(names, Si["S1"])},
        s1_conf={n: float(v) for n, v in zip(names, Si["S1_conf"])},
        st={n: float(v) for n, v in zip(names, Si["ST"])},
        st_conf={n: float(v) for n, v in zip(names, Si["ST_conf"])},
        dG_samples=observed,
        wall_seconds=time.time() - t0,
    )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligand-smiles", required=True)
    ap.add_argument("--center", required=True,
                    help="comma-separated x,y,z Å")
    ap.add_argument("--box", default="22,22,22")
    ap.add_argument("--n-base", type=int, default=8,
                    help="Saltelli base-sample count (total runs = "
                         "N * (2D+2); default 8 → 64 for D=3)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    center = tuple(float(x) for x in args.center.split(","))
    box = tuple(float(x) for x in args.box.split(","))

    r = sobol_dock(
        receptor_pdb=args.receptor,
        ligand_smiles=args.ligand_smiles,
        center_A=center, box_size_A=box,
        n_base=args.n_base, seed=args.seed,
        workers=args.workers)
    print(r.summary())
    sys.exit(0 if r.ok else 1)
