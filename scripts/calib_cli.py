#!/usr/bin/env python3
"""`cellsim calib <yaml>` — run one calibration + render scatter.

Given a calibration-set YAML (e.g.
benchmarks/dock/streptavidin_calibration.yaml), dock each entry,
report Pearson/Spearman/MAE/RMSE/conformal_q95, and save the
scatter plot as PNG.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("yaml", help="calibration-set YAML")
    ap.add_argument("--exh", type=int, default=32,
                    help="Vina exhaustiveness (default 32)")
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--cache", default=None,
                    help="SQLite cache path for Vina reuse")
    ap.add_argument("--save", default=None,
                    help="save scatter PNG to this path "
                         "(default: print stats only)")
    args = ap.parse_args(argv)

    from src.cache import Cache
    from src.uq import run_calibration, render_calibration_result

    cache = Cache(args.cache) if args.cache else None
    r = run_calibration(
        args.yaml, exhaustiveness=args.exh,
        num_modes=args.num_modes,
        seed=args.seed, cpu=args.cpu, cache=cache)
    print()
    print(r.summary())
    if args.save:
        render_calibration_result(
            r, save=args.save, show=False,
            title=f"CellSim calibration  —  {Path(args.yaml).stem}")
    return 0 if r.ok else 1


if __name__ == "__main__":
    sys.exit(main())
