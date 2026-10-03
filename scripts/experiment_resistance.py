#!/usr/bin/env python3
"""Does the schedule change which cells are left, not just how many?

Phase 2, within-line. Treatment does not sample a population at random:
it removes the sensitive cells and leaves the resistant tail behind, so
the survivors are a biased sample of what started. The engine should
show that without being told to, because each representative cell is a
lineage carrying its own anti-apoptotic reserve and drug accumulation,
and a lineage that survives keeps those values as it doubles.

Two things are measured over a two-week course:

1. **Selection.** How far the survivors' mean reserve and accumulation
   drift from the starting population's (both 1.0 by construction).
   Drift in the resistant direction means selection is happening, and
   which of the two drifts says what the drug selects *for*.

2. **Schedule.** The same total exposure given continuously versus in
   pulses with recovery gaps, which is the metronomic-versus-intermittent
   comparison. Population at day 14 and the drift in each arm.

What this engine can and cannot show is worth stating in advance. Its
heterogeneity is drawn once at the start, so it can select from a
pre-existing resistant tail but cannot generate new resistance by
mutation. It also has no carrying capacity — lineages grow
exponentially — so sensitive and resistant cells never compete for
space. Adaptive-therapy results, where keeping a sensitive population
alive suppresses a resistant one, depend on exactly that competition and
are therefore out of reach here; `cellsim/cell/agents.py` is the module
with a lattice and contact inhibition.

Usage:
    python scripts/experiment_resistance.py            # ~10 min
    python scripts/experiment_resistance.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "resistance_experiment.json"
DAYS = 14
# Concentration per drug chosen to kill most but not all of the population
# over the course, so there is something left to select among.
CONC = {"cisplatin": 15.0, "doxorubicin": 0.08, "paclitaxel": 0.04}


def _pulsed(conc: float, on_h: float, off_h: float, total_h: float):
    """Piecewise-constant schedule: `on_h` at `conc`, then `off_h` at zero."""
    segs, t = [], 0.0
    while t < total_h:
        segs.append((t, conc))
        segs.append((min(t + on_h, total_h), 0.0))
        t += on_h + off_h
    return lambda th: next(c for s, c in reversed(segs) if th >= s - 1e-12)


def _run(line, drug, dose_fn, hours, n_cells, seed, p):
    return simulate(line, drug, 0.0, t_end_h=hours, n_cells=n_cells, seed=seed, p=p,
                    k_cyc=calibrate_cycle_scale(line, p), dose_fn=dose_fn,
                    record_every_h=24.0)


def _summary(res):
    pop = res.final
    w = pop.weight * pop.alive
    if w.sum() <= 0:
        return {"alive_weight": 0.0, "mean_bcl2": None, "mean_uptake": None,
                "n_lineages": 0}
    return {"alive_weight": float(w.sum()),
            "mean_bcl2": float(np.average(pop.het_bcl2, weights=w)),
            "mean_uptake": float(np.average(pop.het_uptake, weights=w)),
            "n_lineages": int(pop.alive.sum())}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", default="A549")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    n_cells = 128 if a.quick else 256
    hours = DAYS * 24.0
    line = get_line(a.line)
    p = Params()
    report = {"line": line.name, "days": DAYS, "n_cells": n_cells, "conc_uM": CONC,
              "drugs": {}}

    print(f"{line.name}, {DAYS} days, {n_cells} lineages. Starting population has mean\n"
          f"reserve and accumulation 1.000 by construction; drift shows selection.\n")
    print(f"{'drug':<13}{'arm':<14}{'alive wt':>11}{'lineages':>10}"
          f"{'mean reserve':>14}{'mean uptake':>13}")
    print("-" * 75)
    for name, conc in CONC.items():
        drug = get_drug(name)
        arms = {
            "continuous": (lambda th, c=conc: c),
            "pulsed 1d/2d": _pulsed(conc * 3.0, 24.0, 48.0, hours),   # matched total exposure
        }
        report["drugs"][name] = {}
        for arm, fn in arms.items():
            res = _run(line, drug, fn, hours, n_cells, 1, p)
            s = _summary(res)
            report["drugs"][name][arm] = s
            mb = f"{s['mean_bcl2']:.3f}" if s["mean_bcl2"] is not None else "-"
            mu = f"{s['mean_uptake']:.3f}" if s["mean_uptake"] is not None else "-"
            print(f"{name if arm=='continuous' else '':<13}{arm:<14}"
                  f"{s['alive_weight']:>11.4g}{s['n_lineages']:>10}{mb:>14}{mu:>13}")
        print()

    print("reserve > 1 means survivors are enriched for high anti-apoptotic reserve;\n"
          "uptake < 1 means enriched for low drug accumulation. The pulsed arm gives\n"
          "the same total exposure at 3x the concentration for a third of the time.\n"
          "Note the engine draws its heterogeneity once, so this is selection from a\n"
          "pre-existing tail, not evolution of new resistance, and lineages do not\n"
          "compete for space — see this script's docstring.")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
