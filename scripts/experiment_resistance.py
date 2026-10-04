#!/usr/bin/env python3
"""Does the schedule change which cells are left, not just how many?

Phase 2, within-line. Treatment does not sample a population at random:
it removes the sensitive cells and leaves a resistant tail, so the
survivors are a biased sample of what started. In the engine each
representative cell is a lineage carrying its own anti-apoptotic
reserve (Bcl-2 multiplier) and drug accumulation (uptake multiplier),
and a lineage that survives keeps both as it doubles.

How it is measured, and why
---------------------------
The first version of this experiment drew 128 random lineages, treated
them, and averaged the survivors. Under the harsher schedule that left
ONE surviving lineage for cisplatin and two for doxorubicin, so its
headline (that high-intensity dosing selects on uptake rather than
reserve) was a statement about one or two random draws. Re-measured
exactly, it was backwards. See docs/VALIDATION.md, "Corrections".

This version maps the SURVIVAL LANDSCAPE instead (cellsim.cell.selection): a grid of chosen
(reserve, uptake) values, several cells per grid point at random cycle
phases, each run through the schedule. Survival as a function of trait
is the selection itself; integrating it against the population's trait
distribution (log-normal, Params.het_sigma) gives the expected survivor
means exactly, however rare survival is. A logistic fit of survival on
log reserve and log uptake gives the steepness of the boundary along
each axis; their ratio |b_uptake / b_reserve| says which trait the
schedule selects on.

Arms, at matched total exposure over two weeks:
* continuous, at the drug's own 72 h IC50;
* pulsed, 3 x IC50 for 24 h then 48 h drug-free, repeated.

Limits, stated in advance. Heterogeneity is drawn once, so this is
selection from a pre-existing tail, not evolution of new resistance;
lineages grow without a carrying capacity, so they never compete.
Mutation and competition belong to the spatial dish.

Usage:
    python scripts/experiment_resistance.py            # ~3 min
    python scripts/experiment_resistance.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale, ic50  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.selection import (RESERVE_GRID, UPTAKE_GRID, pulsed,  # noqa: E402
                                    summarize, survival_map)

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "resistance_experiment.json"
DAYS = 14
IC50_GUESS_UM = {"cisplatin": 10.0, "doxorubicin": 0.05, "paclitaxel": 0.03}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", default="A549")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    reps = 3 if a.quick else 6
    hours = DAYS * 24.0
    line = get_line(a.line)
    p = Params()
    k = calibrate_cycle_scale(line, p)
    report = {"line": line.name, "days": DAYS, "reps_per_grid_point": reps,
              "reserve_grid": RESERVE_GRID.tolist(), "uptake_grid": UPTAKE_GRID.tolist(),
              "het_sigma": p.het_sigma, "drugs": {}}

    print(f"{line.name}, {DAYS} days. Survivor means are exact expectations over a\n"
          f"log-normal population (sigma {p.het_sigma}); 1.000 = no selection.\n")
    print(f"{'drug':<13}{'arm':<12}{'survival':>10}{'reserve':>10}{'uptake':>9}"
          f"{'|b_u/b_r|':>11}")
    print("-" * 65)
    for name, guess in IC50_GUESS_UM.items():
        drug = get_drug(name)
        ic = ic50(line, drug, guess_uM=guess, k_cyc=k, p=p)
        arms = {"continuous": lambda t_h, c=ic: c, "pulsed 1d/2d": pulsed(3.0 * ic)}
        report["drugs"][name] = {"ic50_uM": ic}
        for arm, fn in arms.items():
            S = survival_map(line, drug, fn, hours=hours, reps=reps, k_cyc=k, p=p)
            sel = summarize(S, sigma=p.het_sigma)
            report["drugs"][name][arm] = {"survival_map": S.tolist(), **sel}
            mr = sel["survivor_mean_reserve"]
            mu = sel["survivor_mean_uptake"]
            print(f"{name if arm == 'continuous' else '':<13}{arm:<12}"
                  f"{sel['population_survival']:>10.3g}"
                  f"{(f'{mr:.3f}' if mr is not None else '-'):>10}"
                  f"{(f'{mu:.3f}' if mu is not None else '-'):>9}"
                  f"{sel['uptake_over_reserve']:>11.2f}")
        print()
    print("|b_u/b_r| > 1: the survival boundary is set mainly by drug accumulation;\n"
          "< 1: mainly by the apoptotic reserve. Continuous low-dose exposure selects\n"
          "on accumulation; intermittent high-dose exposure shifts selection toward\n"
          "the reserve, which buffers a transient insult but not a sustained one.")

    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
