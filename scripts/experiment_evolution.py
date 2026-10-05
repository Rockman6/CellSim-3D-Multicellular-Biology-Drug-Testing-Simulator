#!/usr/bin/env python3
"""Does resistance EVOLVE, or is it only selected from what was already there?

Phase 2 could not tell these apart. Its heterogeneity was drawn once, so
every surviving lineage had been resistant from the start; the engine
could select, never originate. The dish can do both, and this experiment
separates them by removing the tail:

* **Homogeneous start** (`Params.het_sigma = 0`): every cell begins
  identical, so selection has nothing to act on. Any resistance that
  appears must have arisen during treatment.
* **With the usual spread** (0.3): both routes are available.

Each is run with and without heritable change at division
(`DishParams.mutation_rate`), under a drug-free control, continuous
exposure, and pulsed exposure at matched total exposure.

What to look for, and what would falsify the model:

1. Homogeneous + no mutation must give NO resistance at all. If it does,
   something other than inheritance is drifting, and the measurement is
   wrong.
2. Homogeneous + mutation should give resistance that grows with time.
   That is evolution, not selection.
3. The fold-resistance reached should be in the range derived lines
   actually reach: two- to eight-fold for clinically realistic pulsed
   selection (McDermott et al. 2014 Front Oncol 4:40). Far above that
   means the assumed mutation rate or effect size is too generous.

Resistance is read out two ways: the surviving population's mean drug
accumulation and apoptotic reserve (1.0 = the starting cells), and the
IC50 of the survivors re-measured as a fresh dose-response, which is what
a lab would actually report.

Usage:
    python scripts/experiment_evolution.py            # ~25 min
    python scripts/experiment_evolution.py --quick    # ~6 min
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import DishParams, run_dish  # noqa: E402
from cellsim.cell.engine import (Params, calibrate_cycle_scale, ic50,  # noqa: E402
                                 ic50_from_curve, simulate)
from cellsim.cell.library import get_drug, get_line  # noqa: E402

OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "evolution_experiment.json"
DRUG = "cisplatin"
MUTATION_RATE = 0.05          # per daughter per division, any heritable route
MUTATION_SD = 0.3             # log-sd of the effect, as big as the population spread


def pulsed(conc: float, on_h: float = 24.0, off_h: float = 48.0):
    period = on_h + off_h
    return lambda t_h: conc if (t_h % period) < on_h - 1e-9 else 0.0


def survivor_ic50(dish, drug, p: Params, k_cyc: float, *, guess_uM: float,
                  per_conc: int = 48) -> float:
    """Re-measure the survivors as a lab would: take the traits they
    inherited and run a fresh 72 h dose-response on them.

    Built from `simulate` rather than `dose_response` because every
    concentration must see the SAME trait sample — otherwise the curve
    mixes a dose effect with a sampling difference."""
    n = len(dish.pos)
    if n < 8:
        return float("nan")
    idx = np.resize(np.arange(n), per_conc)          # one sample of survivors
    bcl2 = dish.cells.het_bcl2[idx]
    uptake = dish.cells.het_uptake[idx]
    conc = np.concatenate([[0.0], guess_uM * np.logspace(-1.5, 1.5, 13)])
    groups = np.repeat(np.arange(len(conc)), per_conc)
    traits = {"bcl2": np.tile(bcl2, len(conc)), "uptake": np.tile(uptake, len(conc))}
    res = simulate(dish.line, drug, conc[groups], t_end_h=72.0, n_cells=len(groups),
                   groups=groups, p=p, k_cyc=k_cyc, traits=traits, record_every_h=72.0)
    return ic50_from_curve(conc[1:], res.viability()[1:])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", default="A549")
    ap.add_argument("--days", type=float, default=14.0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)
    days = 7.0 if a.quick else a.days
    n_seed = 200 if a.quick else 400
    grid = 64 if a.quick else 80

    line, drug = get_line(a.line), get_drug(DRUG)
    base = Params()
    k = calibrate_cycle_scale(line, base)
    ic_base = ic50(line, drug, guess_uM=10.0, k_cyc=k, p=base)
    print(f"{line.name} + {DRUG}: baseline 72 h IC50 {ic_base:.3g} uM; {days:g} days, "
          f"{n_seed} cells on a {grid}x{grid} lattice\n")

    report = {"line": line.name, "drug": DRUG, "days": days, "n_seed": n_seed,
              "grid": grid, "baseline_ic50_uM": ic_base, "mutation_rate": MUTATION_RATE,
              "mutation_effect_sd": MUTATION_SD, "arms": {}}
    schedules = {"untreated": (lambda t: 0.0),
                 "continuous 1x IC50": (lambda t: ic_base),
                 "pulsed 3x, 1d/2d": pulsed(3.0 * ic_base)}
    print(f"{'start':<14}{'mutation':<10}{'schedule':<20}{'live':>7}{'uptake':>9}"
          f"{'reserve':>9}{'muts':>7}{'IC50 uM':>10}{'fold':>7}")
    print("-" * 93)
    for het_label, het in (("homogeneous", 0.0), ("spread 0.3", 0.3)):
        p = Params(het_sigma=het, cycle_cv=0.25)
        for mut_label, rate in (("off", 0.0), ("on", MUTATION_RATE)):
            ic_control = float("nan")
            for sched_label, fn in schedules.items():
                dp = DishParams(geometry="monolayer", mutation_rate=rate,
                                mutation_effect_sd=MUTATION_SD if rate else 0.0)
                res = run_dish(line, drug, fn, t_end_h=days * 24.0, n_seed=n_seed,
                               grid_sites=grid, record_every_h=24.0, seed=1, p=p, dp=dp)
                rec, dish = res.record, res.dish
                ic = survivor_ic50(dish, drug, p, k, guess_uM=ic_base)
                # Fold-resistance against the UNTREATED arm of the same
                # condition, measured the same way, which is the parental
                # line a lab keeps in parallel. Comparing against a global
                # baseline would mix in the heterogeneity and cycle-spread
                # settings of the arm rather than the treatment.
                if sched_label == "untreated":
                    ic_control = ic
                fold = ic / ic_control if np.isfinite(ic) and np.isfinite(ic_control) else float("nan")
                key = f"{het_label} | mutation {mut_label} | {sched_label}"
                report["arms"][key] = {
                    "n_live": rec.n_live, "mean_uptake": rec.mean_uptake,
                    "mean_reserve": rec.mean_reserve, "n_mutations": rec.n_mutations,
                    "survivor_ic50_uM": ic, "control_ic50_uM": ic_control,
                    "fold_resistance": fold, "eradicated": rec.n_live[-1] == 0}
                print(f"{het_label:<14}{mut_label:<10}{sched_label:<20}{rec.n_live[-1]:>7}"
                      f"{rec.mean_uptake[-1]:>9.3f}{rec.mean_reserve[-1]:>9.3f}"
                      f"{rec.n_mutations[-1]:>7}"
                      f"{(f'{ic:.3g}' if np.isfinite(ic) else '-'):>10}"
                      f"{(f'{fold:.2f}x' if np.isfinite(fold) else '-'):>7}")
        print()
    print("Fold-resistance is against the untreated arm of the same condition, which is\n"
          "the parental line a lab keeps in parallel. The control that matters is a\n"
          "homogeneous start with mutation off: with no pre-existing spread and no\n"
          "heritable change, no resistance can appear, and an arm that is eradicated\n"
          "reports no IC50 at all rather than a number from a handful of cells.")
    if not a.no_write:
        OUT_JSON.write_text(json.dumps(report, indent=1, default=float) + "\n")
        print(f"\nwrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
