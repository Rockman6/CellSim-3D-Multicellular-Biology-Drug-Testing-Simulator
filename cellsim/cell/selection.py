"""cellsim.cell.selection — which cells does a schedule leave behind?

Treatment does not thin a population at random; it removes sensitive
cells and leaves a resistant tail. In the engine each representative
cell carries two heterogeneity multipliers (mean 1 across the
population): its anti-apoptotic reserve (Bcl-2) and its drug
accumulation (uptake). This module measures selection on them exactly.

Sampling random lineages and averaging the survivors is the obvious
method and the wrong one: under a harsh schedule only one or two
lineages survive, and their mean is a statement about one or two random
draws. (That is how an earlier version of this analysis got its
headline backwards; see docs/VALIDATION.md, "Corrections".)

Instead, `survival_map` runs a GRID of chosen trait values, several
cells per point at random cycle phases, and records the fraction alive
at the end. That map is the selection function itself. `summarize`
integrates it against the population's log-normal trait distribution,
which gives the survivors' expected means however rare survival is, and
fits a logistic boundary whose slopes say which trait the schedule
selects on.
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from cellsim.cell.engine import Params, simulate
from cellsim.cell.library import CellLine, Drug

RESERVE_GRID = np.geomspace(0.5, 3.0, 13)
UPTAKE_GRID = np.geomspace(0.15, 1.5, 13)


def pulsed(conc_uM: float, on_h: float = 24.0, off_h: float = 48.0) -> Callable[[float], float]:
    """`on_h` hours at `conc_uM`, then `off_h` drug-free, repeated."""
    period = on_h + off_h
    return lambda t_h: conc_uM if (t_h % period) < on_h - 1e-9 else 0.0


def survival_map(line: CellLine, drug: Optional[Drug], dose_fn: Callable[[float], float], *,
                 hours: float, reps: int = 6, k_cyc: Optional[float] = None,
                 p: Params = Params(), reserve: np.ndarray = RESERVE_GRID,
                 uptake: np.ndarray = UPTAKE_GRID, seed: int = 7) -> np.ndarray:
    """Fraction of `reps` cells alive after `hours`, for every
    (reserve, uptake) grid point. Shape (len(reserve), len(uptake))."""
    rr, uu = np.meshgrid(reserve, uptake, indexing="ij")
    traits = {"bcl2": np.repeat(rr.ravel(), reps), "uptake": np.repeat(uu.ravel(), reps)}
    res = simulate(line, drug, 0.0, t_end_h=hours, n_cells=len(traits["bcl2"]), k_cyc=k_cyc,
                   p=p, dose_fn=dose_fn, record_every_h=hours, seed=seed, traits=traits)
    return res.final.alive.reshape(len(reserve), len(uptake), reps).mean(-1)


def summarize(S: np.ndarray, reserve: np.ndarray = RESERVE_GRID,
              uptake: np.ndarray = UPTAKE_GRID, sigma: float = Params().het_sigma) -> dict:
    """Expected survival and survivor means under independent log-normal
    traits (mean 1, log-sd `sigma`), and the logistic boundary slopes.

    `uptake_over_reserve` = |b_uptake / b_reserve|: above 1 the survival
    boundary is set mainly by drug accumulation, below 1 mainly by the
    apoptotic reserve. When the map is perfectly separable the slopes
    themselves grow large and only their ratio is meaningful; the
    survivor means are the more robust readout."""
    lr, lu = np.meshgrid(np.log(reserve), np.log(uptake), indexing="ij")
    mu = -0.5 * sigma ** 2
    w = np.exp(-((lr - mu) ** 2 + (lu - mu) ** 2) / (2 * sigma ** 2))
    w /= w.sum()
    surv = float((w * S).sum())
    out = {"population_survival": surv,
           "survivor_mean_reserve": float((w * S * np.exp(lr)).sum() / surv) if surv > 0 else None,
           "survivor_mean_uptake": float((w * S * np.exp(lu)).sum() / surv) if surv > 0 else None,
           # the same expectations with no selection, i.e. the grid's own
           # rendering of a mean-1 population (close to 1 on these grids)
           "baseline_mean_reserve": float((w * np.exp(lr)).sum()),
           "baseline_mean_uptake": float((w * np.exp(lu)).sum())}
    X = np.column_stack([np.ones(lr.size), lr.ravel(), lu.ravel()])
    y = S.ravel()
    beta = np.zeros(3)
    for _ in range(100):                      # Newton-Raphson, lightly ridged
        q = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30, 30)))
        W = q * (1 - q) + 1e-9
        beta += np.linalg.solve(X.T @ (W[:, None] * X) + 1e-4 * np.eye(3), X.T @ (y - q))
    out["boundary_slope_reserve"] = float(beta[1])
    out["boundary_slope_uptake"] = float(beta[2])
    out["uptake_over_reserve"] = (float(abs(beta[2] / beta[1])) if abs(beta[1]) > 1e-6
                                  else float("inf"))
    return out
