"""cellsim.api — the simulator from a notebook, in tidy rows.

Biologists work in notebooks and R, not in a CLI, and they want a table
they can plot or export. Every function here returns one: a list of
dicts, or a `pandas.DataFrame` if pandas is installed (it is not a
dependency). Columns are tidy — one row per observation — so the result
goes straight into seaborn, ggplot via `to_csv`, or a spreadsheet.

The functions map to the questions the engine has actually been checked
against (`docs/VALIDATION.md`), and no further:

    knockout()     the same curve with a gene knocked out, knocked down
                   or overexpressed — what the pathway does without it
    starve()       a spheroid at a given medium glucose, which sets the
                   viable rim as much as oxygen does
    killing()      a cytotoxicity assay: effector cells on a monolayer at
                   a range of E:T ratios -> specific lysis
    curve()        a dose-response curve and its IC50
    ic50_spread()  the same IC50 under several seeds, and how far it
                   moves — run this before quoting one
    gr_curve()     the same run as GR values, which is the axis a
                   measured plate can be compared on
    exposure()     the iso-effect curve: the concentration needed for
                   half kill at each exposure time, which is how a
                   schedule decision is actually made
    washout()      how a colony recovers after the drug is removed
    combination()  two drugs, each arm against independent action
    spheroid()     a 3-D aggregate: size, hypoxia, necrotic core

What this module deliberately does NOT offer is a function that ranks
cell lines by sensitivity. That was measured to be unreachable from the
markers available, and is out of scope (`docs/VALIDATION.md`).

Example
-------
>>> from cellsim.api import curve, exposure
>>> df = curve("A549", "paclitaxel")              # doctest: +SKIP
>>> exposure("A549", "paclitaxel", hours=[3, 24, 72])   # doctest: +SKIP
"""
from __future__ import annotations

from typing import Iterable, Optional, Sequence, Union

import numpy as np

from cellsim.cell.engine import (Params, calibrate_cycle_scale, dose_response, ic50,
                                 ic50_from_curve, simulate)
from cellsim.cell.library import CELL_LINES, DRUGS, CellLine, Drug, get_drug, get_line

__all__ = ["curve", "exposure", "washout", "combination", "spheroid", "gr_curve",
           "ic50_spread", "knockout", "starve", "killing", "genes", "lines", "drugs"]


def _table(rows: list[dict]):
    """A DataFrame when pandas is available, otherwise the rows as given."""
    try:
        import pandas as pd
    except ImportError:
        return rows
    return pd.DataFrame(rows)


def _line(x: Union[str, CellLine]) -> CellLine:
    return x if isinstance(x, CellLine) else get_line(x)


def _drug(x: Union[str, Drug]) -> Drug:
    return x if isinstance(x, Drug) else get_drug(x)


def lines():
    """Every cell line, with the inputs that distinguish it."""
    return _table([{"line": c.name, "tissue": c.tissue,
                    "doubling_time_h": c.doubling_time_h,
                    "p53_functional": c.p53_functional,
                    "efflux_level": c.efflux_level,
                    "bcl2_level": c.bcl2_level,
                    "source": c.source} for c in CELL_LINES.values()])


def drugs():
    """Every drug, its mechanism, and which constant was fitted."""
    return _table([{"drug": d.name, "mechanism": d.mechanism,
                    "pgp_substrate": d.pgp_substrate,
                    "fitted_constant": d.fit_target or "",
                    "fitted_value": getattr(d, d.fit_target) if d.fit_target else None,
                    "source": d.source} for d in DRUGS.values()])


def curve(line: Union[str, CellLine], drug: Union[str, Drug], *,
          concentrations: Optional[Sequence[float]] = None, hours: float = 72.0,
          exposure_h: Optional[float] = None, n_cells: int = 48, seed: int = 1,
          params: Params = Params()):
    """Viability against concentration, one row per concentration.

    `exposure_h` washes the drug out at that time and reads the plate at
    `hours`, which is how an exposure-duration experiment is run.
    Concentrations default to four decades around the drug's own IC50.
    The IC50 of the curve is repeated on every row for convenience."""
    ln, dg = _line(line), _drug(drug)
    k = calibrate_cycle_scale(ln, params)
    if concentrations is None:
        centre = ic50(ln, dg, guess_uM=_guess(dg), t_end_h=hours, n_cells_per_conc=24,
                      k_cyc=k, p=params, exposure_h=exposure_h)
        if not np.isfinite(centre):
            centre = _guess(dg) * 100
        concentrations = centre * np.logspace(-2, 2, 13)
    conc = np.asarray(list(concentrations), float)
    viab, res = dose_response(ln, dg, conc, t_end_h=hours, n_cells_per_conc=n_cells,
                              seed=seed, p=params, k_cyc=k, exposure_h=exposure_h)
    ic = ic50_from_curve(conc, viab[1:])
    growth = res.alive_weight[-1, 0] / max(res.alive_weight[0, 0], 1e-12)
    return _table([{"line": ln.name, "drug": dg.name, "conc_uM": float(c),
                    "viability": float(v), "ic50_uM": float(ic),
                    "exposure_h": exposure_h if exposure_h is not None else hours,
                    "readout_h": hours, "control_fold_change": float(growth)}
                   for c, v in zip(conc, viab[1:])])


def gr_curve(line: Union[str, CellLine], drug: Union[str, Drug], *,
             concentrations: Optional[Sequence[float]] = None, hours: float = 72.0,
             exposure_h: Optional[float] = None, n_cells: int = 48, seed: int = 1,
             params: Params = Params()):
    """The same simulation as `curve()`, reported as GR values.

    A measured plate analysed by `cellsim.plate.gr_metrics` and a
    simulation analysed here are then the same quantity, so they can be
    plotted on one axis and compared. Viability cannot be compared that
    way between a fast and a slow line, which is the whole reason GR
    exists (Hafner et al. 2016 Nat Methods 13:521).

    GR = 1 is untreated growth, 0 is complete cytostasis, below 0 is net
    cell loss. The engine's alive weight is the cell count GR needs, and
    the untreated arm of the same run supplies the control."""
    from cellsim.plate import gr_metrics
    ln, dg = _line(line), _drug(drug)
    k = calibrate_cycle_scale(ln, params)
    if concentrations is None:
        centre = ic50(ln, dg, guess_uM=_guess(dg), t_end_h=hours, n_cells_per_conc=24,
                      k_cyc=k, p=params, exposure_h=exposure_h)
        if not np.isfinite(centre):
            centre = _guess(dg) * 100
        concentrations = centre * np.logspace(-2, 2, 9)
    conc = np.asarray(list(concentrations), float)
    _, res = dose_response(ln, dg, conc, t_end_h=hours, n_cells_per_conc=n_cells,
                           seed=seed, p=params, k_cyc=k, exposure_h=exposure_h)
    # Group 0 is the untreated control that dose_response prepends.
    t0 = float(res.alive_weight[0, 0])
    control = float(res.alive_weight[-1, 0])
    treated = [float(res.alive_weight[-1, i + 1]) for i in range(len(conc))]
    gr = gr_metrics(conc, treated, control=control, t0=t0)
    return _table([{"line": ln.name, "drug": dg.name, "conc_uM": float(c),
                    "gr": float(g), "gr50_uM": gr.gr50_uM, "gr_max": gr.gr_max,
                    "gr_aoc": gr.gr_aoc,
                    "control_doublings": gr.doublings_in_assay,
                    "exposure_h": exposure_h if exposure_h is not None else hours,
                    "readout_h": hours}
                   for c, g in zip(gr.conc_uM, gr.gr)])


def ic50_spread(line: Union[str, CellLine], drug: Union[str, Drug], *,
                seeds: Sequence[int] = (1, 2, 3, 4, 5), hours: float = 72.0,
                exposure_h: Optional[float] = None, n_cells: int = 32,
                params: Params = Params()):
    """The same IC50 measured under several seeds, and how far it moves.

    Worth running before quoting an IC50 from this engine. The number is
    a random variable: repeating an identical simulation under a
    different seed moves it by 10-20 %, because the engine draws a finite
    sample of lineages and a population's IC50 is set by its resistant
    tail rather than by its mean. More cells barely helps — with
    heterogeneity switched off the spread collapses to under 1 %, which
    is how the cause was identified (docs/VALIDATION.md).

    Returns one row per seed plus the median and the spread, so a caller
    can see whether a difference they care about is larger than the
    engine's own. Two numbers from this engine that differ by less than
    the spread do not differ.
    """
    ln, dg = _line(line), _drug(drug)
    k = calibrate_cycle_scale(ln, params)
    vals = [ic50(ln, dg, guess_uM=_guess(dg), t_end_h=hours, n_cells_per_conc=n_cells,
                 k_cyc=k, p=params, seed=int(s), exposure_h=exposure_h) for s in seeds]
    finite = [v for v in vals if np.isfinite(v)]
    median = float(np.median(finite)) if finite else float("nan")
    spread = float(max(finite) / min(finite)) if len(finite) > 1 else float("nan")
    return _table([{"line": ln.name, "drug": dg.name, "seed": int(s),
                    "ic50_uM": float(v), "median_ic50_uM": median,
                    "spread_fold": spread, "n_finite": len(finite),
                    "exposure_h": exposure_h if exposure_h is not None else hours,
                    "readout_h": hours}
                   for s, v in zip(seeds, vals)])


def exposure(line: Union[str, CellLine], drug: Union[str, Drug],
             hours: Iterable[float] = (1, 3, 6, 12, 24, 48, 72), *,
             readout_h: float = 72.0, n_cells: int = 32, seed: int = 1,
             params: Params = Params()):
    """The iso-effect curve: the concentration needed to halve the colony
    when the drug is washed out at each exposure time.

    This is the measurement a schedule decision rests on. If killing were
    governed by total exposure (AUC), `conc_x_hours` would be constant.
    Where no concentration reaches half kill, `c50_uM` is infinite — which
    is the answer for a drug that needs cells to transit a cycle phase."""
    ln, dg = _line(line), _drug(drug)
    k = calibrate_cycle_scale(ln, params)
    base = ic50(ln, dg, guess_uM=_guess(dg), t_end_h=readout_h, n_cells_per_conc=n_cells,
                k_cyc=k, p=params)
    rows = []
    for T in hours:
        T = float(T)
        guess = base * readout_h / T if 0 < T < readout_h else base
        c50 = ic50(ln, dg, guess_uM=guess, t_end_h=readout_h, n_cells_per_conc=n_cells,
                   k_cyc=k, p=params, seed=seed,
                   exposure_h=None if T >= readout_h else T)
        rows.append({"line": ln.name, "drug": dg.name, "exposure_h": T,
                     "readout_h": readout_h, "c50_uM": float(c50),
                     "conc_x_hours_uM_h": float(c50 * T),
                     "reaches_half_kill": bool(np.isfinite(c50))})
    return _table(rows)


def washout(line: Union[str, CellLine], drug: Union[str, Drug], conc_uM: float,
            exposure_h: float, *, hours: float = 168.0, n_cells: int = 96,
            record_every_h: float = 6.0, seed: int = 1, params: Params = Params()):
    """Colony size over time for a pulse of `conc_uM` for `exposure_h`,
    beside an untreated control: how far it falls and whether it regrows.

    `fold_change` is relative to the start; `relative_to_control` is the
    treated colony divided by the untreated one at the same time, which is
    what a plate reader reports."""
    ln, dg = _line(line), _drug(drug)
    k = calibrate_cycle_scale(ln, params)
    dose = lambda t: (conc_uM if t < exposure_h - 1e-9 else 0.0)  # noqa: E731
    treated = simulate(ln, dg, 0.0, t_end_h=hours, n_cells=n_cells, k_cyc=k, p=params,
                       dose_fn=dose, record_every_h=record_every_h, seed=seed)
    control = simulate(ln, None, 0.0, t_end_h=hours, n_cells=n_cells, k_cyc=k, p=params,
                       record_every_h=record_every_h, seed=seed)
    t0, c0 = treated.alive_weight[0, 0], control.alive_weight[0, 0]
    rows = []
    for i, t in enumerate(treated.t_h):
        a, c = treated.alive_weight[i, 0], control.alive_weight[i, 0]
        # The integrator accumulates its clock, so a tick lands on
        # 11.999999999999833 rather than 12. Round to the nearest
        # thousandth of an hour: a table a human groups or plots by time
        # should not carry that noise.
        rows.append({"line": ln.name, "drug": dg.name, "conc_uM": conc_uM,
                     "exposure_h": exposure_h, "t_h": round(float(t), 3),
                     "drug_present": bool(t < exposure_h - 1e-9),
                     "fold_change": float(a / max(t0, 1e-12)),
                     "control_fold_change": float(c / max(c0, 1e-12)),
                     "relative_to_control": float(a / max(c, 1e-12))})
    return _table(rows)


def combination(line: Union[str, CellLine], drug_a: Union[str, Drug],
                drug_b: Union[str, Drug], conc_a_uM: float, conc_b_uM: float, *,
                each_h: float = 36.0, readout_h: float = 120.0, n_cells: int = 128,
                seeds: Sequence[int] = (1, 2, 3), params: Params = Params()):
    """Each drug alone, both together, and each order, against the
    independent-action expectation (the product of the two alone).

    `ratio_to_independent` above 1 is antagonism, below 1 synergy. Every
    arm gives each drug for the same number of hours, so only the timing
    differs. Averaged over `seeds` with a standard deviation, because one
    run carries a few per cent of noise on these ratios — enough to invent
    a sequence effect that is not there."""
    ln = _line(line)
    a, b = _drug(drug_a), _drug(drug_b)
    k = calibrate_cycle_scale(ln, params)
    pair = [a, b]
    half = each_h
    arms = {
        f"{a.name} alone": ([a], lambda t: (conc_a_uM if t < half else 0.0,)),
        f"{b.name} alone": ([b], lambda t: (conc_b_uM if t < half else 0.0,)),
        "simultaneous": (pair, lambda t: (conc_a_uM, conc_b_uM) if t < half else (0.0, 0.0)),
        f"{a.name} then {b.name}": (pair, lambda t: (conc_a_uM, 0.0) if t < half else
                                    ((0.0, conc_b_uM) if t < 2 * half else (0.0, 0.0))),
        f"{b.name} then {a.name}": (pair, lambda t: (0.0, conc_b_uM) if t < half else
                                    ((conc_a_uM, 0.0) if t < 2 * half else (0.0, 0.0))),
    }
    surv: dict[str, list[float]] = {name: [] for name in arms}
    for seed in seeds:
        ctrl = simulate(ln, None, 0.0, t_end_h=readout_h, n_cells=n_cells, k_cyc=k,
                        p=params, record_every_h=readout_h, seed=seed).alive_weight[-1, 0]
        for name, (ds, fn) in arms.items():
            v = simulate(ln, ds, 0.0, t_end_h=readout_h, n_cells=n_cells, k_cyc=k, p=params,
                         dose_fn=fn, record_every_h=readout_h, seed=seed).alive_weight[-1, 0]
            surv[name].append(float(v / max(ctrl, 1e-12)))
    solo = np.array(surv[f"{a.name} alone"]) * np.array(surv[f"{b.name} alone"])
    rows = []
    for name, vals in surv.items():
        v = np.array(vals)
        alone = name.endswith("alone")
        ratio = None if alone else v / np.maximum(solo, 1e-12)
        rows.append({"line": ln.name, "arm": name, "readout_h": readout_h,
                     "surviving_fraction": float(v.mean()),
                     "surviving_sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                     "independent_expectation": float(solo.mean()),
                     "ratio_to_independent": None if alone else float(ratio.mean()),
                     "ratio_sd": None if alone else (float(ratio.std(ddof=1))
                                                     if len(v) > 1 else 0.0),
                     "n_seeds": len(vals)})
    return _table(rows)


def genes():
    """Which genes can be perturbed, and what each one is in the engine."""
    from cellsim.cell.perturb import GENES
    return _table([{"gene": g, "engine_node": v["node"], "kind": v["kind"],
                    "what_it_is": v["note"]} for g, v in sorted(GENES.items())])


def knockout(line: Union[str, CellLine], drug: Union[str, Drug], gene: str, *,
             mode: str = "knockout", x: float = 1.0,
             concentrations: Optional[Sequence[float]] = None, hours: float = 72.0,
             n_cells: int = 32, seed: int = 1, params: Params = Params()):
    """A dose-response with one gene perturbed, beside the unperturbed one.

    Returns both arms in one table so the comparison is the thing you get,
    not something you have to assemble: a `viability` and a
    `viability_wild_type` for every concentration, and the fold change
    between them.

    `cellsim.api.genes()` lists what can be perturbed. A gene the engine
    does not represent is refused rather than silently doing nothing."""
    from cellsim.cell.perturb import Perturbation
    ln, dg = _line(line), _drug(drug)
    pert = [Perturbation(gene, mode, x)]
    k = calibrate_cycle_scale(ln, params)
    if concentrations is None:
        centre = ic50(ln, dg, guess_uM=_guess(dg), n_cells_per_conc=24, k_cyc=k, p=params)
        if not np.isfinite(centre):
            centre = _guess(dg) * 100
        concentrations = centre * np.logspace(-1.5, 1.5, 7)
    conc = np.asarray(list(concentrations), float)
    wt, _ = dose_response(ln, dg, conc, t_end_h=hours, n_cells_per_conc=n_cells,
                          seed=seed, p=params, k_cyc=k)
    ko, _ = dose_response(ln, dg, conc, t_end_h=hours, n_cells_per_conc=n_cells,
                          seed=seed, p=params, k_cyc=k, perturbations=pert)
    return _table([{"line": ln.name, "drug": dg.name, "perturbation": str(pert[0]),
                    "conc_uM": float(c), "viability": float(v),
                    "viability_wild_type": float(w),
                    "fold_vs_wild_type": float(v / w) if w > 0 else float("nan")}
                   for c, v, w in zip(conc, ko[1:], wt[1:])])


def starve(line: Union[str, CellLine], *, glucose_mM: float = 5.5, days: float = 6.0,
           n_seed: int = 2000, grid_sites: int = 48, record_every_h: float = 24.0,
           seed: int = 1, params: Params = Params()):
    """Grow a spheroid in medium carrying a given glucose concentration.

    DMEM is 25 mM and RPMI 11 mM, both far above anything a cell needs;
    physiological blood is about 5.5 mM and a tumour's interior much less.
    Below roughly 0.3 mM the interior cells stop cycling and below 0.06 mM
    they die, so lowering this thins the viable rim the way lowering
    oxygen does (Freyer & Sutherland 1986).

    Returns one row per time point, with the glucose at the centre beside
    the usual size and necrosis columns."""
    from cellsim.cell.dish import DishParams, run_dish
    ln = _line(line)
    dp = DishParams(geometry="spheroid", glucose_medium_mM=glucose_mM)
    res = run_dish(ln, None, None, t_end_h=days * 24.0, n_seed=n_seed,
                   grid_sites=grid_sites, record_every_h=record_every_h, seed=seed,
                   p=params, dp=dp)
    r = res.record
    centre = (float(res.dish.profile_glucose[0])
              if len(res.dish.profile_glucose) else float("nan"))
    return _table([{"line": ln.name, "glucose_medium_mM": glucose_mM,
                    "t_h": round(float(t), 3), "n_live": int(r.n_live[i]),
                    "n_necrotic": int(r.n_necrotic[i]),
                    "radius_um": float(r.radius_um[i]),
                    "necrotic_radius_um": float(r.necrotic_radius_um[i]),
                    "quiescent_frac": float(r.quiescent_frac[i]),
                    "glucose_centre_mM_final": centre}
                   for i, t in enumerate(r.t_h)])


def killing(line: Union[str, CellLine], *, et_ratios: Sequence[float] = (0, 0.1, 0.3, 1, 3),
            hours: float = 24.0, n_tumour: int = 400, grid: int = 40,
            seeds: Sequence[int] = (1, 2), params: Params = Params(), dish_params=None):
    """A cytotoxicity assay: effector cells on a tumour monolayer.

    Effectors crawl on their own layer, engage a tumour cell on contact,
    kill it after about an hour, then rest before engaging again and stop
    after a few kills (Halle et al. 2016). `et_ratios` are effector:target
    ratios; the 0 arm is the untreated control that specific lysis is
    measured against.

    Returns one row per ratio, averaged over `seeds`."""
    from cellsim.cell.dish import Dish, DishParams
    ln = _line(line)
    rows, base = [], None
    for et in et_ratios:
        lives, kills = [], []
        for sd in seeds:
            dp = dish_params or DishParams(geometry="monolayer", spacing_um=20.0,
                                           field_dt_h=0.5, dt_h=0.02)
            d = Dish(ln, n_seed=n_tumour, grid_sites=(grid, grid), seed=int(sd),
                     p=params, dp=dp)
            d.add_effectors(int(round(et * n_tumour)))
            for _ in range(int(round(hours / dp.field_dt_h))):
                d.advance(None)
            lives.append(int(d.cells.alive.sum()))
            kills.append(int(d.n_killed_by_effectors))
        live = float(np.mean(lives))
        if base is None:
            base = live
        rows.append({"line": ln.name, "et_ratio": float(et), "hours": hours,
                     "n_live": live, "n_killed_by_effectors": float(np.mean(kills)),
                     "specific_lysis": float(1.0 - live / base) if base else float("nan")})
    return _table(rows)


def spheroid(line: Union[str, CellLine], drug: Optional[Union[str, Drug]] = None, *,
             conc_uM: float = 0.0, start_h: float = 0.0, stop_h: Optional[float] = None,
             days: float = 6.0, n_seed: int = 2000, grid_sites: int = 48,
             record_every_h: float = 12.0, seed: int = 1, params: Params = Params(),
             dish_params=None):
    """Grow a 3-D aggregate, optionally dosed between `start_h` and
    `stop_h`, and report its size, hypoxic fraction and necrotic core over
    time — one row per time point.

    Oxygen and drug are solved radially with cellular uptake, so the
    interior sees less of both than the rim."""
    from cellsim.cell.dish import DishParams, run_dish
    ln = _line(line)
    dg = _drug(drug) if drug is not None else None
    dp = dish_params or DishParams(geometry="spheroid")
    stop = stop_h if stop_h is not None else days * 24.0
    dose = None if dg is None else (lambda t: conc_uM if start_h <= t < stop else 0.0)
    res = run_dish(ln, dg, dose, t_end_h=days * 24.0, n_seed=n_seed, grid_sites=grid_sites,
                   record_every_h=record_every_h, seed=seed, p=params, dp=dp)
    r = res.record
    rows = []
    for i, t in enumerate(r.t_h):
        rows.append({"line": ln.name, "drug": None if dg is None else dg.name,
                     "conc_uM": conc_uM, "t_h": round(float(t), 3),
                     "day": round(float(t) / 24.0, 4),
                     "n_live": r.n_live[i], "n_apoptotic": r.n_apoptotic[i],
                     "n_necrotic": r.n_necrotic[i],
                     "radius_um": float(r.radius_um[i]),
                     "necrotic_radius_um": float(r.necrotic_radius_um[i]),
                     "hypoxic_fraction": float(r.hypoxic_frac[i]),
                     "quiescent_fraction": float(r.quiescent_frac[i]),
                     "min_o2_mmHg": float(r.min_o2_mmHg[i])})
    return _table(rows)


def _guess(drug: Drug) -> float:
    """A starting concentration for the IC50 search, by mechanism."""
    return {"tubulin": 0.03, "topo2": 0.05, "s_phase": 0.1, "mdm2": 1.0}.get(
        drug.mechanism, 10.0)
