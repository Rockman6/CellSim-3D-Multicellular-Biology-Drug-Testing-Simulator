"""cellsim.calibrate — tune the engine to YOUR cells, and carry the error.

`cellsim.plate` turns a plate-reader export into a dose-response curve
with a measured confidence interval. This module turns that curve into
the two things the engine needs to describe the line it came from:

    doubling_time_from_controls()  the line's own growth, from the
                                   untreated wells of the same plate
    calibrate_potency()            the drug's single potency constant,
                                   fitted so the simulated IC50 matches
                                   the measured one
    predict_band()                 any later prediction run at the low,
                                   central and high ends of that
                                   constant, so the answer is a band

The point of the last one. A calibrated constant is not a number, it is
a number with an interval, and an interval that is not propagated is
decoration. Every prediction here is run three times and reported as a
range, because a schedule comparison whose two arms differ by less than
that range is not a result.

What this does NOT do is fit several constants at once. Each drug in
`cellsim.cell.library` has exactly one fitted quantity, named by
`Drug.fit_target`, and that discipline is what keeps the engine from
absorbing the data instead of being tested by it. A curve constrains one
number; asking it for four would give four numbers that jointly fit and
individually mean nothing.
"""
from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass
from typing import Optional, Sequence, Union

import numpy as np

from cellsim.cell.engine import Params, calibrate_cycle_scale, dose_response, ic50
from cellsim.cell.library import CellLine, Drug, get_drug, get_line

# Half kill: the simulated IC50 equals a concentration exactly when the
# viability there is one half.
HALF = 0.5

__all__ = ["doubling_time_from_controls", "calibrate_potency", "predict_band",
           "DoublingTime", "Calibration"]


@dataclass
class DoublingTime:
    hours: float
    ci: tuple[float, float]
    doublings_observed: float
    note: str = ""

    def summary(self) -> str:
        lo, hi = self.ci
        return (f"doubling time {self.hours:.1f} h (95 % CI {lo:.1f}-{hi:.1f}), "
                f"from {self.doublings_observed:.2f} doublings"
                + (f"  [{self.note}]" if self.note else ""))


def doubling_time_from_controls(assay_hours: float, t0: Sequence[float],
                                final: Sequence[float]) -> DoublingTime:
    """The line's doubling time from the untreated wells of its own plate.

    `t0` and `final` are the untreated readings at the start and end of
    the assay (raw signal is fine if it is proportional to cell number).
    The interval comes from the well-to-well spread of both, propagated
    through ``T_d = t * ln2 / ln(final/t0)``.

    A short assay is the thing to watch: if the control barely doubled,
    the log ratio is small, its relative error is large, and the doubling
    time is poorly determined however many wells were read. That is
    reported rather than hidden.
    """
    t0 = np.asarray(t0, float)
    final = np.asarray(final, float)
    if t0.size == 0 or final.size == 0:
        raise ValueError("need both time-zero and end-point control wells")
    if t0.mean() <= 0 or final.mean() <= 0:
        raise ValueError("control readings must be positive")
    fold = float(final.mean() / t0.mean())
    if fold <= 1.0:
        raise ValueError(f"the control did not grow (fold change {fold:.3f}); "
                         f"a doubling time cannot be read from it")
    doublings = math.log2(fold)
    td = assay_hours / doublings
    # Propagate the standard errors of both means through log(final/t0).
    def rel_se(a):
        return (float(a.std(ddof=1)) / math.sqrt(len(a)) / float(a.mean())
                if len(a) > 1 else 0.0)
    se_log = math.sqrt(rel_se(t0) ** 2 + rel_se(final) ** 2)
    if se_log > 0:
        lo_f, hi_f = fold * math.exp(-1.96 * se_log), fold * math.exp(1.96 * se_log)
        ci = (assay_hours / math.log2(max(hi_f, 1.0 + 1e-9)),
              assay_hours / math.log2(max(lo_f, 1.0 + 1e-9)))
    else:
        ci = (td, td)
    note = ""
    if doublings < 1.0:
        note = (f"the control managed only {doublings:.2f} doublings; a longer "
                f"assay would pin this down far better")
    return DoublingTime(hours=float(td), ci=(float(ci[0]), float(ci[1])),
                        doublings_observed=float(doublings), note=note)


@dataclass
class Calibration:
    """A drug's fitted constant for one line, with its interval.

    A KNOWN OFFSET, measured rather than suspected. The fitted constant
    comes out about 10 % below the truth, consistently rather than on
    average: 0.88x, 0.91x and 0.92x for lines perturbed 0.5x, 1x and 2x,
    with the measurement noise removed entirely.

    The cause is that two definitions of "the IC50" are in play and they
    differ by roughly that much. This module targets the concentration at
    which median viability is one half. `cellsim.cell.engine.ic50`
    instead interpolates a 0.5 crossing on a log grid, and on a curve as
    steep as this engine's it lands about 1.1x low (calibrating to its
    answer and measuring back gives 1.07-1.15x). Targeting the estimator
    instead would cost a full concentration series per iteration — the
    thing that makes calibration fast enough to use — to chase a bias
    that sits well inside the +-20 % run-to-run spread the engine has
    anyway (see `predict_band`). So it is recorded here rather than
    removed, and `predict_band`'s interval is what should be read.
    """
    line: str
    drug: str
    fit_target: str
    value: float
    ci: tuple[float, float]
    measured_ic50_uM: float
    measured_ic50_ci: tuple[float, float]
    simulated_ic50_uM: float
    note: str = ""

    def drug_at(self, which: str = "value") -> Drug:
        """The library drug with this constant substituted in.

        `which` is 'value', 'low' or 'high' — the central estimate and the
        two ends of the interval, which is what `predict_band` runs."""
        v = {"value": self.value, "low": self.ci[0], "high": self.ci[1]}[which]
        return dataclasses.replace(get_drug(self.drug), **{self.fit_target: v})

    def summary(self) -> str:
        lo, hi = self.ci
        mlo, mhi = self.measured_ic50_ci
        gives = (f"; the calibrated engine gives {self.simulated_ic50_uM:.4g} uM"
                 if np.isfinite(self.simulated_ic50_uM) else "")
        return (f"{self.drug} on {self.line}: {self.fit_target} = {self.value:.4g} "
                f"(95 % CI {lo:.3g}-{hi:.3g}), matching a measured IC50 of "
                f"{self.measured_ic50_uM:.4g} uM ({mlo:.3g}-{mhi:.3g}){gives}"
                + (f"\n  note: {self.note}" if self.note else ""))


def _fit_constant(line: CellLine, drug: Drug, target_uM: float, *, n_cells: int,
                  iters: int, p: Params, k_cyc: float, seeds: Sequence[int] = (1,),
                  bracket_decades: float = 1.5, max_decades: float = 3.0,
                  report_ic50: bool = True) -> tuple[float, float]:
    """Bisection in log10 of the drug's one fitted constant, until the
    simulated IC50 matches `target_uM`. Returns (constant, achieved IC50).

    The direction is read off the bracket ends rather than assumed, so a
    potency gain (bigger = more potent) and a partition coefficient
    (bigger = more drug inside) are both handled.

    The search starts narrow and widens only if the target is not
    bracketed, up to `max_decades`. Both bounds matter:

    * Starting narrow is a speed decision. Fixing the bracket wide wasted
      most of the iterations closing a range no real cell line occupies.
    * Stopping at three decades is a correctness decision. Allowed to run
      further, the bisection will happily "reach" any measured IC50 by
      driving the constant somewhere physically absurd — a partition
      coefficient of six million, a cell concentrating drug a
      million-fold. A constant a thousand times from the reference means
      the measurement is not describing this mechanism, and the honest
      answer is to say so rather than to return the number that fits.
    """
    field = drug.fit_target
    if not field:
        raise ValueError(f"{drug.name} has no fit_target; there is no single "
                         f"constant to calibrate")
    start = getattr(drug, field)
    if start <= 0:
        raise ValueError(f"{drug.name}.{field} is {start}; expected a positive constant")

    # Bisect on VIABILITY AT THE TARGET CONCENTRATION, not on the IC50.
    # The two are the same condition — the simulated IC50 equals T exactly
    # when viability at T is 0.5 — but computing an IC50 means running a
    # whole concentration series (about thirty simulations) where this
    # needs one. That is a ~15x saving per iteration, and it is what makes
    # a calibration cheap enough to sit in the test suite.
    #
    # Returned viability is inverted around 0.5 so that the search still
    # reads "too weak / too strong" the same way for both kinds of
    # constant, and the caller's bracketing logic is unchanged.
    cache: dict[float, float] = {}

    def sim(log_v: float) -> float:
        """A monotone stand-in for the simulated IC50: viability at the
        target concentration, mapped so that larger means a higher IC50.

        `seeds` averages the criterion over repeats. It defaults to one,
        because averaging was tried and MEASURED NOT TO HELP: three seeds
        returned byte-identical constants for every case tried, at three
        times the cost, and moved the systematic offset by under half a
        per cent (0.869/0.913/0.946 to 0.878/0.910/0.943). The offset is
        definitional rather than noise — see the note on `Calibration` —
        so averaging the noise cannot touch it. The option is kept
        because a future criterion might need it; the default is not.

        What DOES vary with the seed is the answer itself: the same
        target calibrated under different seeds gives constants spanning
        about 1.14x, which is the engine's own stochasticity showing
        through. That is what `predict_band` spans, and why it must."""
        key = round(log_v, 6)
        if key not in cache:
            d = dataclasses.replace(drug, **{field: 10.0 ** key})
            runs = []
            for sd in seeds:
                viab, _ = dose_response(line, d, np.array([target_uM]),
                                        n_cells_per_conc=n_cells, p=p, k_cyc=k_cyc,
                                        seed=int(sd))
                runs.append(float(viab[1]))
            cache[key] = float(np.median(runs))
        return cache[key]

    centre = math.log10(start)
    width = min(bracket_decades, max_decades)
    while True:
        lo, hi = centre - width, centre + width
        ic_lo, ic_hi = sim(lo), sim(hi)
        falling = ic_lo >= ic_hi                # viability at T falls as potency rises
        weak, strong = (ic_lo, ic_hi) if falling else (ic_hi, ic_lo)
        if strong <= HALF <= weak:
            break                                      # the target is bracketed
        if width >= max_decades:
            # Out of reach inside a physically sensible range. Report the
            # closest the mechanism gets, as a viability, for the caller
            # to turn into a message.
            return (math.nan, strong) if strong > HALF else (math.nan, weak)
        width = min(width * 2.0, max_decades)

    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if (sim(mid) > HALF) == falling:
            lo = mid
        else:
            hi = mid
    final = 0.5 * (lo + hi)
    if not report_ic50:
        return 10.0 ** final, math.nan
    # One real IC50 at the end, on the value actually returned, so the
    # caller can see what the calibrated engine gives rather than what
    # the search was aiming at. It is a whole concentration series, which
    # costs as much as the search itself, so a caller that only wants the
    # constant can turn it off.
    fitted = dataclasses.replace(drug, **{field: 10.0 ** final})
    achieved = ic50(line, fitted, guess_uM=target_uM, n_cells_per_conc=n_cells,
                    p=p, k_cyc=k_cyc)
    return 10.0 ** final, achieved


def calibrate_potency(line: Union[str, CellLine], drug: Union[str, Drug],
                      measured_ic50_uM: float,
                      measured_ic50_ci: Optional[tuple[float, float]] = None, *,
                      params: Params = Params(), n_cells: int = 32,
                      iters: int = 12, seeds: Sequence[int] = (1,),
                      report_ic50: bool = True) -> Calibration:
    """Fit the drug's one constant so the engine reproduces a measured IC50.

    `measured_ic50_ci` should be the interval from `cellsim.plate.fit_4pl`.
    The constant's interval is obtained by repeating the fit at each end
    of it, which is the honest propagation: a measurement that pins the
    IC50 to within a factor of two cannot pin the constant any better.
    Omit it and the result carries a point estimate and says so.
    """
    ln = line if isinstance(line, CellLine) else get_line(line)
    dg = drug if isinstance(drug, Drug) else get_drug(drug)
    k = calibrate_cycle_scale(ln, params)
    value, achieved = _fit_constant(ln, dg, measured_ic50_uM, n_cells=n_cells,
                                    iters=iters, p=params, k_cyc=k, seeds=seeds,
                                    report_ic50=report_ic50)
    if not np.isfinite(value):
        raise ValueError(
            f"no value of {dg.name}.{dg.fit_target} reproduces an IC50 of "
            f"{measured_ic50_uM:.4g} uM on {ln.name}: the closest the engine gets "
            f"is {achieved:.4g} uM. Either the measurement is far outside what this "
            f"mechanism can produce, or the line's other inputs (doubling time, "
            f"p53 status) are wrong for these cells.")
    note = ""
    if measured_ic50_ci is None:
        ci = (value, value)
        note = ("no interval was given for the measured IC50, so this constant has "
                "none either; predictions from it are point estimates")
    else:
        lo_ic, hi_ic = measured_ic50_ci
        ends = []
        for target in (lo_ic, hi_ic):
            # The ends only pin the interval; their own IC50s are not
            # reported, so the extra concentration series is skipped.
            v, _ = _fit_constant(ln, dg, target, n_cells=n_cells, iters=iters,
                                 p=params, k_cyc=k, seeds=seeds, report_ic50=False)
            ends.append(v)
        if not all(np.isfinite(ends)):
            ci = (value, value)
            note = ("one end of the measured interval is outside what this mechanism "
                    "can reach, so the constant's interval is not reported")
        else:
            ci = (float(min(ends)), float(max(ends)))
    return Calibration(line=ln.name, drug=dg.name, fit_target=dg.fit_target,
                       value=float(value), ci=ci,
                       measured_ic50_uM=float(measured_ic50_uM),
                       # When no interval was supplied, the reported one is the
                       # point estimate twice — of the IC50, not of the constant.
                       measured_ic50_ci=(tuple(float(x) for x in measured_ic50_ci)
                                         if measured_ic50_ci
                                         else (float(measured_ic50_uM),) * 2),
                       simulated_ic50_uM=float(achieved), note=note)


def predict_band(calibration: Calibration, fn, *, seeds: Sequence[int] = (1, 2, 3),
                 labels=("low", "value", "high")):
    """Run `fn(drug, seed)` across the calibrated interval AND across
    seeds, and return the envelope.

    `fn` takes a `Drug` and a seed and returns a number — a surviving
    fraction, a fold change, whatever the question is.

    **Both sources of uncertainty belong in the band**, and leaving one
    out was a real error here. The first version varied only the
    constant, at a single seed, and its bands missed the truth in a
    hold-out test more often than they contained it
    (`scripts/validate_calibration.py`). The reason is that the engine is
    stochastic: repeating the same simulation with a different seed moves
    its IC50 by about 1.2x, which is as large as the calibration
    uncertainty it was reporting. A band that counts the measurement's
    error and ignores the model's own is not an honest interval, however
    carefully the first part was computed.

    Two schedules whose bands overlap have not been shown to differ.
    """
    per_label: dict = {}
    for which, label in zip(("low", "value", "high"), labels):
        drug = calibration.drug_at(which)
        vals = [fn(drug, int(s)) for s in seeds]
        finite = [v for v in vals if isinstance(v, (int, float)) and np.isfinite(v)]
        per_label[label] = float(np.median(finite)) if finite else float("nan")
        per_label.setdefault("_all", []).extend(finite)
    allv = per_label.pop("_all")
    out = dict(per_label)
    if allv:
        out["range"] = (float(min(allv)), float(max(allv)))
        out["spread"] = float(max(allv) - min(allv))
        out["n_runs"] = len(allv)
    return out
