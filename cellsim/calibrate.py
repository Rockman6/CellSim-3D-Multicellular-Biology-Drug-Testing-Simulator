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

from cellsim.cell.engine import Params, calibrate_cycle_scale, ic50
from cellsim.cell.library import CellLine, Drug, get_drug, get_line

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
    """A drug's fitted constant for one line, with its interval."""
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
        return (f"{self.drug} on {self.line}: {self.fit_target} = {self.value:.4g} "
                f"(95 % CI {lo:.3g}-{hi:.3g}), matching a measured IC50 of "
                f"{self.measured_ic50_uM:.4g} uM ({mlo:.3g}-{mhi:.3g}); the "
                f"calibrated engine gives {self.simulated_ic50_uM:.4g} uM"
                + (f"\n  note: {self.note}" if self.note else ""))


def _fit_constant(line: CellLine, drug: Drug, target_uM: float, *, n_cells: int,
                  iters: int, p: Params, k_cyc: float) -> tuple[float, float]:
    """Bisection in log10 of the drug's one fitted constant, until the
    simulated IC50 matches `target_uM`. Returns (constant, achieved IC50).

    The direction is read off the bracket ends rather than assumed, so a
    potency gain (bigger = more potent) and a partition coefficient
    (bigger = more drug inside) are both handled.
    """
    field = drug.fit_target
    if not field:
        raise ValueError(f"{drug.name} has no fit_target; there is no single "
                         f"constant to calibrate")
    start = getattr(drug, field)
    if start <= 0:
        raise ValueError(f"{drug.name}.{field} is {start}; expected a positive constant")
    lo, hi = math.log10(start) - 4, math.log10(start) + 4

    def sim(log_v: float) -> float:
        d = dataclasses.replace(drug, **{field: 10.0 ** log_v})
        return ic50(line, d, guess_uM=target_uM, n_cells_per_conc=n_cells, p=p, k_cyc=k_cyc)

    ic_lo, ic_hi = sim(lo), sim(hi)
    falling = ic_lo >= ic_hi
    weak, strong = (ic_lo, ic_hi) if falling else (ic_hi, ic_lo)
    if strong > target_uM:
        return math.nan, strong          # cannot be made potent enough
    if weak < target_uM:
        return math.nan, weak            # cannot be made weak enough
    ic = math.nan
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ic = sim(mid)
        if (ic > target_uM) == falling:
            lo = mid
        else:
            hi = mid
    return 10.0 ** (0.5 * (lo + hi)), ic


def calibrate_potency(line: Union[str, CellLine], drug: Union[str, Drug],
                      measured_ic50_uM: float,
                      measured_ic50_ci: Optional[tuple[float, float]] = None, *,
                      params: Params = Params(), n_cells: int = 32,
                      iters: int = 12) -> Calibration:
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
                                    iters=iters, p=params, k_cyc=k)
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
            v, _ = _fit_constant(ln, dg, target, n_cells=n_cells, iters=iters,
                                 p=params, k_cyc=k)
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


def predict_band(calibration: Calibration, fn, *, labels=("low", "value", "high")):
    """Run `fn(drug)` at each end of the calibrated interval and in the
    middle, and return {label: result} plus the spread.

    `fn` takes a `Drug` and returns a number — a surviving fraction, a
    fold change, whatever the question is. The spread it comes back with
    is the part that matters: two schedules that differ by less than the
    calibration's own width have not been shown to differ.
    """
    out = {}
    for which, label in zip(("low", "value", "high"), labels):
        out[label] = fn(calibration.drug_at(which))
    vals = [v for v in out.values() if isinstance(v, (int, float)) and np.isfinite(v)]
    if vals:
        out["spread"] = float(max(vals) - min(vals))
        out["range"] = (float(min(vals)), float(max(vals)))
    return out
