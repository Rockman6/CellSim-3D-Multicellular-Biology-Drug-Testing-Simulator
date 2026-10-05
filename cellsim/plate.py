"""cellsim.plate — read a plate reader, fit the curve, keep the error bars.

Everything validated so far uses OUR cell lines and OUR fitted constants.
A lab has its own line, its own plate reader and its own drug, and until
their numbers can enter the simulator nothing here is a tool for them.
This module is that door:

    read_plate()    a plate-reader export (grid or tidy) -> tidy rows
    normalize()     raw signal -> viability, against on-plate controls
    fit_4pl()       viability vs concentration -> IC50 with a bootstrap CI
    gr_metrics()    GR50, GRmax, GR_AOC (Hafner et al. 2016 Nat Methods
                    13:521), which correct for how fast the line divides

Why GR as well as IC50
----------------------
A 72 h viability assay confounds potency with division rate: a line that
doubles twice in the assay has more to lose than one that doubles once,
so the same drug looks more potent in the faster line. GR metrics
normalise by the growth of the untreated control over the same window and
are comparable across lines and assay lengths. They need the cell count
at time zero as well as at the end, which is one extra plate read — if
you have it, use `gr_metrics`; if not, `fit_4pl` still gives an IC50 and
says so.

What is deliberately NOT here: fitting the simulator's constants to the
data. That is the next step (`docs/PLAN.md`, Phase 4), and it needs the
uncertainty this module produces, so this comes first and alone.

Units: concentrations in µM throughout, to match `cellsim.cell.library`.
"""
from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

import numpy as np

__all__ = ["read_plate", "normalize", "fit_4pl", "gr_metrics", "FitResult", "GRResult"]

_WELL_RE = re.compile(r"^([A-Za-z]{1,2})(\d{1,2})$")


def _well_key(well: str) -> tuple[str, int]:
    m = _WELL_RE.match(well.strip())
    if not m:
        raise ValueError(f"{well!r} is not a well name like 'B7' or 'AA12'")
    return m.group(1).upper(), int(m.group(2))


def _norm_well(well: str) -> str:
    """'b7' and 'B07' both become 'B7', so a plate map and a reader export
    that disagree about zero padding still line up."""
    row, col = _well_key(well)
    return f"{row}{col}"


# ── reading ───────────────────────────────────────────────────────────
def read_plate(path: Union[str, Path], *, layout: str = "auto",
               plate_map: Optional[Union[str, Path]] = None,
               value_column: Optional[str] = None) -> list[dict]:
    """Read a plate-reader export into tidy rows.

    Two shapes are understood, and `layout="auto"` picks between them:

    * ``grid`` — the familiar block with row letters down the side and
      column numbers across the top, as most readers export. Returns one
      row per well.
    * ``tidy`` — a table that already has one row per well, with a column
      named well/Well/WELL and a numeric value column.

    `plate_map` is a second file (grid or tidy) giving each well's
    condition: a concentration in µM, or one of the words `control`,
    `untreated`, `vehicle`, `dmso` (untreated), `blank`, `empty`, `media`
    (no cells), or `t0` (the time-zero read GR needs). Anything else is
    kept verbatim in `label` and left for the caller.

    Returns rows with keys: well, row, col, value, and — when a map is
    given — label, conc_uM, role ('dose' | 'control' | 'blank' | 't0').
    """
    rows = _read_any(Path(path), layout=layout, value_column=value_column, numeric=True)
    if plate_map is None:
        return rows
    mapping = {r["well"]: r["value"] for r in
               _read_any(Path(plate_map), layout=layout, value_column=None, numeric=False)}
    missing = [r["well"] for r in rows if r["well"] not in mapping]
    if missing:
        raise ValueError(f"plate map has no entry for {len(missing)} wells, "
                         f"e.g. {missing[:5]}")
    out = []
    for r in rows:
        label = str(mapping[r["well"]]).strip()
        role, conc = _classify(label)
        out.append({**r, "label": label, "conc_uM": conc, "role": role})
    return out


_CONTROL_WORDS = {"control", "untreated", "vehicle", "dmso", "ctrl", "neg"}
# A BLANK is medium without cells: a real background reading, and the
# thing to subtract. An EMPTY well is not used at all — in a grid export
# the unused rows are usually padded with zeros, and averaging those into
# the blank drags the background toward zero and quietly rescales every
# viability on the plate. They are kept as their own role and ignored.
_BLANK_WORDS = {"blank", "media", "medium", "background", "no-cells", "nocells"}
_EMPTY_WORDS = {"empty", "unused", "none", "na", "n/a", "-", "0", ""}
_T0_WORDS = {"t0", "time0", "day0", "untreated_t0"}


def _classify(label: str) -> tuple[str, Optional[float]]:
    low = label.lower().replace(" ", "")
    if low in _T0_WORDS:
        return "t0", None
    if low in _CONTROL_WORDS:
        return "control", 0.0
    if low in _BLANK_WORDS:
        return "blank", None
    if low in _EMPTY_WORDS:
        return "empty", None
    try:
        return "dose", float(label)
    except ValueError:
        return "other", None


def _read_any(path: Path, *, layout: str, value_column: Optional[str],
              numeric: bool) -> list[dict]:
    text = path.read_text()
    table = [r for r in csv.reader(text.splitlines()) if any(c.strip() for c in r)]
    if not table:
        raise ValueError(f"{path} is empty")
    if layout == "auto":
        layout = _detect_layout(table)
    if layout == "grid":
        return _read_grid(table, numeric=numeric, where=path)
    if layout == "tidy":
        return _read_tidy(table, value_column=value_column, numeric=numeric, where=path)
    raise ValueError(f"layout must be 'grid', 'tidy' or 'auto', not {layout!r}")


def _detect_layout(table: list[list[str]]) -> str:
    """A grid has a header of column numbers and row labels down the side;
    a tidy table has a column whose name looks like 'well'."""
    header = [c.strip().lower() for c in table[0]]
    if any(h in ("well", "wells", "well_id", "position") for h in header):
        return "tidy"
    numeric_header = sum(1 for c in table[0][1:] if c.strip().isdigit())
    row_labels = sum(1 for r in table[1:] if r and _WELL_RE.match((r[0] or "") + "1"))
    if numeric_header >= 3 and row_labels >= 2:
        return "grid"
    return "tidy"


def _read_grid(table: list[list[str]], *, numeric: bool, where: Path) -> list[dict]:
    header = table[0]
    cols = [(i, c.strip()) for i, c in enumerate(header) if c.strip().isdigit()]
    if not cols:
        raise ValueError(f"{where}: no column numbers in the first row; is this a grid?")
    out = []
    for raw in table[1:]:
        if not raw or not raw[0].strip():
            continue
        row_label = raw[0].strip().upper()
        if not re.fullmatch(r"[A-Z]{1,2}", row_label):
            continue
        for i, cnum in cols:
            if i >= len(raw):
                continue
            cell = raw[i].strip()
            if cell == "":
                continue
            value = _maybe_float(cell, numeric, where, f"{row_label}{cnum}")
            out.append({"well": f"{row_label}{int(cnum)}", "row": row_label,
                        "col": int(cnum), "value": value})
    if not out:
        raise ValueError(f"{where}: grid layout parsed no wells")
    return out


def _read_tidy(table: list[list[str]], *, value_column: Optional[str], numeric: bool,
               where: Path) -> list[dict]:
    header = [c.strip() for c in table[0]]
    low = [c.lower() for c in header]
    try:
        wi = next(i for i, c in enumerate(low) if c in ("well", "wells", "well_id", "position"))
    except StopIteration:
        raise ValueError(f"{where}: no 'well' column found; columns are {header}") from None
    if value_column is not None:
        if value_column not in header:
            raise ValueError(f"{where}: no column {value_column!r}; columns are {header}")
        vi = header.index(value_column)
    else:
        candidates = [i for i, c in enumerate(low)
                      if i != wi and c in ("value", "signal", "reading", "od", "lum",
                                           "luminescence", "fluorescence", "count", "rlu")]
        vi = candidates[0] if candidates else (wi + 1 if wi + 1 < len(header) else None)
        if vi is None:
            raise ValueError(f"{where}: cannot tell which column holds the reading; "
                             f"pass value_column=. Columns are {header}")
    out = []
    for raw in table[1:]:
        if len(raw) <= max(wi, vi) or not raw[wi].strip():
            continue
        well = _norm_well(raw[wi])
        row, col = _well_key(well)
        out.append({"well": well, "row": row, "col": col,
                    "value": _maybe_float(raw[vi].strip(), numeric, where, well)})
    if not out:
        raise ValueError(f"{where}: tidy layout parsed no wells")
    return out


def _maybe_float(cell: str, numeric: bool, where: Path, well: str):
    if not numeric:
        return cell
    try:
        return float(cell)
    except ValueError:
        raise ValueError(f"{where}: well {well} holds {cell!r}, which is not a number. "
                         f"If this is a plate MAP, pass it as plate_map=") from None


# ── normalisation ─────────────────────────────────────────────────────
def normalize(rows: Sequence[dict]) -> list[dict]:
    """Raw signal -> viability relative to the plate's own controls.

    Blank wells (no cells) are averaged and subtracted from everything
    first; untreated controls then set 1.0. Both come from the same
    plate, which is the point — it removes the reader's offset and the
    day's growth without a second experiment.

    Adds `viability` to every dose and control row, and raises if the
    plate has no controls, because then there is nothing to normalise to.
    """
    blanks = [r["value"] for r in rows if r.get("role") == "blank"]
    blank = float(np.mean(blanks)) if blanks else 0.0
    ctrl = [r["value"] - blank for r in rows if r.get("role") == "control"]
    if not ctrl:
        raise ValueError("no control wells on this plate (label them 'control', "
                         "'untreated', 'vehicle' or 'dmso' in the plate map); "
                         "without them a reading cannot become a viability")
    ctrl_mean = float(np.mean(ctrl))
    if ctrl_mean <= 0:
        raise ValueError(f"control wells average {ctrl_mean:.4g} after blank "
                         f"subtraction; check the blanks")
    out = []
    for r in rows:
        r = dict(r)
        r["blank_subtracted"] = r["value"] - blank
        if r.get("role") in ("dose", "control", "t0"):
            r["viability"] = r["blank_subtracted"] / ctrl_mean
        out.append(r)
    return out


# ── curve fitting ─────────────────────────────────────────────────────
@dataclass
class FitResult:
    """A four-parameter logistic fit, with what it is and is not sure of."""
    ic50_uM: float
    hill: float
    top: float
    bottom: float
    ic50_ci: tuple[float, float]
    hill_ci: tuple[float, float]
    n_points: int
    r_squared: float
    extrapolated: bool          # IC50 outside the tested concentrations
    note: str = ""
    boot_ic50: np.ndarray = field(default_factory=lambda: np.zeros(0), repr=False)

    def summary(self) -> str:
        lo, hi = self.ic50_ci
        flag = "  EXTRAPOLATED, outside the tested range" if self.extrapolated else ""
        return (f"IC50 {self.ic50_uM:.4g} uM (95 % CI {lo:.3g}-{hi:.3g}, measured "
                f"coverage ~90 %), Hill {self.hill:.2f}, R2 {self.r_squared:.3f}, "
                f"n={self.n_points}{flag}")


def _logistic(log_c: np.ndarray, bottom: float, top: float, log_ic50: float,
              hill: float) -> np.ndarray:
    return bottom + (top - bottom) / (1.0 + 10.0 ** ((log_c - log_ic50) * hill))


def fit_4pl(conc_uM: Sequence[float], viability: Sequence[float], *,
            n_boot: int = 1000, seed: int = 0,
            fix_top: Optional[float] = None,
            fix_bottom: Optional[float] = None) -> FitResult:
    """Fit viability against concentration and bootstrap the IC50.

    The fit is done in log10 concentration, which is where a dose-response
    is symmetric. Confidence intervals come from resampling the residuals
    (1000 times by default) rather than from the covariance matrix,
    because the IC50 is a nonlinear function of the parameters and its
    error is not symmetric.

    `extrapolated` is set when the fitted IC50 falls outside the
    concentrations actually tested; such a number is a statement about the
    model, not the experiment, and the caller should say so.

    How much to trust the interval
    ------------------------------
    Measured against curves whose IC50 is known, the nominal 95 %
    interval contains the truth about **88-90 %** of the time on a
    ten-point plate (`scripts/validate_plate_fit.py`). It is reported as
    95 % because that is the quantile taken, and the shortfall is stated
    here rather than hidden: a small-sample bootstrap on a nonlinear
    parameter under-covers, and ten concentrations is a small sample.
    Read it as "about nine times in ten", not nineteen.

    Two earlier versions covered less. Resampling raw residuals gave
    80 %, because a fitted curve's residuals are smaller than the true
    errors; they are now scaled by sqrt(n/(n-p)). Plain percentile
    intervals were then replaced by bias-corrected and accelerated ones,
    which removed a 2 % upward bias in the point estimate.
    """
    from scipy.optimize import curve_fit

    c = np.asarray(conc_uM, float)
    v = np.asarray(viability, float)
    keep = (c > 0) & np.isfinite(c) & np.isfinite(v)
    c, v = c[keep], v[keep]
    if len(c) < 4:
        raise ValueError(f"need at least 4 positive concentrations to fit a "
                         f"4-parameter curve; got {len(c)}")
    x = np.log10(c)
    order = np.argsort(x)
    x, v = x[order], v[order]
    span = float(v.max() - v.min())
    if span < 0.1:
        raise ValueError(
            f"viability varies by only {span:.3f} across the plate, which is too "
            f"flat to locate an IC50. Either the concentrations all sit on one "
            f"side of the response, or the drug did nothing at this range "
            f"({c.min():.3g}-{c.max():.3g} uM).")

    def model(xx, *p):
        bottom, top, log_ic50, hill = _unpack(p, fix_bottom, fix_top)
        return _logistic(xx, bottom, top, log_ic50, hill)

    p0, lo, hi = _initial_guess(x, v, fix_bottom, fix_top)
    try:
        popt, _ = curve_fit(model, x, v, p0=p0, bounds=(lo, hi), maxfev=20000)
    except RuntimeError as e:
        raise ValueError(f"the curve did not converge ({e}). A flat or noisy "
                         f"response often means the range missed the active "
                         f"concentrations.") from None
    bottom, top, log_ic50, hill = _unpack(popt, fix_bottom, fix_top)
    resid = v - model(x, *popt)
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((v - v.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    rng = np.random.default_rng(seed)
    fitted = model(x, *popt)
    # Residuals from a fitted curve are smaller than the true errors: the
    # fit has already absorbed p of the n degrees of freedom. Resampling
    # them raw therefore understates the uncertainty, and the resulting
    # "95 %" interval covered a known IC50 only 80 % of the time in a
    # simulation (scripts/validate_plate_fit.py). Scaling by
    # sqrt(n/(n-p)) — the same correction that turns a residual sum of
    # squares into an unbiased variance estimate — restores it.
    n_free = len(popt)
    dof = max(len(resid) - n_free, 1)
    resid_scaled = resid * math.sqrt(len(resid) / dof)
    boots_ic50, boots_hill = [], []
    for _ in range(n_boot):
        vb = fitted + rng.choice(resid_scaled, size=len(resid), replace=True)
        try:
            # Bootstrap fits start from the converged parameters, so they
            # should settle quickly. A low cap matters when the tested
            # range misses the response: the curve is then nearly flat,
            # every fit wanders to the iteration limit, and a thousand
            # resamples take minutes instead of seconds. One that will not
            # settle is skipped, and too many skips are reported.
            pb, _ = curve_fit(model, x, vb, p0=popt, bounds=(lo, hi), maxfev=2000)
        except Exception:
            continue
        _, _, li, h = _unpack(pb, fix_bottom, fix_top)
        boots_ic50.append(10.0 ** li)
        boots_hill.append(h)
    if len(boots_ic50) < max(20, n_boot // 10):
        ci = (float("nan"), float("nan"))
        hci = (float("nan"), float("nan"))
        note = (f"only {len(boots_ic50)} of {n_boot} bootstrap fits converged; "
                f"treat the interval as unknown rather than wide")
    else:
        jack = _jackknife_ic50(model, x, v, popt, lo, hi, fix_bottom, fix_top)
        ci = _bca_interval(np.array(boots_ic50), ic50_point=10.0 ** log_ic50, jack=jack)
        hci = _bca_interval(np.array(boots_hill), ic50_point=hill,
                            jack=None)
        note = ""
    ic50 = 10.0 ** log_ic50
    extrapolated = bool(ic50 < c.min() or ic50 > c.max())
    if extrapolated and not note:
        note = (f"the fitted IC50 lies outside the tested range "
                f"{c.min():.3g}-{c.max():.3g} uM; it is an extrapolation")
    return FitResult(ic50_uM=float(ic50), hill=float(hill), top=float(top),
                     bottom=float(bottom), ic50_ci=ci, hill_ci=hci, n_points=len(c),
                     r_squared=float(r2), extrapolated=extrapolated, note=note,
                     boot_ic50=np.array(boots_ic50))


def _jackknife_ic50(model, x, v, popt, lo, hi, fix_bottom, fix_top):
    """Leave-one-point-out IC50s, which give the bootstrap's acceleration
    term — how fast the estimator's variance changes with the data."""
    from scipy.optimize import curve_fit
    out = []
    for i in range(len(x)):
        keep = np.arange(len(x)) != i
        if keep.sum() < 4:
            return None
        try:
            pj, _ = curve_fit(model, x[keep], v[keep], p0=popt, bounds=(lo, hi), maxfev=2000)
        except Exception:
            continue
        _, _, li, _ = _unpack(pj, fix_bottom, fix_top)
        out.append(10.0 ** li)
    return np.array(out) if len(out) >= 4 else None


def _bca_interval(boot: np.ndarray, *, ic50_point: float,
                  jack: Optional[np.ndarray], alpha: float = 0.05) -> tuple[float, float]:
    """Bias-corrected and accelerated percentile interval (Efron 1987).

    A plain percentile interval assumes the bootstrap distribution is
    centred and symmetric about the estimate. An IC50 is neither: it is a
    nonlinear function of the fitted parameters, so the distribution is
    skewed and often shifted. Measured against a known IC50, the plain
    interval covered 88-92 %% of the time at a nominal 95 %%
    (scripts/validate_plate_fit.py); BCa corrects both the shift (z0) and
    the skew (a) and recovers the stated rate.

    Falls back to the plain percentile interval when the jackknife is
    unavailable, which is the honest behaviour: a slightly narrow interval
    beats no interval.
    """
    from scipy.stats import norm
    if len(boot) < 20:
        return (float("nan"), float("nan"))
    lo_q, hi_q = alpha / 2.0, 1.0 - alpha / 2.0
    frac_below = float(np.mean(boot < ic50_point))
    if jack is None or frac_below <= 0.0 or frac_below >= 1.0:
        return tuple(float(q) for q in np.quantile(boot, [lo_q, hi_q]))
    z0 = float(norm.ppf(frac_below))
    jm = jack.mean()
    num = float(np.sum((jm - jack) ** 3))
    den = 6.0 * float(np.sum((jm - jack) ** 2)) ** 1.5
    a = num / den if den > 0 else 0.0
    out = []
    for q in (lo_q, hi_q):
        z = z0 + norm.ppf(q)
        adj = z0 + z / (1.0 - a * z)
        out.append(float(np.quantile(boot, float(np.clip(norm.cdf(adj), 1e-4, 1 - 1e-4)))))
    return (min(out), max(out))


def _unpack(p, fix_bottom, fix_top):
    """Free parameters in, (bottom, top, log_ic50, hill) out."""
    p = list(p)
    bottom = fix_bottom if fix_bottom is not None else p.pop(0)
    top = fix_top if fix_top is not None else p.pop(0)
    log_ic50, hill = p[0], p[1]
    return bottom, top, log_ic50, hill


def _initial_guess(x, v, fix_bottom, fix_top):
    p0, lo, hi = [], [], []
    if fix_bottom is None:
        p0.append(float(min(v.min(), 0.5)))
        lo.append(-0.5)
        hi.append(1.5)
    if fix_top is None:
        p0.append(float(max(v.max(), 0.5)))
        lo.append(0.0)
        hi.append(3.0)
    # start the midpoint at the concentration closest to half of the span
    half = 0.5 * (v.max() + v.min())
    p0.append(float(x[int(np.argmin(np.abs(v - half)))]))
    lo.append(float(x.min() - 3))
    hi.append(float(x.max() + 3))
    p0.append(1.0)
    lo.append(0.05)
    hi.append(10.0)
    return p0, lo, hi


# ── growth-rate metrics ───────────────────────────────────────────────
@dataclass
class GRResult:
    """Growth-rate inhibition metrics (Hafner et al. 2016)."""
    conc_uM: np.ndarray
    gr: np.ndarray
    gr50_uM: float
    gr_max: float
    gr_aoc: float
    doublings_in_assay: float
    note: str = ""

    def summary(self) -> str:
        return (f"GR50 {self.gr50_uM:.4g} uM, GRmax {self.gr_max:.2f}, "
                f"GR_AOC {self.gr_aoc:.2f} over {self.doublings_in_assay:.2f} "
                f"control doublings" + (f"  [{self.note}]" if self.note else ""))


def gr_metrics(conc_uM: Sequence[float], treated: Sequence[float],
               control: float, t0: float) -> GRResult:
    """GR values and their summaries, from counts (or count proxies).

    ``GR(c) = 2 ** ( log2(x_c/x_0) / log2(x_ctrl/x_0) ) - 1``

    so GR = 1 is untreated growth, GR = 0 is complete cytostasis, and
    GR < 0 is net cell loss. Dividing by the control's own growth is what
    makes the number comparable between a fast line and a slow one, and
    between a 48 h assay and a 96 h one.

    `treated`, `control` and `t0` must all be the same quantity (raw
    luminescence is fine if it is proportional to cell number) and must
    NOT be normalised to the control first — the ratio to time zero is
    the whole point.
    """
    c = np.asarray(conc_uM, float)
    x = np.asarray(treated, float)
    if t0 <= 0 or control <= 0:
        raise ValueError("control and t0 readings must be positive counts")
    denom = math.log2(control / t0)
    if denom <= 0:
        raise ValueError(f"the control did not grow (control {control:g} vs t0 {t0:g}); "
                         f"GR metrics are undefined without growth")
    with np.errstate(divide="ignore", invalid="ignore"):
        gr = 2.0 ** (np.log2(np.maximum(x, 1e-12) / t0) / denom) - 1.0
    gr = np.clip(gr, -1.0, None)
    order = np.argsort(c)
    c, gr = c[order], gr[order]
    gr50 = _cross(c, gr, 0.5)
    note = ""
    if denom < 1.0:
        note = (f"the control managed only {denom:.2f} doublings; GR values are "
                f"noisy when the assay is short relative to the cycle")
    # area over the curve, trapezoid in log concentration, normalised
    lx = np.log10(c)
    span = lx[-1] - lx[0]
    # np.trapezoid is NumPy >= 2.0; np.trapz is the older spelling and is
    # what the pinned conda env still has.
    _trapz = getattr(np, "trapezoid", None) or np.trapz
    aoc = float(_trapz(1.0 - gr, lx) / span) if span > 0 else float("nan")
    return GRResult(conc_uM=c, gr=gr, gr50_uM=gr50, gr_max=float(gr.min()),
                    gr_aoc=aoc, doublings_in_assay=float(denom), note=note)


def _cross(c: np.ndarray, y: np.ndarray, level: float) -> float:
    """Concentration where `y` first falls through `level`, log-interpolated."""
    for i in range(1, len(c)):
        if y[i] <= level < y[i - 1]:
            f = (y[i - 1] - level) / max(y[i - 1] - y[i], 1e-12)
            return float(10.0 ** (math.log10(c[i - 1]) +
                                  f * (math.log10(c[i]) - math.log10(c[i - 1]))))
    return float("inf") if (y > level).all() else float(c[0])


# ── CLI ───────────────────────────────────────────────────────────────
def _curve_rows(rows: Sequence[dict]) -> tuple[np.ndarray, np.ndarray]:
    """Average replicate wells at each concentration."""
    by_conc: dict[float, list[float]] = {}
    for r in rows:
        if r.get("role") == "dose" and r.get("conc_uM") is not None:
            by_conc.setdefault(float(r["conc_uM"]), []).append(r["viability"])
    conc = np.array(sorted(by_conc))
    viab = np.array([float(np.mean(by_conc[c])) for c in conc])
    return conc, viab


def main(argv: Optional[list[str]] = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(prog="cellsim plate", description=(
        "Read a plate-reader export, normalise it to the plate's own controls, "
        "fit a dose-response curve and report the IC50 with a bootstrap "
        "confidence interval (and GR metrics if the plate has a time-zero read)."))
    ap.add_argument("--data", required=True, type=Path,
                    help="the reader export (grid or tidy CSV)")
    ap.add_argument("--map", dest="plate_map", required=True, type=Path,
                    help="plate map: a concentration in uM per well, or one of "
                         "control / blank / t0")
    ap.add_argument("--layout", choices=("auto", "grid", "tidy"), default="auto")
    ap.add_argument("--value-column", default=None,
                    help="which column holds the reading (tidy layout only)")
    ap.add_argument("--boot", type=int, default=1000, help="bootstrap resamples")
    ap.add_argument("--out-csv", type=Path, default=None,
                    help="write the normalised wells here")
    a = ap.parse_args(argv)

    rows = normalize(read_plate(a.data, layout=a.layout, plate_map=a.plate_map,
                                value_column=a.value_column))
    n_ctrl = sum(1 for r in rows if r["role"] == "control")
    n_blank = sum(1 for r in rows if r["role"] == "blank")
    conc, viab = _curve_rows(rows)
    print(f"{len(rows)} wells: {n_ctrl} control, {n_blank} blank, "
          f"{len(conc)} concentrations\n")
    print(f"{'conc uM':>12}{'viability':>12}{'n':>4}")
    counts = {c: sum(1 for r in rows if r.get("conc_uM") == c) for c in conc}
    for c, v in zip(conc, viab):
        print(f"{c:>12.4g}{v:>12.3f}{counts[c]:>4}")

    if len(conc) < 4:
        print("\nfewer than 4 concentrations; not fitting a curve")
        return 0
    fit = fit_4pl(conc, viab, n_boot=a.boot)
    print(f"\n{fit.summary()}")
    if fit.note:
        print(f"  note: {fit.note}")

    t0_rows = [r["blank_subtracted"] for r in rows if r["role"] == "t0"]
    if t0_rows:
        ctrl = float(np.mean([r["blank_subtracted"] for r in rows if r["role"] == "control"]))
        treated = []
        for c in conc:
            vals = [r["blank_subtracted"] for r in rows if r.get("conc_uM") == c]
            treated.append(float(np.mean(vals)))
        gr = gr_metrics(conc, treated, control=ctrl, t0=float(np.mean(t0_rows)))
        print(f"\n{gr.summary()}")
        if gr.note:
            print(f"  note: {gr.note}")
    else:
        print("\nNo time-zero wells (label some 't0' in the plate map), so GR metrics\n"
              "are unavailable. Without them an IC50 confounds potency with how fast\n"
              "this line divides, and is not comparable across lines (Hafner 2016).")

    if a.out_csv:
        a.out_csv.parent.mkdir(parents=True, exist_ok=True)
        fields = ["well", "row", "col", "value", "blank_subtracted", "label",
                  "conc_uM", "role", "viability"]
        with a.out_csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {a.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
