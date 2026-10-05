#!/usr/bin/env python3
"""Phase-4 gates: reading a real plate, and fitting it honestly.

This is the module a lab touches first, so the gates are about the ways
a plate actually arrives rather than about the maths alone: a reader's
grid export, a tidy table, zero-padded well names, blanks that must be
subtracted, and controls that must set the scale.

The fitting gates check the two claims `fit_4pl` makes — that it recovers
an IC50 it was not told, and that its interval means what it says — on
curves whose truth is known. The coverage claim is measured properly by
`scripts/validate_plate_fit.py`; here it is sampled just enough to catch
a regression.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.plate import fit_4pl, gr_metrics, normalize, read_plate  # noqa: E402

TRUE_IC50, TRUE_HILL = 1.7, 1.3
DOSES = [0.0137, 0.0411, 0.123, 0.370, 1.111, 3.333, 10.0, 30.0]


def _viability(c: float) -> float:
    return 0.03 + 0.97 / (1.0 + 10.0 ** ((np.log10(c) - np.log10(TRUE_IC50)) * TRUE_HILL))


def _write_plate(d: Path, *, padded: bool = False, tidy: bool = False) -> tuple[Path, Path]:
    """A 96-well plate: row A blanks, row B controls, rows C-E a dose
    series in triplicate. Signal = 1000 * viability + 50 background."""
    signal: dict[str, float] = {}
    label: dict[str, str] = {}
    for col in range(1, 13):
        signal[f"A{col}"] = 50.0
        label[f"A{col}"] = "blank"
        signal[f"B{col}"] = 50.0 + 1000.0
        label[f"B{col}"] = "control"
    for r, row in enumerate("CDE"):
        for i, dose in enumerate(DOSES):
            well = f"{row}{i + 1}"
            jitter = 1.0 + 0.01 * ((r + i) % 3 - 1)      # deterministic, tiny
            signal[well] = 50.0 + 1000.0 * _viability(dose) * jitter
            label[well] = f"{dose:g}"
    data, mapping = d / "reading.csv", d / "map.csv"
    if tidy:
        def fmt(w):
            return f"{w[0]}{int(w[1:]):02d}" if padded else w
        data.write_text("Well,Luminescence\n" +
                        "".join(f"{fmt(w)},{v:.2f}\n" for w, v in signal.items()))
        mapping.write_text("well,condition\n" +
                           "".join(f"{fmt(w)},{s}\n" for w, s in label.items()))
    else:
        def grid(values, default=""):
            out = ["," + ",".join(str(c) for c in range(1, 13))]
            for row in "ABCDEFGH":
                cells = [str(values.get(f"{row}{c}", default)) for c in range(1, 13)]
                out.append(row + "," + ",".join(cells))
            return "\n".join(out) + "\n"
        data.write_text(grid({k: f"{v:.2f}" for k, v in signal.items()}, default="0"))
        mapping.write_text(grid(label, default="empty"))
    return data, mapping


def test_reads_a_reader_grid_export():
    with tempfile.TemporaryDirectory() as d:
        data, mapping = _write_plate(Path(d))
        rows = read_plate(data, plate_map=mapping)
    assert len(rows) == 96, len(rows)
    # A grid export pads the unused wells; those must come back as
    # 'empty', NOT as blanks, or the padding is averaged into the
    # background and silently rescales every viability on the plate.
    assert {r["role"] for r in rows} == {"blank", "control", "dose", "empty"}
    assert sum(1 for r in rows if r["role"] == "blank") == 12
    doses = sorted({r["conc_uM"] for r in rows if r["role"] == "dose"})
    assert len(doses) == len(DOSES)
    assert sum(1 for r in rows if r["role"] == "control") == 12


def test_reads_a_tidy_table_with_padded_well_names():
    """'B07' from the reader and 'B7' from the plate map must line up, or
    every well silently fails to match.

    A tidy export lists only the wells that exist (48 here: 24 blanks and
    controls, 24 dosed), where the grid export pads the unused rows. Both
    are normal, which is why the reader keys on well names rather than on
    position in a block."""
    with tempfile.TemporaryDirectory() as d:
        data, mapping = _write_plate(Path(d), padded=True, tidy=True)
        rows = read_plate(data, plate_map=mapping)
    assert len(rows) == 48, len(rows)
    assert all(r["well"] == f"{r['row']}{r['col']}" for r in rows)
    assert any(r["well"] == "B7" for r in rows), "zero padding was not normalised"
    assert sum(1 for r in rows if r["role"] == "dose") == 24


def test_an_unmapped_well_is_an_error_not_a_silent_drop():
    with tempfile.TemporaryDirectory() as d:
        data, mapping = _write_plate(Path(d), tidy=True)
        short = Path(d) / "short_map.csv"
        short.write_text("\n".join(mapping.read_text().splitlines()[:40]) + "\n")
        try:
            read_plate(data, plate_map=short)
        except ValueError as e:
            assert "no entry" in str(e), e
        else:
            raise AssertionError("a plate map missing wells must raise")


def test_normalisation_subtracts_blanks_and_sets_controls_to_one():
    with tempfile.TemporaryDirectory() as d:
        data, mapping = _write_plate(Path(d))
        rows = normalize(read_plate(data, plate_map=mapping))
    ctrl = [r["viability"] for r in rows if r["role"] == "control"]
    assert abs(float(np.mean(ctrl)) - 1.0) < 1e-9, float(np.mean(ctrl))
    top = [r["viability"] for r in rows if r.get("conc_uM") == min(DOSES)]
    bot = [r["viability"] for r in rows if r.get("conc_uM") == max(DOSES)]
    assert float(np.mean(top)) > 0.9, float(np.mean(top))
    assert float(np.mean(bot)) < 0.15, float(np.mean(bot))


def test_a_plate_without_controls_is_refused():
    with tempfile.TemporaryDirectory() as d:
        data, mapping = _write_plate(Path(d), tidy=True)
        rows = [r for r in read_plate(data, plate_map=mapping) if r["role"] != "control"]
        try:
            normalize(rows)
        except ValueError as e:
            assert "control" in str(e)
        else:
            raise AssertionError("normalising without controls must raise")


def test_the_fit_recovers_an_ic50_it_was_not_told():
    with tempfile.TemporaryDirectory() as d:
        data, mapping = _write_plate(Path(d))
        rows = [r for r in normalize(read_plate(data, plate_map=mapping))
                if r["role"] == "dose"]
    f = fit_4pl([r["conc_uM"] for r in rows], [r["viability"] for r in rows], n_boot=200)
    assert 0.7 * TRUE_IC50 < f.ic50_uM < 1.3 * TRUE_IC50, f.summary()
    assert f.r_squared > 0.98, f.summary()
    assert not f.extrapolated
    assert f.ic50_ci[0] <= f.ic50_uM <= f.ic50_ci[1]


def test_a_range_that_misses_the_ic50_is_flagged_as_extrapolated():
    """The failure mode worth catching: the fit still returns a number and
    a narrow interval, and both are about the model, not the plate."""
    conc = np.geomspace(0.001, 0.05, 8)          # entirely below the IC50
    viab = np.array([_viability(c) for c in conc])
    f = fit_4pl(conc, viab, n_boot=200)
    assert f.extrapolated, f.summary()
    assert "extrapolation" in f.note.lower(), f.note


def test_the_interval_covers_a_known_ic50_about_as_often_as_claimed():
    """Sampled coverage. The proper measurement is
    scripts/validate_plate_fit.py; this catches a regression that would
    make the interval badly wrong."""
    rng = np.random.default_rng(3)
    conc = np.geomspace(0.01, 100, 10)
    clean = np.array([_viability(c) for c in conc])
    hits = n = 0
    for t in range(40):
        f = fit_4pl(conc, clean + rng.normal(0, 0.05, len(conc)), n_boot=200, seed=t)
        if not np.isfinite(f.ic50_ci[0]):
            continue
        n += 1
        hits += f.ic50_ci[0] <= TRUE_IC50 <= f.ic50_ci[1]
    assert n >= 35, f"only {n} of 40 fits produced an interval"
    assert hits / n >= 0.80, f"95 % interval covered {hits}/{n} = {hits / n:.0%}"


def test_gr_metrics_separate_cytostasis_from_cell_killing():
    """GR = 0 is a drug that stops growth; GR < 0 is one that kills. A
    viability readout alone cannot tell those apart."""
    stasis = gr_metrics([1.0, 10.0], [1000.0, 1000.0], control=4000.0, t0=1000.0)
    killing = gr_metrics([1.0, 10.0], [1000.0, 250.0], control=4000.0, t0=1000.0)
    assert abs(stasis.gr[-1]) < 1e-9, stasis.gr
    assert killing.gr[-1] < -0.4, killing.gr
    assert abs(stasis.doublings_in_assay - 2.0) < 1e-9


def test_gr_corrects_for_how_fast_the_line_divides():
    """The point of GR. Two lines, same drug effect per division, but one
    doubles twice in the assay and the other once. Viability makes the
    fast line look more sensitive; GR does not."""
    fast = gr_metrics([10.0], [2000.0], control=4000.0, t0=1000.0)   # 2 doublings
    slow = gr_metrics([10.0], [1414.0], control=2000.0, t0=1000.0)   # 1 doubling
    viab_fast, viab_slow = 2000.0 / 4000.0, 1414.0 / 2000.0
    assert viab_fast < viab_slow - 0.1, "the setup should differ on viability"
    assert abs(fast.gr[0] - slow.gr[0]) < 0.02, (fast.gr[0], slow.gr[0])


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} plate gates pass")
    raise SystemExit(1 if failed else 0)
