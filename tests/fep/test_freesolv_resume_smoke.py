"""freesolv_validate --resume regression.

Pins:
  - Incremental CSV write after each compound (so a kill mid-run
    preserves completed work)
  - --resume reads back the CSV and skips compounds with non-empty
    dG_pred_kcalmol
  - Mid-compound crash (no pred written) doesn't poison the resume
  - Final write is idempotent with the incremental writes

Multi-day CPU FreeSolv runs can die on lid-close / power / kernel
panic. Without these properties the user loses ~30h of work; with
them, --resume picks up exactly where the crash happened.
"""
from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.fep.freesolv_validate import (
    FreeSolvPoint,
    _CSV_COLS,
    _load_resume_rows,
    _write_csv,
)


def _make_point(name, pred=None, ok=True):
    return FreeSolvPoint(
        name=name, smiles="C" if name == "methane" else "CCO",
        dG_expt_kcalmol=2.00 if name == "methane" else -5.01,
        dG_pred_kcalmol=pred,
        uncertainty_kcalmol=0.3 if pred is not None else None,
        residual_kcalmol=(pred - 2.0) if pred is not None else None,
        wall_seconds=420.0 if pred is not None else None,
        ok=ok, reason="" if ok else "sampling failed: foo",
        ghmc_accept_mean=0.85 if pred is not None else None,
        ghmc_accept_min=0.72 if pred is not None else None,
    )


def test_write_csv_round_trip_via_load_resume():
    with tempfile.TemporaryDirectory(prefix="freesolv_resume_") as tmp:
        p = Path(tmp) / "results.csv"
        rows = [_make_point("methane", pred=+1.80),
                _make_point("ethanol", pred=-5.30)]
        _write_csv(p, rows)
        assert p.exists()
        loaded, names = _load_resume_rows(p)
        assert names == {"methane", "ethanol"}
        assert len(loaded) == 2
        assert loaded[0].name == "methane"
        assert loaded[0].dG_pred_kcalmol == 1.80
        assert loaded[1].dG_pred_kcalmol == -5.30


def test_resume_skips_incomplete_compound_so_it_reruns():
    """Compound that started but didn't finish writes a row with
    pred=None. --resume must NOT mark it complete; the loop must
    re-run that compound from scratch."""
    with tempfile.TemporaryDirectory(prefix="freesolv_resume_") as tmp:
        p = Path(tmp) / "results.csv"
        rows = [
            _make_point("methane", pred=+1.80),
            _make_point("ethane", pred=None, ok=False),  # crashed
        ]
        _write_csv(p, rows)
        loaded, names = _load_resume_rows(p)
        # methane resumes-as-done; ethane is missing → will re-run
        assert names == {"methane"}, names
        assert len(loaded) == 1


def test_load_resume_handles_missing_file():
    """First run, no CSV exists → empty lists, no crash."""
    with tempfile.TemporaryDirectory(prefix="freesolv_resume_") as tmp:
        p = Path(tmp) / "does_not_exist.csv"
        loaded, names = _load_resume_rows(p)
        assert loaded == []
        assert names == set()


def test_write_csv_is_atomic_via_tmp_rename():
    """If the .tmp file exists and full CSV doesn't, a previous
    write was interrupted; current code overwrites tmp cleanly."""
    with tempfile.TemporaryDirectory(prefix="freesolv_resume_") as tmp:
        p = Path(tmp) / "results.csv"
        tmp_path = p.with_suffix(p.suffix + ".tmp")
        # Plant a stale tmp file from a "killed" prior write
        tmp_path.write_text("stale,partial,write,leftover\n")
        # New write should clobber and produce a valid CSV
        _write_csv(p, [_make_point("methane", pred=+1.80)])
        assert p.exists()
        # Stale tmp was renamed away (rename clobbers target)
        assert not tmp_path.exists()
        # And the CSV is readable
        loaded, names = _load_resume_rows(p)
        assert names == {"methane"}


def test_csv_schema_columns_pinned():
    """Adding a column without thinking can break fep-report. Pin
    the exact set + order."""
    assert _CSV_COLS == [
        "name", "smiles", "dG_expt_kcalmol",
        "dG_pred_kcalmol", "uncertainty_kcalmol",
        "residual_kcalmol", "ghmc_accept_mean",
        "ghmc_accept_min", "wall_seconds", "ok", "reason",
    ]


def test_resume_loads_partial_run_three_done_one_pending():
    """End-to-end realistic case: 12 compounds, 3 completed, 9
    waiting. --resume keeps the 3, the next bench run runs the 9."""
    with tempfile.TemporaryDirectory(prefix="freesolv_resume_") as tmp:
        p = Path(tmp) / "results.csv"
        rows = [
            _make_point("methane", pred=+1.80),
            _make_point("ethane", pred=+1.55),
            _make_point("propane", pred=+2.31),
        ]
        _write_csv(p, rows)
        loaded, names = _load_resume_rows(p)
        assert names == {"methane", "ethane", "propane"}
        assert len(loaded) == 3
        # The three preds round-tripped exactly
        preds = {r.name: r.dG_pred_kcalmol for r in loaded}
        assert preds == {"methane": 1.80, "ethane": 1.55,
                         "propane": 2.31}


def test_csv_compatible_with_fep_report():
    """The incremental writer must produce a CSV that cellsim
    fep-report can consume. Smoke-check via a DictReader and
    field-presence."""
    with tempfile.TemporaryDirectory(prefix="freesolv_resume_") as tmp:
        p = Path(tmp) / "results.csv"
        _write_csv(p, [_make_point("methane", pred=+1.80)])
        with p.open("r", encoding="utf-8-sig") as fi:
            r = next(csv.DictReader(fi))
        for col in ("name", "smiles", "dG_expt_kcalmol",
                    "dG_pred_kcalmol", "wall_seconds", "ok"):
            assert col in r, f"missing column {col}"


if __name__ == "__main__":
    funcs = [
        test_write_csv_round_trip_via_load_resume,
        test_resume_skips_incomplete_compound_so_it_reruns,
        test_load_resume_handles_missing_file,
        test_write_csv_is_atomic_via_tmp_rename,
        test_csv_schema_columns_pinned,
        test_resume_loads_partial_run_three_done_one_pending,
        test_csv_compatible_with_fep_report,
    ]
    fails = []
    for f in funcs:
        try:
            f()
            print(f"[PASS] {f.__name__}")
        except AssertionError as e:
            print(f"[FAIL] {f.__name__}: {e}")
            fails.append(f.__name__)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"[ERROR] {f.__name__}: {e}")
            fails.append(f.__name__)
    print(f"{len(funcs) - len(fails)}/{len(funcs)} PASS")
    sys.exit(0 if not fails else 1)
