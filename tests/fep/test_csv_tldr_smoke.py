"""csv_tldr.py regression tests.

Auto-detect kind, format summary, set verdict, pick exit code. The
whole point of the script is that it fits on one Slack line AND
that the exit code is honest (CI can gate on it). The tests pin
both surface behaviours at once."""
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))


def _run(csv_path: Path, *extra: str) -> tuple[int, str]:
    r = subprocess.run(
        ["python", "scripts/csv_tldr.py", str(csv_path), *extra],
        cwd=REPO, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def _write_csv(dir_path: Path, rows: list[dict],
               name: str = "results.csv") -> Path:
    p = dir_path / name
    with p.open("w", encoding="utf-8", newline="") as fo:
        writer = csv.DictWriter(fo, fieldnames=[
            "name", "smiles", "dG_expt_kcalmol",
            "dG_pred_kcalmol", "uncertainty_kcalmol",
            "residual_kcalmol", "wall_seconds", "ok", "reason"])
        writer.writeheader()
        for r in rows:
            r.setdefault("uncertainty_kcalmol", "0.4")
            r.setdefault("residual_kcalmol", "0.0")
            r.setdefault("wall_seconds", "100")
            r.setdefault("ok", "True")
            r.setdefault("reason", "")
            writer.writerow(r)
    return p


def test_hydration_pass_via_bundled_fixture():
    rc, out = _run(REPO / "tests/fep/fixtures/ok_case/"
                          "freesolv_results.csv")
    assert rc == 0, out
    assert "FreeSolv hydration FEP" in out
    assert "12/12 ok" in out
    assert "PASS" in out
    assert "Pearson r" in out


def test_hydration_fail_via_bundled_fixture():
    """fail_case has methane sign-flipped + MAE blown out."""
    rc, out = _run(REPO / "tests/fep/fixtures/fail_case/"
                          "freesolv_results.csv")
    assert rc == 1, out
    assert "FAIL" in out
    assert "sign flip" in out


def test_binding_streptavidin_pass_label_correct():
    with tempfile.TemporaryDirectory(prefix="tldr_") as tmp:
        p = _write_csv(Path(tmp), [
            {"name": "biotin", "smiles": "O=C1N...",
             "dG_expt_kcalmol": "-18.30",
             "dG_pred_kcalmol": "-17.90"},
            {"name": "desthiobiotin", "smiles": "...",
             "dG_expt_kcalmol": "-13.20",
             "dG_pred_kcalmol": "-12.50"},
            {"name": "2-iminobiotin", "smiles": "...",
             "dG_expt_kcalmol": "-12.00",
             "dG_pred_kcalmol": "-11.40"},
            {"name": "biotin_methyl_ester", "smiles": "...",
             "dG_expt_kcalmol": "-10.00",
             "dG_pred_kcalmol": "-9.50"},
        ])
        rc, out = _run(p)
    assert rc == 0, out
    assert "streptavidin binding FEP" in out, out
    assert "gate ≤2.0" in out, out
    assert "Kendall τ" in out
    assert "PASS" in out


def test_binding_egfr_label_selected():
    with tempfile.TemporaryDirectory(prefix="tldr_") as tmp:
        p = _write_csv(Path(tmp), [
            {"name": "erlotinib", "smiles": "...",
             "dG_expt_kcalmol": "-11.86",
             "dG_pred_kcalmol": "-11.42"},
            {"name": "gefitinib", "smiles": "...",
             "dG_expt_kcalmol": "-10.20",
             "dG_pred_kcalmol": "-10.01"},
            {"name": "lapatinib", "smiles": "...",
             "dG_expt_kcalmol": "-12.50",
             "dG_pred_kcalmol": "-12.90"},
        ])
        rc, out = _run(p)
    assert rc == 0, out
    assert "EGFR kinase binding FEP" in out, out


def test_binding_non_binder_flagged():
    with tempfile.TemporaryDirectory(prefix="tldr_") as tmp:
        p = _write_csv(Path(tmp), [
            {"name": "biotin", "smiles": "...",
             "dG_expt_kcalmol": "-18.30",
             "dG_pred_kcalmol": "-17.90"},
            {"name": "desthiobiotin", "smiles": "...",
             "dG_expt_kcalmol": "-13.20",
             # Predicted non-binder (+0.5) — must FAIL.
             "dG_pred_kcalmol": "+0.50"},
            {"name": "2-iminobiotin", "smiles": "...",
             "dG_expt_kcalmol": "-12.00",
             "dG_pred_kcalmol": "-11.40"},
        ])
        rc, out = _run(p)
    assert rc == 1, out
    assert "FAIL" in out
    assert "non-binder" in out, out


def test_inconclusive_when_some_rows_not_ok():
    """Partial: one compound errored out mid-run — verdict must be
    inconclusive (exit 2) not PASS."""
    with tempfile.TemporaryDirectory(prefix="tldr_") as tmp:
        p = _write_csv(Path(tmp), [
            {"name": "methane", "smiles": "C",
             "dG_expt_kcalmol": "2.00",
             "dG_pred_kcalmol": "1.72"},
            {"name": "acetamide", "smiles": "CC(=O)N",
             "dG_expt_kcalmol": "-9.71",
             "dG_pred_kcalmol": "",
             "ok": "False",
             "reason": "MBAR failed"},
        ])
        rc, out = _run(p)
    assert rc == 2, out
    assert "inconclusive" in out, out
    assert "1/2" in out


def test_no_such_file_exits_3():
    rc, out = _run(Path("/nonexistent/path/asdf.csv"))
    assert rc == 3, out
    assert "no such file" in out.lower() or "not found" in out.lower()


def test_gate_override_turns_fail_into_pass():
    """Same fail_case fixture, but with a relaxed gate the MAE
    passes — except sign-flip still FAILs the binding/hydration
    sign rule, so verdict is still FAIL. This pins that the gate
    override doesn't silently mask the sign rule."""
    rc, out = _run(REPO / "tests/fep/fixtures/fail_case/"
                          "freesolv_results.csv",
                   "--gate", "5.0")
    assert rc == 1, out
    assert "FAIL" in out
    assert "sign flip" in out


def test_output_is_single_line():
    """TL;DR = one line. Pin that no stray multi-line output
    sneaks in."""
    rc, out = _run(REPO / "tests/fep/fixtures/ok_case/"
                          "freesolv_results.csv")
    # stdout should be one line (+trailing newline). stderr empty.
    lines = [line for line in out.splitlines() if line.strip()]
    assert len(lines) == 1, (
        f"expected 1 non-empty output line; got {len(lines)}:\n{out}")


if __name__ == "__main__":
    funcs = [
        test_hydration_pass_via_bundled_fixture,
        test_hydration_fail_via_bundled_fixture,
        test_binding_streptavidin_pass_label_correct,
        test_binding_egfr_label_selected,
        test_binding_non_binder_flagged,
        test_inconclusive_when_some_rows_not_ok,
        test_no_such_file_exits_3,
        test_gate_override_turns_fail_into_pass,
        test_output_is_single_line,
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
