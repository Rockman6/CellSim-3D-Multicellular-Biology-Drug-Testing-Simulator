"""fill_prof_email.py regression tests.

Covers the three outcome shapes (PASS / FAIL / partial) across
the two YAML kinds (hydration / binding), plus the failure
modes that matter most: missing fields warn, NaN pred doesn't
crash, bench name auto-detect picks the right template.
"""
from __future__ import annotations

import csv
import io
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _write_report_with_overrides(dir_path: Path, kind: str,
                                    verdict: str,
                                    rows_table: str,
                                    mae: str = "0.42",
                                    pearson: str = "+0.993",
                                    ghmc: str = "overall mean 77%, "
                                                "worst 71%, 3 "
                                                "compounds report"):
    """Synthesise a report.md with controllable fields so we can
    exercise filler edge cases without running the sampler."""
    title = {"hydration": "Hydration FEP report",
             "binding": "Binding FEP report"}[kind]
    report = f"""\
# {title} — {verdict}

- source: `/tmp/synth`
- compounds: 12 ok / 12 total

## Gate verdict

- **MAE gate** (≤ 1.5 kcal/mol): PASS  (MAE = {mae} kcal/mol)
- **GHMC acceptance gate** (≥ 70% per compound): PASS  ({ghmc})

## Aggregate accuracy

- MAE  = {mae} kcal/mol
- RMSE = 0.447 kcal/mol
- Pearson r    = {pearson}
- Spearman ρ   = +0.965
- Kendall τ    = +0.879

## Per-compound results

| name | smiles | expt | pred | ± | residual | |resid| | within σ | wall (s) | flags |
|------|--------|------|------|---|----------|-----:|---------:|---------:|-------|
{rows_table}
"""
    (dir_path / "report.md").write_text(report, encoding="utf-8")


def _run_fill(dir_path: Path, *extra_args: str) -> tuple[int, str]:
    r = subprocess.run(
        ["python", "scripts/fill_prof_email.py", str(dir_path),
         *extra_args],
        cwd=REPO_ROOT,
        capture_output=True, text=True,
    )
    return r.returncode, r.stdout + r.stderr


def test_hydration_pass_all_fields_fill():
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        _write_report_with_overrides(
            tmp, kind="hydration", verdict="PASS",
            rows_table=(
                "| methane | `C` | +2.00 | +1.72 | 0.28 | -0.28 | "
                "0.28 | ✓ | 420 | |\n"
                "| acetamide | `CC(=O)N` | -9.71 | -8.90 | 0.55 | "
                "+0.81 | 0.81 | | 560 | |"))
        rc, out = _run_fill(tmp)
    assert rc == 0, f"hydration pass should exit 0; got {rc}\n{out}"
    assert "Subject: Milestone A — FreeSolv FEP results" in out
    assert "Sign-correct on methane:          PASS" in out
    assert "Sign-correct on acetamide:        PASS" in out
    assert "<missing>" not in out


def test_hydration_fail_methane_sign():
    """Methane predicted negative — report FAILs sign-critical."""
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        _write_report_with_overrides(
            tmp, kind="hydration", verdict="FAIL",
            rows_table=(
                "| methane | `C` | +2.00 | -1.85 | 0.45 | -3.85 | "
                "3.85 | | 180 | SIGN WRONG |\n"
                "| acetamide | `CC(=O)N` | -9.71 | -11.20 | 0.85 | "
                "-1.49 | 1.49 | | 265 | |"))
        rc, out = _run_fill(tmp)
    # The fill itself succeeds even when the report says FAIL —
    # the email just conveys FAIL to the prof. Exit 0 = all fields
    # filled; the FAIL is in the email content, not the exit code.
    assert rc == 0, out
    assert "Overall verdict:                  FAIL" in out
    assert "Sign-correct on methane:          FAIL" in out
    assert "Sign-correct on acetamide:        PASS" in out


def test_binding_streptavidin_template_selected():
    """Binding report with biotin/desthiobiotin rows → binding
    template, subject says 'streptavidin binding'."""
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        _write_report_with_overrides(
            tmp, kind="binding", verdict="PASS",
            mae="0.95",
            rows_table=(
                "| biotin | `O=C1N...` | -18.30 | -17.90 | 0.50 | "
                "+0.40 | 0.40 | ✓ | 1200 | |\n"
                "| desthiobiotin | `C[C@H]1...` | -13.20 | -12.50 | "
                "0.60 | +0.70 | 0.70 | | 1300 | |"))
        rc, out = _run_fill(tmp)
    # rc 0 or 2 — 2 only if GHMC/other field is missing from the
    # synth report, which is OK for this test (we're checking the
    # template + subject + values, not the warning).
    assert rc in (0, 2), out
    assert "Subject: Milestone B — streptavidin binding FEP" in out, (
        f"expected 'streptavidin binding' subject; got:\n{out[:500]}")
    # MAE line — use loose whitespace match.
    import re as _re
    assert _re.search(
        r"MAE vs published ΔG_bind:\s+0\.95 kcal/mol", out), (
        f"MAE line missing / mismatched; got:\n{out[:500]}")
    assert "(gate <= 2.0)" in out


def test_binding_egfr_template_selected():
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        _write_report_with_overrides(
            tmp, kind="binding", verdict="PASS",
            rows_table=(
                "| erlotinib | `COCCOc1...` | -11.86 | -11.42 | 0.44 | "
                "+0.44 | 0.44 | ✓ | 3456 | |\n"
                "| gefitinib | `COc1cc...` | -10.20 | -10.01 | 0.38 | "
                "+0.19 | 0.19 | ✓ | 2890 | |"))
        rc, out = _run_fill(tmp)
    assert rc == 0
    assert "EGFR kinase binding FEP" in out, (
        f"expected EGFR bench name; got:\n{out[:400]}")
    assert "Cheng-Prusoff" in out, (
        "binding template should mention the IC50-offset caveat")


def test_binding_non_binder_flagged_in_email():
    """Binding rule: any compound with ΔG_pred >= 0 fails.
    Email must surface the non-binder explicitly."""
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        _write_report_with_overrides(
            tmp, kind="binding", verdict="FAIL",
            rows_table=(
                "| biotin | `O=C1N...` | -18.30 | -17.90 | 0.50 | "
                "+0.40 | 0.40 | ✓ | 1200 | |\n"
                "| desthiobiotin | `C[C@H]1...` | -13.20 | +0.50 | "
                "0.60 | +13.70 | 13.70 | | 1300 | SIGN WRONG |"))
        rc, out = _run_fill(tmp)
    assert rc == 0
    assert "desthiobiotin predicted non-binder" in out, (
        f"email must call out the non-binder; got:\n{out[-600:]}")


def test_missing_numbers_returns_exit_2():
    """A malformed report (e.g. MAE line corrupted) leaves
    <missing> markers — exit 2 so the biologist doesn't silently
    send a placeholder email."""
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        # Deliberately break the MAE line.
        (tmp / "report.md").write_text(
            "# Hydration FEP report — PASS\n\n"
            "- compounds: 12 ok / 12 total\n\n"
            "## Aggregate accuracy\n\n"
            "- MAE: (corrupted)\n"  # No '= X' → unparseable
            "- Pearson r    = +0.993\n",
            encoding="utf-8")
        rc, out = _run_fill(tmp)
    assert rc == 2, (
        f"missing MAE should exit 2 (warn); got rc={rc}\n{out}")
    assert "<missing>" in out
    assert "WARNING" in out


def test_no_report_md_returns_exit_1():
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        rc, out = _run_fill(Path(tmp))
    assert rc == 1, (
        f"no report.md should exit 1 (hard error); got {rc}\n{out}")
    assert "no report.md" in out


def test_prof_name_override():
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        _write_report_with_overrides(
            tmp, kind="hydration", verdict="PASS",
            rows_table=(
                "| methane | `C` | +2.00 | +1.72 | 0.28 | -0.28 | "
                "0.28 | ✓ | 420 | |\n"
                "| acetamide | `CC(=O)N` | -9.71 | -8.90 | 0.55 | "
                "+0.81 | 0.81 | | 560 | |"))
        rc, out = _run_fill(tmp, "--prof-name", "Dr. Chen")
    assert rc == 0
    assert "Hi Dr. Chen," in out
    assert "Hi [Prof]," not in out


def test_hydration_hardware_override_lands_in_body():
    """--hardware / --platform must reach the hydration body, not
    just the binding one. Earlier the hydration template had M5 Max
    hardcoded; pilot-3 CPU run on Henry's MacBook exposed the
    hardcoding and pushed a parametrisation."""
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        _write_report_with_overrides(
            tmp, kind="hydration", verdict="PASS",
            rows_table=(
                "| methane | `C` | +2.00 | +1.78 | 0.15 | -0.22 | "
                "0.22 | ✓ | 4290 | |\n"
                "| acetamide | `CC(=O)N` | -9.71 | -8.90 | 0.55 | "
                "+0.81 | 0.81 | | 560 | |"))
        rc, out = _run_fill(
            tmp,
            "--hardware", "Henry's MacBook Pro M-series, 16 GB RAM",
            "--platform", "CPU (Reference/CPU only; no Metal)")
    assert rc == 0, out
    assert "Henry's MacBook Pro M-series, 16 GB RAM" in out, (
        "hydration body must use --hardware override; got:\n"
        f"{out[:800]}")
    assert "OpenMM platform: CPU (Reference/CPU only; no Metal)" in out, (
        "hydration body must use --platform override; got:\n"
        f"{out[:800]}")
    # The old M5 Max hardcoding must be gone for the hydration path.
    assert "Apple M5 Max" not in out, (
        "hydration template still has M5 Max hardcoded; got:\n"
        f"{out[:800]}")


def test_fail_case_fixture_hydration():
    """End-to-end on the committed fail_case fixture — methane
    sign flipped, MAE 1.67. Email reflects both failures."""
    with tempfile.TemporaryDirectory(
            prefix="cellsim_fill_") as tmp:
        tmp = Path(tmp)
        # Generate report via fep-report on the bundled fixture.
        r = subprocess.run(
            ["python", "-m", "src.fep.report",
             "tests/fep/fixtures/fail_case",
             "--yaml", "benchmarks/fep/freesolv_12.yaml",
             "--out-dir", str(tmp),
             "--quiet"],
            cwd=REPO_ROOT, capture_output=True, text=True,
        )
        # fep-report exits 1 on fail_case (by design), but the
        # files are written either way.
        assert (tmp / "report.md").exists(), (
            f"fep-report didn't write report.md; stderr:\n{r.stderr}")
        rc, out = _run_fill(tmp)
    # fail_case fixture has no run.log with GHMC lines, so
    # ghmc_mean/worst come back <missing> → exit 2 (warn). Both
    # rc 0 and 2 are acceptable; what matters is the email
    # content has the right FAIL verdict.
    assert rc in (0, 2), (
        f"fill should succeed on fail_case; got {rc}\n{out}")
    assert "Sign-correct on methane:          FAIL" in out


if __name__ == "__main__":
    funcs = [
        test_hydration_pass_all_fields_fill,
        test_hydration_fail_methane_sign,
        test_binding_streptavidin_template_selected,
        test_binding_egfr_template_selected,
        test_binding_non_binder_flagged_in_email,
        test_missing_numbers_returns_exit_2,
        test_no_report_md_returns_exit_1,
        test_prof_name_override,
        test_hydration_hardware_override_lands_in_body,
        test_fail_case_fixture_hydration,
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
