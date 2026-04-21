"""fep-report analyser regression tests.

Pins the pass/fail/inconclusive verdict logic on three canonical
inputs that live in `tests/fep/fixtures/`:

  ok_case/    12-row FreeSolv CSV with MAE 0.42, Pearson 0.99,
              3 compounds reporting GHMC acceptance → PASS.
  fail_case/  12-row FreeSolv CSV with methane sign flipped
              (−1.85 vs expt +2.00), MAE 1.67 → FAIL.
  scaffold    (synthesised inline) 2-row CSV with all predictions
              blank → inconclusive (exit 0).

The analyser is the Milestone-A-to-Campaign-2 gate — the instant
the M5 Max tarball lands, biologists run one command and get a
pass/fail verdict. This test file protects that end-to-end path
from silent drift: if the MAE formula, sign check, or tri-state
logic regresses, CI will catch it.

Fixtures chosen over live-running fep-freesolv because:
1. This has to run on laptops without AM1-BCC / OpenMM in < 1 s,
   not a 2-hour FEP gate.
2. We want to test the analyser, not the sampler. Testing the
   analyser on real data means the test can't distinguish an
   analyser bug from a sampler bug.
"""

from __future__ import annotations

import csv
import io
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.fep.report import (  # noqa: E402
    analyse, format_markdown, _write_table_csv,
    GATE_MAE_KCALMOL,
)

FIXTURES = Path(__file__).parent / "fixtures"


def test_ok_case_passes_gate():
    r = analyse(FIXTURES / "ok_case")
    assert r.pass_overall is True, (
        f"ok_case must PASS; pass_mae={r.pass_mae} "
        f"pass_ghmc={r.pass_ghmc} "
        f"pass_sign_critical={r.pass_sign_critical}")
    assert r.pass_mae is True
    assert r.mae_kcalmol is not None
    assert r.mae_kcalmol < GATE_MAE_KCALMOL
    # GHMC stats should be parsed from run.log (3 compounds
    # emit GHMC lines in the fixture).
    assert r.ghmc_accept_compound_count == 3, (
        f"expected 3 GHMC-reporting compounds; got "
        f"{r.ghmc_accept_compound_count}")
    assert r.pass_ghmc is True
    # Sign check: methane (+1.72) and acetamide (−8.90) signs
    # match expt (+2.00 / −9.71).
    assert r.pass_sign_critical is True
    # Pearson r on 12 compounds with MAE 0.42 should be > 0.95.
    assert r.pearson_r is not None and r.pearson_r > 0.95


def test_fail_case_fails_sign_critical():
    r = analyse(FIXTURES / "fail_case")
    # MAE 1.67 fails the default 1.5 gate.
    assert r.pass_mae is False
    # Methane predicted −1.85 vs expt +2.00 → sign WRONG.
    assert r.pass_sign_critical is False, (
        f"fail_case has methane sign flipped; should FAIL sign "
        f"check. got {r.pass_sign_critical}")
    # Overall verdict: FAIL.
    assert r.pass_overall is False


def test_fail_case_mae_gate_override_still_fails_on_sign():
    """With --mae-gate 2.0 the MAE check passes (1.67 < 2.0), but
    the sign-critical failure must still make overall FAIL. This
    is the "don't let permissive thresholds mask physics errors"
    guard."""
    r = analyse(FIXTURES / "fail_case", mae_gate_kcalmol=2.0)
    assert r.pass_mae is True
    assert r.pass_sign_critical is False
    assert r.pass_overall is False


def test_scaffold_only_csv_is_inconclusive_not_fail():
    """Generate a 2-row CSV with ok=True but no predictions (the
    `fep-binding bench` scaffold-only output). Should verdict as
    'inconclusive' and exit 0 — CI must not fail pre-sample runs."""
    scaffold_rows = [
        {"name": "methane", "smiles": "C",
         "dG_expt_kcalmol": -2.0, "dG_pred_kcalmol": "",
         "uncertainty_kcalmol": "", "residual_kcalmol": "",
         "ghmc_accept_mean": "", "ghmc_accept_min": "",
         "wall_seconds": 16.4, "ok": True, "reason": ""},
        {"name": "ethane", "smiles": "CC",
         "dG_expt_kcalmol": -1.5, "dG_pred_kcalmol": "",
         "uncertainty_kcalmol": "", "residual_kcalmol": "",
         "ghmc_accept_mean": "", "ghmc_accept_min": "",
         "wall_seconds": 12.1, "ok": True, "reason": ""},
    ]
    with tempfile.TemporaryDirectory(
            prefix="cellsim_report_test_") as tmp:
        csv_path = Path(tmp) / "results.csv"
        with csv_path.open("w", newline="",
                           encoding="utf-8-sig") as fo:
            w = csv.DictWriter(
                fo, fieldnames=list(scaffold_rows[0].keys()))
            w.writeheader()
            for row in scaffold_rows:
                w.writerow(row)
        r = analyse(csv_path)

    assert r.n_ok == 2, (
        f"n_ok should count scaffold-built rows; got {r.n_ok}")
    assert r.pass_mae is None
    assert r.pass_sign_critical is None
    assert r.pass_overall is None, (
        "all-inconclusive should give None overall, not True/False")


def test_markdown_render_stable():
    """Regression guard on the markdown format. If someone edits
    format_markdown and the header / verdict lines change shape,
    CI will fail — protects the format downstream tooling parses."""
    r = analyse(FIXTURES / "ok_case")
    md = format_markdown(r)
    # Verdict header is the first line and must say PASS.
    assert md.startswith("# FreeSolv FEP report — PASS"), (
        f"header line changed:\n{md[:200]}")
    # Gate-verdict section must reference the 1.5 gate.
    assert f"≤ {GATE_MAE_KCALMOL} kcal/mol" in md
    # Per-compound table header must be present and stable.
    assert "| name | smiles | expt | pred | residual |" in md
    # Sign-critical compounds must be labelled explicitly.
    assert "acetamide/methane" in md


def test_tarball_path_unwraps_cleanly():
    """`cellsim fep-report freesolv_m5max_<stamp>.tar.gz` is the
    friend-workflow use case. Verify we can tar the ok_case and
    get the same verdict."""
    import tarfile
    with tempfile.TemporaryDirectory(
            prefix="cellsim_report_tar_") as tmp:
        tb = Path(tmp) / "ok_case.tar.gz"
        with tarfile.open(tb, "w:gz") as tf:
            tf.add(FIXTURES / "ok_case", arcname="run")
        r = analyse(tb)
    assert r.pass_overall is True
    assert r.n_ok == 12


def test_partial_run_flagged_even_if_per_row_looks_ok():
    """Friend's M5 Max dies at 6/12. Analyser must refuse to PASS
    even if those 6 rows look fine per-row, because the other 6
    might be the hard ones (methane / acetamide included). Matches
    the professor's 'do not pass on incomplete data' rule."""
    r = analyse(FIXTURES / "ok_case", expected_rows=24)
    assert r.is_partial is True
    assert r.pass_overall is False, (
        "partial run must fail overall even if MAE on rows "
        "present looks good")
    # n_ok still reports the actual rows present — don't hide the
    # count, just make the verdict partial.
    assert r.n_ok == 12
    assert r.expected_rows == 24


def test_full_run_matches_expected_passes():
    """When n_total == expected_rows, no partial flag; PASS if
    gates pass."""
    r = analyse(FIXTURES / "ok_case", expected_rows=12)
    assert r.is_partial is False
    assert r.pass_overall is True


def test_partial_run_markdown_has_partial_header():
    """The verdict header must scream 'partial run' so a skimming
    biologist doesn't mistakenly treat the report as a pass."""
    r = analyse(FIXTURES / "ok_case", expected_rows=24)
    md = format_markdown(r)
    assert "partial" in md.lower(), (
        f"header should contain 'partial' for n_total < expected; "
        f"got: {md[:200]}")
    assert "PARTIAL — Milestone A cannot pass" in md


if __name__ == "__main__":
    funcs = [
        test_ok_case_passes_gate,
        test_fail_case_fails_sign_critical,
        test_fail_case_mae_gate_override_still_fails_on_sign,
        test_scaffold_only_csv_is_inconclusive_not_fail,
        test_markdown_render_stable,
        test_tarball_path_unwraps_cleanly,
        test_partial_run_flagged_even_if_per_row_looks_ok,
        test_full_run_matches_expected_passes,
        test_partial_run_markdown_has_partial_header,
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
