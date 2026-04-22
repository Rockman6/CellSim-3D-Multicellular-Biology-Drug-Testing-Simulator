"""Regression tests for fep-report provenance extraction.

The shell wrapper scripts evolved over time in where they tee the
`git commit:` line — earlier versions sent it to stdout only, then
into run.log, then into env.log. `cellsim fep-report` has to find
it wherever it landed, because a biologist handed a tarball from
someone else's machine cannot re-run git commands: the commit hash
IS the provenance.

Pins: _parse_run_log scans run.log, env.log, and doctor.log and
returns the first git_commit / platform it finds.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.fep.report import (
    _parse_run_log,
    _scan_sibling_logs_for_provenance,
)


def test_git_commit_picked_up_from_run_log():
    with tempfile.TemporaryDirectory(prefix="fep_prov_") as tmp:
        d = Path(tmp)
        (d / "run.log").write_text(
            "CellSim\n"
            "  git commit:  abc1234567890abcdef\n"
            "FEP sampling platform: CUDA\n"
            "[freesolv] methane C\n")
        meta = _parse_run_log(d / "run.log", rows=[])
    assert meta["git_commit"] == "abc1234567890abcdef", meta
    assert meta["platform"] == "CUDA", meta


def test_git_commit_fallback_to_env_log():
    """Old scripts sent the header to stdout and env.log. run.log
    contains only bench output. fep-report must still populate
    git_commit by falling back to env.log."""
    with tempfile.TemporaryDirectory(prefix="fep_prov_") as tmp:
        d = Path(tmp)
        (d / "env.log").write_text(
            "CellSim — FreeSolv FEP gate\n"
            "  git commit:  deadbeefcafe0000111122223333\n"
            "  openmm 8.2.0\n")
        (d / "run.log").write_text(
            "[freesolv] methane C\n"
            "...bench output without provenance...\n"
            "FEP sampling platform: Metal\n")
        meta = _parse_run_log(d / "run.log", rows=[])
    assert meta["git_commit"] == "deadbeefcafe0000111122223333", meta
    assert meta["platform"] == "Metal", meta


def test_git_commit_fallback_when_no_run_log():
    """Malformed tarball (no run.log) — only env.log has provenance.
    Parser must still recover git_commit via sibling scan."""
    with tempfile.TemporaryDirectory(prefix="fep_prov_") as tmp:
        d = Path(tmp)
        (d / "env.log").write_text(
            "  git commit:  1111222233334444aaaabbbbccccdddd\n")
        # Intentionally no run.log.
        meta = _parse_run_log(d / "run.log", rows=[])
    assert meta["git_commit"] == "1111222233334444aaaabbbbccccdddd", meta


def test_doctor_log_as_last_resort():
    """Some runs only landed provenance in doctor.log."""
    with tempfile.TemporaryDirectory(prefix="fep_prov_") as tmp:
        d = Path(tmp)
        (d / "doctor.log").write_text(
            "cellsim doctor\n"
            "git commit: feedface99887766\n"
            "OK\n")
        meta = _parse_run_log(d / "run.log", rows=[])
    assert meta["git_commit"] == "feedface99887766", meta


def test_no_logs_at_all_returns_none():
    """Nothing to find — git_commit stays None, doesn't crash."""
    with tempfile.TemporaryDirectory(prefix="fep_prov_") as tmp:
        d = Path(tmp)
        meta = _parse_run_log(d / "run.log", rows=[])
    assert meta["git_commit"] is None
    assert meta["platform"] is None


def test_run_log_wins_over_sibling():
    """When both run.log and env.log have git commit, run.log wins
    (most local to the bench process)."""
    with tempfile.TemporaryDirectory(prefix="fep_prov_") as tmp:
        d = Path(tmp)
        (d / "run.log").write_text(
            "git commit: 0000000aaaaaaaabbbbbbbbbbcccccc\n")
        (d / "env.log").write_text(
            "git commit: ffffffffffffffffffffffffffffffff\n")
        meta = _parse_run_log(d / "run.log", rows=[])
    assert meta["git_commit"] == "0000000aaaaaaaabbbbbbbbbbcccccc", meta


def test_sibling_scan_direct_helper():
    """_scan_sibling_logs_for_provenance is a public-ish helper;
    cover it directly so the fallback policy (env > doctor > header)
    is pinned."""
    with tempfile.TemporaryDirectory(prefix="fep_prov_") as tmp:
        d = Path(tmp)
        (d / "env.log").write_text("git commit: aaaaaaaaaaaaaaaa\n")
        (d / "doctor.log").write_text("git commit: bbbbbbbbbbbbbbbb\n")
        meta = _scan_sibling_logs_for_provenance(
            d, {"git_commit": None, "platform": None,
                "ghmc_means": [], "ghmc_mins": []})
    # env.log scanned first.
    assert meta["git_commit"] == "aaaaaaaaaaaaaaaa", meta


if __name__ == "__main__":
    funcs = [
        test_git_commit_picked_up_from_run_log,
        test_git_commit_fallback_to_env_log,
        test_git_commit_fallback_when_no_run_log,
        test_doctor_log_as_last_resort,
        test_no_logs_at_all_returns_none,
        test_run_log_wins_over_sibling,
        test_sibling_scan_direct_helper,
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
