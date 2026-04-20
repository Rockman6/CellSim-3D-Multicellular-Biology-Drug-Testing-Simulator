"""Layer 1.3 kinase-receptor heads-up warning smoke.

Verifies:
  - bundled kinase PDB (1m17) triggers the warning
  - serine-protease PDB (3ptb) does NOT trigger it
  - biotin-binder PDB (1stp) does NOT trigger it
  - a non-existent receptor path is silently skipped (no crash)

The warning is print-to-stdout; we capture stdout and check for
the expected sentinel phrase.
"""

from __future__ import annotations

import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.batch import _warn_kinase_receptor  # noqa: E402


def _fires(pdb: Path) -> bool:
    buf = io.StringIO()
    with redirect_stdout(buf):
        _warn_kinase_receptor(pdb)
    return "kinase" in buf.getvalue().lower()


def test_kinase_warning_fires_on_bundled_kinase_pdb():
    assert _fires(REPO_ROOT / "benchmarks" / "dock" / "1m17.pdb")


def test_kinase_warning_silent_on_protease():
    assert not _fires(REPO_ROOT / "benchmarks" / "dock" / "3ptb.pdb")


def test_kinase_warning_silent_on_biotin_binder():
    assert not _fires(REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb")


def test_kinase_warning_silent_on_cyp_heme():
    # 1tqn is a CYP3A4 crystal with heme — not a kinase.
    assert not _fires(REPO_ROOT / "benchmarks" / "dock" / "1tqn.pdb")


def test_kinase_warning_silent_on_missing_file():
    # Should not raise.
    _warn_kinase_receptor(
        Path("/does/not/exist/fake_receptor.pdb"))


if __name__ == "__main__":
    tests = [
        test_kinase_warning_fires_on_bundled_kinase_pdb,
        test_kinase_warning_silent_on_protease,
        test_kinase_warning_silent_on_biotin_binder,
        test_kinase_warning_silent_on_cyp_heme,
        test_kinase_warning_silent_on_missing_file,
    ]
    for t in tests:
        try:
            t()
            print(f"  {t.__name__} PASS")
        except AssertionError as e:
            print(f"  {t.__name__} FAIL: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"  {t.__name__} ERROR: {e}")
            sys.exit(2)
    print("PASS")
