"""Layer 1.3 --receptor PDB-ID auto-fetch resolver smoke.

Covers non-network behaviour:
  - existing file path returns unchanged
  - non-PDB-ID string (contains '/', too long, has punctuation)
    returns unchanged (no network hit)
  - a 4-char alphanumeric token that ISN'T an existing path is
    the only trigger for a fetch attempt (verified indirectly by
    patching urlopen)

We don't exercise the actual RCSB fetch here — CI should not
depend on the network. A user-side smoke that the path resolver
wiring works is enough; the fetch itself is simple urllib and
visible in the diff.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.receptor_resolve import resolve_receptor as _resolve_receptor  # noqa: E402


def test_existing_path_passthrough():
    p = "benchmarks/dock/1stp.pdb"
    out = _resolve_receptor(p)
    assert out == p, out


def test_non_pdb_id_passthrough():
    # Contains non-alnum → can't be a PDB ID.
    p = "nonexistent-path/receptor.pdb"
    out = _resolve_receptor(p)
    assert out == p


def test_longer_than_4_chars_passthrough():
    p = "abcdef"          # alnum but wrong length
    out = _resolve_receptor(p)
    assert out == p


def test_exactly_3_chars_passthrough():
    p = "abc"             # too short
    out = _resolve_receptor(p)
    assert out == p


def test_already_downloaded_no_refetch(monkeypatch=None):
    # If data/receptors/1stp.pdb exists (as it will after a
    # fresh fetch) the resolver should return the cached path
    # without hitting the network. We simulate this by writing
    # a dummy file and patching urlopen to raise.
    import tempfile
    from unittest import mock
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        cache = td / "data" / "receptors"
        cache.mkdir(parents=True)
        (cache / "xxxx.pdb").write_text("HEADER   dummy")
        # Run with cwd = td so data/receptors/xxxx.pdb exists.
        import os
        old_cwd = os.getcwd()
        os.chdir(td)
        try:
            with mock.patch("urllib.request.urlopen",
                             side_effect=RuntimeError("network")):
                out = _resolve_receptor("xxxx")
            # Resolver returns relative 'data/receptors/xxxx.pdb'
            # from the td cwd; resolve both sides before compare.
            assert Path(out).resolve() == (cache / "xxxx.pdb").resolve(), out
            print(f"  cached PDB reused: {out}")
        finally:
            os.chdir(old_cwd)


if __name__ == "__main__":
    tests = [
        test_existing_path_passthrough,
        test_non_pdb_id_passthrough,
        test_longer_than_4_chars_passthrough,
        test_exactly_3_chars_passthrough,
        test_already_downloaded_no_refetch,
    ]
    for t in tests:
        try:
            t()
            print(f"  {t.__name__} PASS")
        except AssertionError as e:
            print(f"  {t.__name__} FAIL: {e}")
            sys.exit(1)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  {t.__name__} ERROR: {e}")
            sys.exit(2)
    print("PASS")
