"""Regression test for the opt-in Hamiltonian replica exchange
path in sample_alchemical_windows.

Why: the hand-rolled independent-replica path fails MBAR overlap
on tight intramolecular charge networks (acetic_acid, acetamide,
biotin/streptavidin). The replica-exchange path (opt-in via
use_replica_exchange=True) is the fix. Without a regression test
pinning that path, a future refactor could silently break it and
we wouldn't notice until the next overnight bench failed."""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


def test_replica_exchange_flag_wired_through_compute_hydration_dg():
    """compute_hydration_dg must accept use_replica_exchange and
    pass it down. Code-level inspection, no MD."""
    import inspect
    from src.fep import compute_hydration_dg
    sig = inspect.signature(compute_hydration_dg)
    assert "use_replica_exchange" in sig.parameters
    assert sig.parameters["use_replica_exchange"].default is False


def test_replica_exchange_flag_wired_through_compute_binding_dg():
    """compute_absolute_binding_dg must accept use_replica_exchange."""
    import inspect
    from src.fep.binding import compute_absolute_binding_dg
    sig = inspect.signature(compute_absolute_binding_dg)
    assert "use_replica_exchange" in sig.parameters
    assert sig.parameters["use_replica_exchange"].default is False


def test_sample_alchemical_windows_accepts_replica_exchange_flag():
    import inspect
    from src.fep.sampling import sample_alchemical_windows
    sig = inspect.signature(sample_alchemical_windows)
    assert "use_replica_exchange" in sig.parameters
    assert sig.parameters["use_replica_exchange"].default is False


def test_replica_exchange_runs_methane_vacuum_end_to_end():
    """Run sample_alchemical_windows with use_replica_exchange=True
    on methane vacuum. Tiny sampling (3 windows × 5 iter × 50 steps,
    ~30 s CPU). Just verifies the pipeline runs end-to-end on the
    opt-in path; numerical accuracy is irrelevant at this budget."""
    import glob
    import os
    import tempfile

    from src.fep import _build_alchemical_legs
    from src.fep.sampling import sample_alchemical_windows

    (vac_alch, _solv_alch, vac_top, _solv_top,
     vac_pos, _solv_pos, _n) = _build_alchemical_legs("C")

    # BUG_AUDIT.md #13: the RE path used to leak a tmpdir + .nc per
    # call. Snapshot the scratch dirs before/after and assert no leak.
    rex_glob = os.path.join(tempfile.gettempdir(), "cellsim_rex_*")
    before = set(glob.glob(rex_glob))

    r = sample_alchemical_windows(
        vac_alch, vac_top, vac_pos,
        n_windows=3,
        n_equilibration_steps=20,
        n_production_steps=250,
        sample_stride=50,
        seed=1,
        use_replica_exchange=True,
    )
    print(f"  {r.summary()}")

    leaked = set(glob.glob(rex_glob)) - before
    assert not leaked, (
        f"replica-exchange leaked scratch dirs (BUG_AUDIT #13): "
        f"{sorted(leaked)}")

    assert r.ok, (
        f"replica-exchange path must run end-to-end on methane "
        f"vacuum; got: {r.reason}")
    assert r.dG_kcalmol is not None
    assert math.isfinite(r.dG_kcalmol)
    assert math.isfinite(r.dG_uncertainty_kcalmol)
    # Methane vacuum (no intramolecular nonbondeds) → dG ≈ 0 at
    # ANY sampling budget. Generous bound for noise.
    assert abs(r.dG_kcalmol) < 5.0, (
        f"methane vacuum dG should be ~0 (no intra nonbondeds); "
        f"got {r.dG_kcalmol:+.3f}")


def test_cli_bench_has_replica_exchange_flag():
    """`cellsim fep-binding bench --help` must expose
    --replica-exchange so biologists can enable the fix without
    editing code."""
    import subprocess
    r = subprocess.run(
        ["python", "-m", "src.fep.binding", "bench", "--help"],
        cwd=REPO_ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "--replica-exchange" in r.stdout, (
        f"--replica-exchange missing from bench --help; got:\n"
        f"{r.stdout[:1000]}")


if __name__ == "__main__":
    funcs = [
        test_replica_exchange_flag_wired_through_compute_hydration_dg,
        test_replica_exchange_flag_wired_through_compute_binding_dg,
        test_sample_alchemical_windows_accepts_replica_exchange_flag,
        test_cli_bench_has_replica_exchange_flag,
        test_replica_exchange_runs_methane_vacuum_end_to_end,
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
