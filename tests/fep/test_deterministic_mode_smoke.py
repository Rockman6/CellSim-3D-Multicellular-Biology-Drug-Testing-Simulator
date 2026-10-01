"""Layer 1.3 FEP — deterministic mode reproducibility (BUG_AUDIT.md #2).

The #2 seed fix alone does not make runs reproducible (seeding reaches
the integrator, but the fast platforms are non-deterministic and the
openmmtools GHMC integrator is seed-independent on the deterministic
Reference platform). `deterministic=True` delivers real reproducibility
by pinning the Reference platform on EVERY context (per-window minimise,
GHMC propagation, cross-λ energy eval) and using a fresh, un-shared
context cache so no state leaks between runs.

This test is the reproducibility witness. It runs on the Reference
platform (single-threaded, slow) so it deliberately uses tiny sampling.
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


def test_deterministic_flag_is_plumbed_through_entry_points():
    from src.fep import compute_hydration_dg
    from src.fep.binding import compute_absolute_binding_dg
    from src.fep.sampling import (
        sample_alchemical_windows, sample_restraint_coupling_dg)
    for fn in (sample_alchemical_windows, sample_restraint_coupling_dg,
               compute_hydration_dg, compute_absolute_binding_dg):
        p = inspect.signature(fn).parameters.get("deterministic")
        assert p is not None and p.default is False, fn.__name__


def test_deterministic_mode_is_bitwise_reproducible():
    """Two deterministic runs of the same input must give bitwise-equal
    ΔG. Before deterministic mode, two seed=1 runs differed by ~1
    kcal/mol on every platform."""
    from src.fep import _build_alchemical_legs
    from src.fep.sampling import sample_alchemical_windows

    # Use the VACUUM leg of ethanol: 9 atoms, no water box, so the
    # single-threaded Reference platform stays fast (~1 min), yet it has
    # real intramolecular LJ/charge so the decoupling ΔG is NONZERO —
    # a meaningful reproducibility witness (methane vacuum would be
    # identically 0 and prove nothing). The same machinery (minimise +
    # GHMC + eval contexts) is exercised as for a solvated system.
    vac, _solv, vtop, _stop, vpos, _spos, _n = _build_alchemical_legs("CCO")
    kw = dict(n_windows=3, n_equilibration_steps=30,
              n_production_steps=100, sample_stride=50, seed=1)
    r1 = sample_alchemical_windows(vac, vtop, vpos,
                                   deterministic=True, **kw)
    r2 = sample_alchemical_windows(vac, vtop, vpos,
                                   deterministic=True, **kw)
    assert r1.ok, r1.reason
    assert r2.ok, r2.reason
    assert r1.dG_kcalmol is not None
    assert r1.dG_kcalmol == r2.dG_kcalmol, (
        f"deterministic mode not reproducible: "
        f"{r1.dG_kcalmol!r} vs {r2.dG_kcalmol!r}")


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
