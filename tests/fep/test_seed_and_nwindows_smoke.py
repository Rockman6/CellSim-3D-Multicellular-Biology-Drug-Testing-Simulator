"""Layer 1.3 FEP — seed plumbing (#2) + n_windows guard (#8).

Two BUG_AUDIT.md regressions, both testable without running MD:

- #2: the `seed` parameter never reached any integrator. `_seed_ghmc_
  move` now sets it on the integrator the GHMCMove builds. But seeding
  alone does not give reproducibility — the real cause is platform
  non-determinism (and openmmtools GHMC is seed-independent on the
  deterministic Reference platform). Reproducibility is delivered by
  `sample_alchemical_windows(..., deterministic=True)`, which pins the
  Reference platform on every context and is bitwise-reproducible
  (see test_deterministic_mode_smoke). These seed tests assert only
  that the seed lands on the built integrator.
- #8: `n_windows < 3` silently clamped inside `_split_lambda_schedule`
  while the caller kept the original `n_windows`, desyncing the u_kn
  matrix and crashing with an out-of-bounds write. We verify the
  invalid value is now rejected with a clear error at every entry
  point.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


# --------------------------------------------------------------------
# #2 — seed plumbing
# --------------------------------------------------------------------

def _ghmc_move(n_steps=5):
    from openmmtools import mcmc
    from openmm import unit
    return mcmc.GHMCMove(
        timestep=1.0 * unit.femtosecond,
        collision_rate=20 / unit.picosecond,
        n_steps=n_steps)


def _thermo_state():
    from openmmtools import states, testsystems
    from openmm import unit
    sysobj = testsystems.HarmonicOscillator()
    return states.ThermodynamicState(
        system=sysobj.system, temperature=298.15 * unit.kelvin)


def test_seed_reaches_ghmc_integrator():
    """The integrator built by a seeded move must carry the requested
    RNG seed. Before the fix, setRandomNumberSeed was never called."""
    from src.fep.sampling import _seed_ghmc_move

    move = _seed_ghmc_move(_ghmc_move(), 4242)
    integ = move._get_integrator(_thermo_state())
    assert integ.getRandomNumberSeed() == 4242


def test_seed_zero_maps_to_nonzero():
    """OpenMM treats seed 0 as 'nondeterministic'; the helper must map
    it to a fixed positive value so seed=0 is still reproducible."""
    from src.fep.sampling import _seed_ghmc_move

    move = _seed_ghmc_move(_ghmc_move(), 0)
    integ = move._get_integrator(_thermo_state())
    assert integ.getRandomNumberSeed() != 0


def test_seed_helper_updatable_not_double_wrapped():
    """Re-seeding an already-hooked move updates the seed in place
    rather than stacking wrappers."""
    from src.fep.sampling import _seed_ghmc_move

    move = _ghmc_move()
    _seed_ghmc_move(move, 1)
    hooked = move._get_integrator
    _seed_ghmc_move(move, 99)
    assert move._get_integrator is hooked, "wrapper was re-installed"
    integ = move._get_integrator(_thermo_state())
    assert integ.getRandomNumberSeed() == 99


# --------------------------------------------------------------------
# #8 — n_windows guard
# --------------------------------------------------------------------

@pytest.mark.parametrize("bad", [0, 1, 2])
def test_split_lambda_schedule_rejects_small_nwindows(bad):
    from src.fep.sampling import _split_lambda_schedule
    with pytest.raises(ValueError, match="n_windows must be >= 3"):
        _split_lambda_schedule(bad)


def test_split_lambda_schedule_ok_at_three():
    from src.fep.sampling import _split_lambda_schedule
    sched = _split_lambda_schedule(3)
    assert len(sched) == 3
    assert sched[0] == (0.0, 0.0) and sched[-1] == (1.0, 1.0)


@pytest.mark.parametrize("bad", [1, 2])
def test_sampler_rejects_small_nwindows_before_any_build(bad):
    """The guard fires before touching the (dummy) system, so we can
    pass None for the system/positions."""
    from src.fep.sampling import sample_alchemical_windows
    with pytest.raises(ValueError, match="n_windows must be >= 3"):
        sample_alchemical_windows(None, None, None, n_windows=bad)


def test_hydration_entry_returns_clean_error_on_small_nwindows():
    """Public API returns a non-ok envelope with a clear reason — no
    crash, no expensive OpenFF build."""
    from src.fep import compute_hydration_dg
    r = compute_hydration_dg("C", n_windows=2)
    assert not r.ok
    assert "n_windows must be >= 3" in r.reason


def test_binding_entry_returns_clean_error_on_small_nwindows():
    from src.fep.binding import compute_absolute_binding_dg
    # Bogus receptor path is never read — the n_windows guard returns
    # first.
    r = compute_absolute_binding_dg(
        "C", "/nonexistent/receptor.pdb", n_windows=1, sample=True)
    assert not r.ok
    assert "n_windows must be >= 3" in r.reason


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
