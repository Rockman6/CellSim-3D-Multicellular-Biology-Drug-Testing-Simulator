"""Full-system transient — anchored to the unified steady state.

The decisive check: the full composed ODE (permeation + pH + efflux +
sink) must converge to `solve_steady_state`. Also pins that the binding
sink acts as a capacitor (slows the approach, same steady state) and that
occupancy tracks the free concentration.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    Permeability, spherical_cell_geometry,
    EffluxPump, BindingSink, COMPARTMENT_PH,
    solve_steady_state, simulate_disposition,
)

_PERM = Permeability(P_cm_per_s=1e-6)
_GEOM = spherical_cell_geometry(10.0)


def test_ode_converges_to_unified_steady_state_full_stack():
    """THE anchor: transient t→∞ == solve_steady_state, all effects on."""
    pump = EffluxPump(Vmax_M_per_s=2e-9, Km_M=1e-7)
    sink = BindingSink(capacity_M=1e-4, Kd_M=1e-6)
    kw = dict(pKa=9.0, ion_type="base", pH_in=COMPARTMENT_PH["lysosome"],
              pump=pump, sink=sink)
    tc = simulate_disposition(1e-7, permeability=_PERM, geometry=_GEOM, **kw)
    ss = solve_steady_state(1e-7, permeability=_PERM, geometry=_GEOM, **kw)
    assert ss.free_in_M > 0
    assert abs(tc.converged_free_M - ss.free_in_M) / ss.free_in_M < 1e-3
    assert abs(tc.steady_state_free_M - ss.free_in_M) / ss.free_in_M < 1e-12


def test_free_starts_at_zero_and_is_monotone_up():
    tc = simulate_disposition(1e-7, permeability=_PERM, geometry=_GEOM)
    assert tc.free_M[0] < 1e-12
    for a, b in zip(tc.free_M, tc.free_M[1:]):
        assert b >= a - 1e-12


def test_binding_sink_slows_approach_same_steady_state():
    """Capacitor behaviour: bigger sink → longer t½, identical steady state."""
    no_sink = simulate_disposition(1e-7, permeability=_PERM, geometry=_GEOM)
    big_sink = simulate_disposition(
        1e-7, permeability=_PERM, geometry=_GEOM,
        sink=BindingSink(capacity_M=1e-3, Kd_M=1e-7))
    # Same steady-state free (sink doesn't change it)...
    assert abs(big_sink.steady_state_free_M - no_sink.steady_state_free_M) \
        / no_sink.steady_state_free_M < 1e-9
    # ...but a much longer time to get halfway there.
    assert big_sink.time_to_fraction_of_steady_state(0.5) > \
        5 * no_sink.time_to_fraction_of_steady_state(0.5)


def test_occupancy_trace_tracks_free():
    prior = binding_to_hill(-9.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    tc = simulate_disposition(1e-7, permeability=_PERM, geometry=_GEOM,
                              prior=prior)
    assert tc.occupancy is not None
    assert tc.occupancy[0] < 1e-6           # nothing bound at t=0
    assert tc.occupancy[-1] > tc.occupancy[0]


def test_no_sink_matches_total_equals_free():
    tc = simulate_disposition(1e-7, permeability=_PERM, geometry=_GEOM)
    for f, t in zip(tc.free_M, tc.total_M):
        assert abs(f - t) < 1e-18
