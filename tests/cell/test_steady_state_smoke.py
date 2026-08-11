"""Unified steady state — and the cross-checks that it composes correctly.

The composed solver is trustworthy only if it reduces to each standalone
module in its limit. These tests pin exactly that:
  * no pump, no pH gradient        → free = C_out;
  * no pump (pH only)              → free = C_out · R (partitioning module);
  * no pH gradient, no sink        → free matches efflux_steady_state;
  * the sink changes TOTAL but NOT the steady-state free (the physics
    insight the composition makes explicit);
  * occupancy is taken at the free concentration.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    Permeability, spherical_cell_geometry,
    EffluxPump, efflux_steady_state,
    BindingSink,
    partition_across_membrane, COMPARTMENT_PH,
    occupancy_from_prior,
    solve_steady_state,
)

_PERM = Permeability(P_cm_per_s=1e-6)
_GEOM = spherical_cell_geometry(10.0)


def test_reduces_to_no_gradient_without_pump_or_pH():
    ss = solve_steady_state(1e-7, permeability=_PERM, geometry=_GEOM)
    assert abs(ss.free_in_M - 1e-7) / 1e-7 < 1e-9
    assert abs(ss.free_accumulation - 1.0) < 1e-9


def test_reduces_to_partitioning_module_without_pump():
    """pH trapping only → free = C_out · R from the partitioning module."""
    C_out = 1e-7
    ss = solve_steady_state(
        C_out, permeability=_PERM, geometry=_GEOM,
        pKa=9.0, ion_type="base", pH_in=COMPARTMENT_PH["lysosome"],
        pH_out=COMPARTMENT_PH["blood"])
    part = partition_across_membrane(
        C_out, pKa=9.0, pH_in=COMPARTMENT_PH["lysosome"],
        pH_out=COMPARTMENT_PH["blood"], ion_type="base")
    assert abs(ss.free_in_M - part.total_in_M) / part.total_in_M < 1e-9


def test_reduces_to_efflux_module_without_pH_or_sink():
    """Pump only → free matches the standalone efflux steady state."""
    pump = EffluxPump(Vmax_M_per_s=5e-9, Km_M=1e-7)
    C_out = 1e-7
    ss = solve_steady_state(C_out, permeability=_PERM, geometry=_GEOM,
                            pump=pump)
    ref = efflux_steady_state(C_out, _PERM, _GEOM, pump)
    assert abs(ss.free_in_M - ref.C_in_M) / ref.C_in_M < 1e-9
    assert ss.residual < 1e-18


def test_sink_changes_total_not_steady_state_free():
    """The binding sink stores drug (raises total) but leaves free fixed."""
    pump = EffluxPump(Vmax_M_per_s=5e-9, Km_M=1e-7)
    C_out = 1e-7
    no_sink = solve_steady_state(C_out, permeability=_PERM, geometry=_GEOM,
                                 pump=pump)
    with_sink = solve_steady_state(
        C_out, permeability=_PERM, geometry=_GEOM, pump=pump,
        sink=BindingSink(capacity_M=1e-4, Kd_M=1e-7))
    # Free unchanged...
    assert abs(with_sink.free_in_M - no_sink.free_in_M) / no_sink.free_in_M < 1e-9
    # ...but total is much larger (drug stored on the sink).
    assert with_sink.total_in_M > 5 * no_sink.total_in_M
    assert with_sink.total_accumulation > with_sink.free_accumulation


def test_occupancy_uses_free_concentration():
    prior = binding_to_hill(-9.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb", method="vina")
    pump = EffluxPump(Vmax_M_per_s=5e-9, Km_M=1e-7)
    ss = solve_steady_state(1e-7, permeability=_PERM, geometry=_GEOM, pump=pump)
    expected = occupancy_from_prior(prior, ss.free_in_M).theta
    assert abs(ss.occupancy(prior) - expected) < 1e-12


def test_trust_is_limiting_input():
    """Any uncalibrated physics input makes the composed state uncalibrated."""
    ss = solve_steady_state(
        1e-7, permeability=Permeability(P_cm_per_s=1e-6, trust="calibrated"),
        geometry=_GEOM, pump=EffluxPump(Vmax_M_per_s=1e-9, Km_M=1e-7,
                                        trust="uncalibrated"))
    assert not ss.inputs_calibrated
    assert ss.limiting_trust == "uncalibrated"
