"""Active efflux (P-gp): the pump holds the drug out — until it saturates.

Pins the steady-state influx=efflux balance (residual ~0), the limits
(no pump → C_in=C_out; strong pump → C_in≪C_out), saturation recovery at
high dose, and construction from a Michaelis bridge prior.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    Permeability,
    spherical_cell_geometry,
    EffluxPump,
    efflux_steady_state,
    pump_from_michaelis_prior,
)

_PERM = Permeability(P_cm_per_s=1e-6)
_GEOM = spherical_cell_geometry(10.0)


def _ss(C_out, Vmax, Km):
    pump = EffluxPump(Vmax_M_per_s=Vmax, Km_M=Km)
    return efflux_steady_state(C_out, _PERM, _GEOM, pump)


def test_no_pump_gives_no_gradient():
    ss = _ss(1e-7, Vmax=0.0, Km=1e-6)
    assert abs(ss.C_in_M - ss.C_out_M) < 1e-18
    assert abs(ss.accumulation_ratio - 1.0) < 1e-12


def test_steady_state_balances_influx_and_efflux():
    ss = _ss(1e-7, Vmax=1e-9, Km=1e-7)
    # The solved C_in must make passive influx equal pump efflux.
    assert ss.residual < 1e-18
    assert 0.0 < ss.C_in_M < ss.C_out_M


def test_strong_pump_clamps_intracellular_low():
    weak = _ss(1e-9, Vmax=1e-10, Km=1e-7)
    strong = _ss(1e-9, Vmax=1e-7, Km=1e-7)
    assert strong.C_in_M < weak.C_in_M
    assert strong.fold_reduction_vs_passive > weak.fold_reduction_vs_passive


def test_saturating_dose_overwhelms_the_pump():
    """At C_out >> Km the pump saturates, so C_in/C_out recovers upward."""
    Vmax, Km = 1e-9, 1e-9
    low = _ss(1e-11, Vmax, Km)     # well below Km
    high = _ss(1e-5, Vmax, Km)     # well above Km
    assert high.accumulation_ratio > low.accumulation_ratio
    assert high.pump_saturated_fraction > 0.9   # nearly maxed out
    assert low.pump_saturated_fraction < high.pump_saturated_fraction


def test_pump_from_michaelis_prior():
    from src.bridge import affinity_to_michaelis
    prior = affinity_to_michaelis(kcat_per_s=10.0, KM_M=2e-6)
    pump = pump_from_michaelis_prior(prior, pump_copies=1e5,
                                     cell_volume_L=4e-12)
    assert pump.Km_M == 2e-6
    assert pump.Vmax_M_per_s > 0
    # Michaelis pass-through prior is uncalibrated → pump uncalibrated.
    assert pump.trust == "uncalibrated"


def test_rejects_bad_inputs():
    try:
        EffluxPump(Vmax_M_per_s=1e-9, Km_M=0.0)
        assert False
    except ValueError:
        pass
    try:
        _ss(-1.0, Vmax=1e-9, Km=1e-7)
        assert False
    except ValueError:
        pass
