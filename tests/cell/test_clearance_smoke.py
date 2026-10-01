"""Metabolic clearance and the PK/PD coupling.

Anchors:
  * AUC = Dose/CL — the fundamental PK identity, checked against a
    numerically integrated trajectory;
  * hepatic clearance can never exceed hepatic blood flow (flow limit),
    and reduces to fu·CL_int in the capacity-limited regime;
  * accumulation on repeat dosing matches 1/(1−e^{−kτ});
  * the payoff: a potent drug that is cleared fast never sustains the
    critical occupancy, so it FAILS despite the potency.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402
from src.cell import (  # noqa: E402
    HepaticClearance,
    single_dose_exposure,
    repeat_dose_exposure,
    dosing_accumulation_ratio,
    evaluate_regimen,
)


def _prior(dG):
    return binding_to_hill(dG, uncertainty_kcalmol=0.2,
                           receptor="benchmarks/dock/3ptb.pdb", method="vina")


def test_auc_equals_dose_over_clearance():
    """THE PK identity: AUC = Dose/CL, independent of V_d."""
    cl = HepaticClearance(CLint_L_per_h=20.0, fu_plasma=0.5, Vd_L=50.0)
    dose = 1e-3   # mol
    exp = single_dose_exposure(dose, cl, duration_h=400.0, n_points=40000)
    expected = dose / cl.CL_hepatic_L_per_h
    assert abs(exp.auc_M_h() - expected) / expected < 1e-3


def test_clearance_never_exceeds_hepatic_blood_flow():
    """Flow limit: however fast the enzyme, blood must deliver the drug."""
    huge = HepaticClearance(CLint_L_per_h=1e9, fu_plasma=1.0)
    assert huge.CL_hepatic_L_per_h < huge.Q_hepatic_L_per_h
    assert huge.CL_hepatic_L_per_h / huge.Q_hepatic_L_per_h > 0.99
    assert huge.extraction_ratio < 1.0
    assert huge.bioavailability_first_pass > 0.0


def test_capacity_limited_regime_reduces_to_fu_times_clint():
    small = HepaticClearance(CLint_L_per_h=0.05, fu_plasma=0.4)
    approx = small.fu_plasma * small.CLint_L_per_h
    assert abs(small.CL_hepatic_L_per_h - approx) / approx < 1e-3


def test_half_life_and_first_order_decay():
    cl = HepaticClearance(CLint_L_per_h=10.0, fu_plasma=1.0, Vd_L=70.0)
    t_half = cl.half_life_h
    exp = single_dose_exposure(1e-6, cl, duration_h=t_half, n_points=2001)
    # Concentration halves over one half-life.
    assert abs(exp.C_M[-1] / exp.C_M[0] - 0.5) < 1e-3


def test_repeat_dose_accumulation_matches_formula():
    cl = HepaticClearance(CLint_L_per_h=8.0, fu_plasma=1.0, Vd_L=70.0)
    tau = 12.0
    ratio = dosing_accumulation_ratio(cl, tau)
    assert ratio > 1.0
    exp = repeat_dose_exposure(1e-6, cl, interval_h=tau, n_doses=40)
    C0_single = 1e-6 / cl.Vd_L
    # Peak at steady state ≈ C0 × accumulation ratio.
    assert abs(exp.Cmax_M / (C0_single * ratio) - 1.0) < 0.02


def test_fast_clearance_defeats_a_potent_drug():
    """The payoff: potency without sustained exposure is not efficacy."""
    prior = _prior(-11.0)                      # very potent, K_d ~8 nM
    dose = 1e-4

    slow = HepaticClearance(CLint_L_per_h=0.5, fu_plasma=1.0, Vd_L=70.0)
    fast = HepaticClearance(CLint_L_per_h=5000.0, fu_plasma=1.0, Vd_L=70.0)

    out_slow = evaluate_regimen(
        single_dose_exposure(dose, slow, duration_h=168.0), prior)
    out_fast = evaluate_regimen(
        single_dose_exposure(dose, fast, duration_h=168.0), prior)

    # Same drug, same dose — only the clearance differs.
    assert out_slow.hours_above_threshold > out_fast.hours_above_threshold
    assert out_slow.log10_kill > out_fast.log10_kill
    assert out_slow.net_effect == "kills"
    assert out_fast.net_effect in ("holds", "fails")


def test_rejects_bad_inputs():
    for bad in (dict(CLint_L_per_h=-1.0), dict(CLint_L_per_h=1.0, Vd_L=0.0),
                dict(CLint_L_per_h=1.0, fu_plasma=0.0)):
        try:
            HepaticClearance(**bad)
            assert False, bad
        except ValueError:
            pass
