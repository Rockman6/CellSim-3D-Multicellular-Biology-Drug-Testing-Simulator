"""src/bridge — Campaign-1 → Campaign-2 rate-law emitter smoke.

Pin the closed-form thermodynamic conversions so a future
refactor of the Hill / Michaelis primitives can't silently
drift the K_d / CI numbers Campaign-2 will consume.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import (  # noqa: E402
    RateLawPrior,
    binding_to_hill,
    affinity_to_michaelis,
)


def test_binding_to_hill_kd_matches_closed_form():
    """ΔG = -10 kcal/mol at 298.15 K should give K_d ≈ 50 nM
    (the canonical 'good drug' affinity)."""
    prior = binding_to_hill(-10.0)
    Kd_M = prior.parameters["Kd_M"]
    # K_d = exp(-10 / (R*T)) at T=298.15
    expected = math.exp(-10.0 / (0.0019872041 * 298.15))
    assert abs(Kd_M - expected) < 1e-15
    # Sanity: ~50 nM (5e-8 M)
    assert 1e-8 < Kd_M < 1e-7, f"K_d {Kd_M} M not in nM range"


def test_binding_to_hill_uncertainty_propagates_multiplicatively():
    """1σ on ΔG → multiplicative CI on K_d (since K_d is
    exponential in ΔG)."""
    prior = binding_to_hill(-10.0, uncertainty_kcalmol=0.5)
    Kd_M = prior.parameters["Kd_M"]
    Kd_lo, Kd_hi = prior.parameter_ci95["Kd_M"]

    # 95% CI is ±1.96σ in ΔG → exp(±1.96σ/RT) in K_d.
    sigma_ln = 0.5 / (0.0019872041 * 298.15)
    expected_lo = Kd_M * math.exp(-1.96 * sigma_ln)
    expected_hi = Kd_M * math.exp(+1.96 * sigma_ln)
    assert abs(Kd_lo - expected_lo) / expected_lo < 1e-12
    assert abs(Kd_hi - expected_hi) / expected_hi < 1e-12

    # CI must be asymmetric on a multiplicative scale.
    # log10(Kd_hi/Kd_M) should equal log10(Kd_M/Kd_lo).
    geom_lo = math.log10(Kd_M / Kd_lo)
    geom_hi = math.log10(Kd_hi / Kd_M)
    assert abs(geom_lo - geom_hi) < 1e-12


def test_binding_to_hill_no_uncertainty_no_ci():
    """No σ → no CI populated. Campaign-2 consumer can branch
    on `bool(prior.parameter_ci95)` to know if estimate is
    calibrated."""
    prior = binding_to_hill(-10.0)
    assert prior.parameter_ci95 == {}
    assert prior.source == "FEP"  # default method starts with 'FEP'
    assert prior.method == "FEP-DDM-MBAR"


def test_binding_to_hill_phenomenological_source_tag():
    """When the input came from the OLD/ Tier-0 BindingMatcher
    (descriptor-fit), the consumer needs to know the prior is
    NOT physics-grounded — propagate that via source='phenomenological'."""
    prior = binding_to_hill(
        -7.0, method="phenomenological-Lipinski",
        uncertainty_kcalmol=2.5)  # inflated σ for the descriptor fit
    assert prior.source == "phenomenological"
    assert "phenomenological" in prior.method


def test_affinity_to_michaelis_passes_through():
    """Pass-through impl for now. Just verify the dataclass shape
    is sane so Campaign-2 loaders can unpack consistently."""
    prior = affinity_to_michaelis(
        kcat_per_s=100.0, KM_M=1e-5,
        uncertainty_kcat=10.0, uncertainty_KM=2e-6)
    assert prior.type == "michaelis"
    assert prior.parameters == {"kcat_per_s": 100.0, "KM_M": 1e-5}
    assert "kcat_per_s" in prior.parameter_ci95
    assert "KM_M" in prior.parameter_ci95
    # Default method is 'literature'
    assert prior.source == "literature"


def test_summary_one_line():
    prior = binding_to_hill(-12.0, uncertainty_kcalmol=0.4)
    s = prior.summary()
    assert s.startswith("[hill]")
    assert "src=FEP" in s
    assert "T=298.1K" in s


if __name__ == "__main__":
    funcs = [
        test_binding_to_hill_kd_matches_closed_form,
        test_binding_to_hill_uncertainty_propagates_multiplicatively,
        test_binding_to_hill_no_uncertainty_no_ci,
        test_binding_to_hill_phenomenological_source_tag,
        test_affinity_to_michaelis_passes_through,
        test_summary_one_line,
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
