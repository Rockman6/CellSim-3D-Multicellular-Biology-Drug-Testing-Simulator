"""Regression tests for the biologist-readable MBAR-failure
translator in src/fep/sampling.py.

The underlying problem: pymbar's built-in error strings ("Column
sum W_nk = 0 for state 3 and 8 other columns") are meaningful to
statistical mechanics experts and meaningless to a wet-lab biologist
running `cellsim fep-binding` for the first time. After a 6-hour
GPU burn, the biologist needs to know (a) *why* it failed and (b)
*what parameter to change* — not what a reweighting matrix is.

These tests pin the translator behaviour so we never regress to
the raw pymbar string."""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.fep.sampling import _biologist_reason_for_mbar_error


def _msg(exc_text: str, **kwargs) -> str:
    defaults = dict(n_windows=11, n_production_steps=1500)
    defaults.update(kwargs)
    return _biologist_reason_for_mbar_error(
        Exception(exc_text), **defaults)


def test_column_sum_triggers_overlap_message():
    out = _msg(
        "Column sum W_nk = 0 for state 0 and 4 other columns; "
        "check that the states have sufficient overlap")
    assert "λ-windows don't overlap" in out, out
    assert "tight binders" in out, out
    assert "--n-windows" in out, out
    assert "--n-production-steps" in out, out


def test_w_nk_phrase_triggers_overlap_message():
    # Some pymbar versions lowercase the matrix name.
    out = _msg("w_nk matrix has columns that do not overlap")
    assert "λ-windows don't overlap" in out


def test_singular_matrix_triggers_convergence_message():
    out = _msg("Singular matrix encountered in MBAR iteration")
    assert "did not converge" in out or "convergence" in out.lower(), out
    assert "--n-windows" in out


def test_nan_reduced_potentials_mentions_integrator():
    out = _msg("Reduced potential contains NaN values")
    assert "NaN" in out or "integrator" in out.lower(), out
    assert "timestep" in out.lower(), out


def test_bounds_error_mentions_samples():
    out = _msg("BoundsError: squared uncertainty is negative")
    assert "uncertainty" in out.lower()
    assert "samples per window" in out.lower() or (
        "samples" in out.lower())


def test_fallback_preserves_original_text():
    out = _msg("Mystery error the translator hasn't seen")
    # Fallback must include a truncated copy of the original so
    # debug info isn't lost, but also still suggest bumped params.
    assert "Mystery error" in out
    assert "--n-windows" in out


def test_suggested_windows_bumped_at_least_21():
    # Tight-binder-tier floor: 21 windows (Milestone-A practice).
    out = _msg("column sum 0", n_windows=11, n_production_steps=1500)
    assert "--n-windows 22" in out or "--n-windows 21" in out, out
    assert "--n-production-steps" in out


def test_suggested_prod_capped_at_100k():
    # Don't print scare numbers — 100k is a reasonable ceiling
    # for the suggestion text regardless of user input.
    out = _msg("column sum 0", n_windows=11,
               n_production_steps=50000)
    # Bump would be 3× = 150000, but capped at 100000.
    assert "--n-production-steps 100000" in out, out


def test_suggested_windows_capped_at_31():
    out = _msg("column sum 0", n_windows=21,
               n_production_steps=1500)
    # 2× = 42, capped at 31.
    assert "--n-windows 31" in out, out


def test_gpu_recommended_for_tight_binders():
    out = _msg("Column sum W_nk = 0")
    assert "GPU" in out


if __name__ == "__main__":
    funcs = [
        test_column_sum_triggers_overlap_message,
        test_w_nk_phrase_triggers_overlap_message,
        test_singular_matrix_triggers_convergence_message,
        test_nan_reduced_potentials_mentions_integrator,
        test_bounds_error_mentions_samples,
        test_fallback_preserves_original_text,
        test_suggested_windows_bumped_at_least_21,
        test_suggested_prod_capped_at_100k,
        test_suggested_windows_capped_at_31,
        test_gpu_recommended_for_tight_binders,
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
