"""Sampled-phase binding smoke — first end-to-end test exercising
the full amber14 → sampler → MBAR → ΔG_bind pipeline.

Runs for 5-15 minutes on CPU, so it's OPT-IN: set environment
variable CELLSIM_RUN_SAMPLED_SMOKE=1 to enable. Otherwise it
prints 'SKIP' and returns 0 so CI doesn't burn a 15-minute slot
per PR.

Intended use: quick sanity-check after touching the binding
builder, sampler, or restraint correction. At these tiny params
the returned ΔG_bind will be numerical noise — the gate here is
'pipeline runs to completion without NaN', not 'ΔG accuracy'.

    CELLSIM_RUN_SAMPLED_SMOKE=1 \\
        python tests/fep/test_sampled_binding_smoke.py

For a real sampled run use:
    cellsim fep-binding bench benchmarks/fep/binding_egfr.yaml \\
        --padding 1.2 --sample --out-csv run.csv
    cellsim fep-report run.csv \\
        --yaml benchmarks/fep/binding_egfr.yaml
"""

from __future__ import annotations

import math
import os
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


def test_sampled_binding_pipeline_runs_to_completion():
    if not os.environ.get("CELLSIM_RUN_SAMPLED_SMOKE"):
        print("[SKIP] set CELLSIM_RUN_SAMPLED_SMOKE=1 to enable")
        return

    from src.fep.binding import compute_absolute_binding_dg

    t0 = time.time()
    r = compute_absolute_binding_dg(
        "C",
        REPO_ROOT / "benchmarks/md/1ubq.pdb",
        padding_nm=0.8,
        force_field_path="amber14",
        n_windows=3,
        n_equilibration_steps=100,
        n_production_steps=100,
        sample_stride=10,
        sample=True,
        seed=1,
    )
    wall = time.time() - t0

    assert r.ok, f"sampled pipeline failed: {r.reason}"
    assert r.phase == "sampled"
    assert r.dG_bind_kcalmol is not None
    assert not math.isnan(r.dG_bind_kcalmol), (
        "ΔG_bind is NaN — sampler hit an integrator blow-up")
    assert r.uncertainty_kcalmol is not None
    assert 0 < r.uncertainty_kcalmol < 50.0, (
        f"uncertainty {r.uncertainty_kcalmol} outside sane range")
    # At these toy parameters the ΔG will be noise (methane has no
    # real pocket on ubiquitin), but it should not be wildly huge
    # — a sign the restraint + restraint-correction have paired.
    assert abs(r.dG_bind_kcalmol) < 50.0, (
        f"ΔG_bind {r.dG_bind_kcalmol} unreasonably large; restraint "
        "or correction likely decoupled")
    print(f"[PASS] sampled binding pipeline  "
          f"wall={wall:.1f}s  "
          f"ΔG_bind={r.dG_bind_kcalmol:+.2f} ± "
          f"{r.uncertainty_kcalmol:.2f} kcal/mol  "
          f"(ΔG_dec_cx={r.dG_dec_complex_kcalmol:+.2f}, "
          f"ΔG_dec_sv={r.dG_dec_solvent_kcalmol:+.2f}, "
          f"ΔG_R={r.dG_restraint_correction_kcalmol:+.2f})")


def test_sampled_binding_with_restraint_on_real_leg():
    """Opt-in: exercise the restraint-on-real-bound leg (BUG_AUDIT #4)
    end-to-end. Validates that the 3-leg cycle runs to completion and
    populates the field / provenance flag. The NUMBER is noise at these
    toy params — production magnitude / k-invariance is a GPU gate."""
    if not os.environ.get("CELLSIM_RUN_SAMPLED_SMOKE"):
        print("[SKIP] set CELLSIM_RUN_SAMPLED_SMOKE=1 to enable")
        return

    from src.fep.binding import compute_absolute_binding_dg

    r = compute_absolute_binding_dg(
        "C",
        REPO_ROOT / "benchmarks/md/1ubq.pdb",
        padding_nm=0.8,
        force_field_path="amber14",
        n_windows=3,
        n_equilibration_steps=100,
        n_production_steps=100,
        sample_stride=10,
        sample=True,
        seed=1,
        include_restraint_on_real_leg=True,
    )
    assert r.ok, f"3-leg pipeline failed: {r.reason}"
    assert r.phase == "sampled"
    assert r.restraint_on_real_included is True, (
        "restraint-on-real leg should be marked included")
    assert r.dG_restraint_on_real_kcalmol is not None
    assert not math.isnan(r.dG_restraint_on_real_kcalmol)
    assert not math.isnan(r.dG_bind_kcalmol)
    print(f"[PASS] 3-leg sampled binding  "
          f"ΔG_bind={r.dG_bind_kcalmol:+.2f}  "
          f"ΔG_restraint_on_real={r.dG_restraint_on_real_kcalmol:+.2f} "
          f"(included={r.restraint_on_real_included})")


def test_default_leaves_restraint_on_real_provisional():
    """Without the flag, the result must flag itself provisional so no
    downstream consumer silently trusts the absolute value."""
    if not os.environ.get("CELLSIM_RUN_SAMPLED_SMOKE"):
        print("[SKIP] set CELLSIM_RUN_SAMPLED_SMOKE=1 to enable")
        return

    from src.fep.binding import compute_absolute_binding_dg

    r = compute_absolute_binding_dg(
        "C", REPO_ROOT / "benchmarks/md/1ubq.pdb",
        padding_nm=0.8, force_field_path="amber14",
        n_windows=3, n_equilibration_steps=100,
        n_production_steps=100, sample_stride=10,
        sample=True, seed=1)
    assert r.ok
    assert r.restraint_on_real_included is False
    assert r.dG_restraint_on_real_kcalmol is None
    assert "PROVISIONAL" in (r.reason or "")
    print(f"[PASS] default provisional: {r.reason[:80]}")


if __name__ == "__main__":
    try:
        test_sampled_binding_pipeline_runs_to_completion()
        test_sampled_binding_with_restraint_on_real_leg()
        test_default_leaves_restraint_on_real_provisional()
    except AssertionError as e:
        print(f"[FAIL] {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"[ERROR] {e}")
        sys.exit(2)
