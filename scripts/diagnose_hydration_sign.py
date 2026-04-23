#!/usr/bin/env python3
"""diagnose_hydration_sign.py — reproduce the Milestone A failure
mode on a laptop in ~5 min.

Background: the friend's M5 Max run of milestone-a-pilot-1 returned
a systematic over-negative bias across 7 compounds. We expected
this might be OpenCL-specific, but methane + ethanol reproduce the
failure on CPU/HEAD with 7 windows × 1000 prod steps, ruling out
the platform hypothesis.

This script makes the failure mode visible in ~5 min of laptop
compute by printing both legs of the hydration cycle plus the
final ΔG_hyd. Anyone reviewing the diagnosis can run it without
needing a GPU or the friend's tarball.

Observed pattern:
  methane:  vac=+0.000, solv=+1.813 → pred -1.81 (expt +2.00)
    sign-flipping the composition rescues methane (+1.81 ≈ +2.00).
  ethanol:  vac=+2.36,  solv=+15.55 → pred -13.19 (expt -5.01)
    sign-flip alone does NOT rescue ethanol — the solvent-leg
    magnitude is ~3× what the physics predicts for a polar solute.

Two bugs, one composition-level sign flip and one solvent-leg
magnitude issue for polar/aromatic compounds. See BENCHMARKS.md
§ "Milestone A post-mortem" for the full analysis.

Usage:
    python scripts/diagnose_hydration_sign.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

CASES = [
    # (name, smiles, expt_kcalmol, expected_sign)
    ("methane",  "C",    +2.00, "positive (hydrophobic)"),
    ("ethanol",  "CCO",  -5.01, "negative (H-bond favored)"),
]

PARAMS = dict(
    n_windows=7, n_production_steps=1000,
    n_equilibration_steps=200, sample_stride=100, seed=1)


def main() -> int:
    from src.fep import compute_hydration_dg

    print("=" * 70)
    print("Milestone A hydration sign/magnitude diagnostic")
    print(f"Params: {PARAMS}")
    print("=" * 70)
    fails = 0
    for name, smiles, expt, expected_sign in CASES:
        print(f"\n[{name}]  SMILES={smiles}  expt={expt:+.2f} kcal/mol")
        r = compute_hydration_dg(smiles, **PARAMS)
        if not r.ok:
            print(f"  ERROR: {r.reason}")
            fails += 1
            continue
        pred = r.dG_hydration_kcalmol
        vac = r.dG_vacuum_decouple_kcalmol
        solv = r.dG_solvent_decouple_kcalmol
        print(f"  vac_decouple   = {vac:+.3f}")
        print(f"  solv_decouple  = {solv:+.3f}")
        print(f"  pred (as-is)   = {pred:+.3f}")
        print(f"  pred (sign-flip) = {-pred:+.3f}")
        print(f"  expected sign  = {expected_sign}")
        ok_as_is = (pred > 0) == (expt > 0)
        ok_flipped = (-pred > 0) == (expt > 0)
        print(f"  sign match as-is   : {'OK' if ok_as_is else 'FAIL'}")
        print(f"  sign match flipped : {'OK' if ok_flipped else 'FAIL'}")
        resid_flipped = abs((-pred) - expt)
        if resid_flipped > 3.0:
            print(f"  ! sign-flip alone doesn't rescue "
                  f"{name}: |resid|={resid_flipped:.2f} > 3 kcal/mol")
        if not ok_as_is:
            fails += 1
    print()
    print("=" * 70)
    print(f"{len(CASES) - fails}/{len(CASES)} compounds produced the "
          "correct sign under the current composition formula.")
    print("If that's < 2/2 you're seeing the Milestone A bug. See "
          "BENCHMARKS.md §'Milestone A post-mortem'.")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
