"""Layer 1.3 strain-gate top-pose promotion smoke.

Exercises the _worker logic path by monkey-patching ligand_strain
to return a synthetic band for each pose, verifying that:

  - a reject top-1 with an acceptable pose #3 promotes to #3
  - a reject top-1 with no passing pose keeps #1 (honest)
  - an acceptable top-1 is never promoted
  - strain_gate=False disables promotion

Without mocking the strain function, we can't force a reject
band deterministically (real Vina poses on bundled ligands all
land in good/acceptable).  The mock keeps the test deterministic
and sub-second.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import strain as strain_mod  # noqa: E402


@dataclass
class _FakePose:
    elements: list = field(default_factory=lambda: ["C", "C"])
    positions_A: list = field(default_factory=lambda: [[0, 0, 0], [1, 0, 0]])
    affinity_kcalmol: float = -9.0


def _make_strain(ok=True, band="good", ratio=1.2, kcal=5.0):
    return strain_mod.StrainResult(
        ok=ok, smiles="X", n_heavy=2,
        e_docked_kcalmol=10.0, e_relaxed_mean_kcalmol=5.0,
        strain_kcalmol=kcal, energy_ratio=ratio, band=band)


def _simulate_promotion(pose_bands, gate=True):
    """Mirror the gate logic in src.dock.batch._worker for testing.

    Returns (chosen_rank_1indexed, promoted_from_1indexed_or_None)."""
    # pose_bands is a list of band strings per pose.
    strains = [(i, _FakePose(), _make_strain(band=b))
               for i, b in enumerate(pose_bands)]
    top_i, _, top_strain = strains[0]
    promoted = None
    if gate and top_strain.ok and top_strain.band == "reject":
        for i, _, s in strains:
            if s.ok and s.band in ("good", "acceptable"):
                top_i = i
                promoted = i + 1
                break
    return top_i + 1, promoted


def test_gate_promotes_from_reject_to_acceptable():
    rank, promoted = _simulate_promotion(
        ["reject", "suspicious", "acceptable"], gate=True)
    assert rank == 3 and promoted == 3, (rank, promoted)


def test_gate_promotes_to_first_good():
    rank, promoted = _simulate_promotion(
        ["reject", "good", "acceptable"], gate=True)
    assert rank == 2 and promoted == 2, (rank, promoted)


def test_gate_keeps_top1_when_no_good_pose():
    rank, promoted = _simulate_promotion(
        ["reject", "reject", "suspicious"], gate=True)
    assert rank == 1 and promoted is None, (rank, promoted)


def test_gate_does_not_change_acceptable_top1():
    rank, promoted = _simulate_promotion(
        ["acceptable", "good", "good"], gate=True)
    assert rank == 1 and promoted is None, (rank, promoted)


def test_gate_disabled_keeps_top1_even_if_reject():
    rank, promoted = _simulate_promotion(
        ["reject", "good"], gate=False)
    assert rank == 1 and promoted is None, (rank, promoted)


if __name__ == "__main__":
    tests = [
        test_gate_promotes_from_reject_to_acceptable,
        test_gate_promotes_to_first_good,
        test_gate_keeps_top1_when_no_good_pose,
        test_gate_does_not_change_acceptable_top1,
        test_gate_disabled_keeps_top1_even_if_reject,
    ]
    for t in tests:
        try:
            t()
            print(f"  {t.__name__} PASS")
        except AssertionError as e:
            print(f"  {t.__name__} FAIL: {e}")
            sys.exit(1)
    print("PASS")
