"""Layer 1.3 CYP3A4 inhibition strain-aware risk downgrade smoke.

Unit-tests the classification-plus-downgrade logic by directly
constructing CypInhibitionResult objects with known (ΔG, Fe-dist,
strain_band) combinations and calling the downgrade step.

Verifies:
  - strain=acceptable with high-risk geometry keeps high
  - strain=reject with high-risk geometry drops to low
  - strain=reject with low-risk geometry stays low
  - strain=suspicious does NOT downgrade (only reject does;
    suspicious is a review signal, not a trust-rejection)
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.cyp_inhibition import _classify  # noqa: E402


def _apply_downgrade(dG, fe, strain_band, strain_ratio=3.0):
    """Mirror the downgrade step in cyp3a4_inhibition."""
    risk, reason = _classify(dG, fe)
    if strain_band == "reject" and risk in ("high", "medium"):
        reason = (f"risk downgraded from {risk} → low: strain="
                  f"{strain_band} (ratio {strain_ratio:.2f}); "
                  "non-physical pose, ΔG not trustworthy")
        risk = "low"
    return risk, reason


def test_acceptable_strain_keeps_high():
    risk, _ = _apply_downgrade(
        dG=-9.85, fe=2.13, strain_band="acceptable")
    assert risk == "high"


def test_reject_strain_downgrades_high_to_low():
    risk, reason = _apply_downgrade(
        dG=-9.0, fe=3.0, strain_band="reject", strain_ratio=8.5)
    assert risk == "low"
    assert "downgraded" in reason
    assert "reject" in reason


def test_reject_strain_downgrades_medium_to_low():
    risk, _ = _apply_downgrade(
        dG=-7.2, fe=6.0, strain_band="reject", strain_ratio=7.4)
    assert risk == "low"


def test_reject_on_already_low_stays_low():
    risk, reason = _apply_downgrade(
        dG=-5.0, fe=10.0, strain_band="reject")
    assert risk == "low"
    assert "downgraded" not in reason


def test_suspicious_does_not_downgrade():
    risk, _ = _apply_downgrade(
        dG=-9.0, fe=2.0, strain_band="suspicious", strain_ratio=5.0)
    assert risk == "high"


def test_no_strain_data_keeps_classification():
    risk, _ = _apply_downgrade(
        dG=-9.0, fe=2.0, strain_band=None)
    assert risk == "high"


if __name__ == "__main__":
    tests = [
        test_acceptable_strain_keeps_high,
        test_reject_strain_downgrades_high_to_low,
        test_reject_strain_downgrades_medium_to_low,
        test_reject_on_already_low_stays_low,
        test_suspicious_does_not_downgrade,
        test_no_strain_data_keeps_classification,
    ]
    for t in tests:
        try:
            t()
            print(f"  {t.__name__} PASS")
        except AssertionError as e:
            print(f"  {t.__name__} FAIL: {e}")
            sys.exit(1)
    print("PASS")
