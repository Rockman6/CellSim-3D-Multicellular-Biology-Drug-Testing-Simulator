"""Layer 1.3 batch-triage verdict smoke.

Verifies the _triage_call rules that synthesise (ΔG, strain,
pocket, ADMET) into a single follow_up / review / deprioritise /
drop verdict.

These are pure-python rules over the batch record dict, no
external state. Gate ensures the logic doesn't drift.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.batch import _triage_call  # noqa: E402


CASES = [
    # Strong ΔG, clean everything → follow up.
    (dict(dG_kcalmol=-9.2, strain_band="good",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="low", ro5_violations=0),
     "follow_up"),
    # Strong ΔG but suspicious strain → review.
    (dict(dG_kcalmol=-9.0, strain_band="suspicious",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="low", ro5_violations=0),
     "review"),
    # Strong ΔG but high mutagenicity → drop.
    (dict(dG_kcalmol=-9.5, strain_band="good",
          pocket_ok=True, mutagenic_risk="high",
          herg_risk="low", ro5_violations=0),
     "drop"),
    # Strain=reject → drop regardless.
    (dict(dG_kcalmol=-10.0, strain_band="reject",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="low", ro5_violations=0),
     "drop"),
    # Weak ΔG > -6 → drop.
    (dict(dG_kcalmol=-5.2, strain_band="good",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="low", ro5_violations=0),
     "drop"),
    # Borderline ΔG -7 with hERG flag → deprioritise.
    (dict(dG_kcalmol=-7.0, strain_band="good",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="high", ro5_violations=0),
     "deprioritise"),
    # Strong ΔG with hERG flag → review (worth human time).
    (dict(dG_kcalmol=-9.0, strain_band="good",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="high", ro5_violations=0),
     "review"),
    # Borderline ΔG, nothing flagged → deprioritise by ΔG tier.
    (dict(dG_kcalmol=-7.2, strain_band="good",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="low", ro5_violations=0),
     "deprioritise"),
    # 3 Ro5 violations → drop.
    (dict(dG_kcalmol=-9.0, strain_band="good",
          pocket_ok=True, mutagenic_risk="low",
          herg_risk="low", ro5_violations=3),
     "drop"),
    # No ΔG at all → drop.
    (dict(dG_kcalmol=None, strain_band=None), "drop"),
]


def test_triage_rules():
    for i, (rec, expected) in enumerate(CASES):
        call, why = _triage_call(rec)
        assert call == expected, (
            f"case {i}: {rec}\n  expected={expected}, got={call} "
            f"(reason={why!r})")
        print(f"  case {i}: {call:<13s} — {why}")


if __name__ == "__main__":
    try:
        test_triage_rules()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
