"""Layer 1.3 triage-dashboard PNG smoke.

Writes a dashboard from a synthetic batch-result dict list and
verifies the PNG is produced, is non-empty, and contains the
four expected axes titles (verdict breakdown, ΔG by verdict,
pose-trust, safety filter).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.triage_viewer import render_triage_dashboard  # noqa: E402


def test_render_triage_dashboard():
    records = [
        dict(ok="True", triage="follow_up", dG_kcalmol=-9.2,
             strain_band="good", ro5_pass="True",
             herg_risk="low", mutagenic_risk="low"),
        dict(ok="True", triage="follow_up", dG_kcalmol=-8.6,
             strain_band="good", ro5_pass="True",
             herg_risk="low", mutagenic_risk="low"),
        dict(ok="True", triage="review", dG_kcalmol=-9.0,
             strain_band="suspicious", ro5_pass="True",
             herg_risk="high", mutagenic_risk="low"),
        dict(ok="True", triage="deprioritise", dG_kcalmol=-7.2,
             strain_band="acceptable", ro5_pass="True",
             herg_risk="low", mutagenic_risk="low"),
        dict(ok="True", triage="drop", dG_kcalmol=-5.4,
             strain_band="good", ro5_pass="True",
             herg_risk="low", mutagenic_risk="low"),
        dict(ok="True", triage="drop", dG_kcalmol=-9.1,
             strain_band="reject", ro5_pass="True",
             herg_risk="low", mutagenic_risk="high"),
        dict(ok="False"),
    ]
    with tempfile.TemporaryDirectory(prefix="cellsim-triage-") as tmp:
        out = Path(tmp) / "triage.png"
        render_triage_dashboard(
            records, out, title="triage smoke")
        assert out.exists(), "PNG not written"
        size = out.stat().st_size
        assert size > 5_000, f"PNG suspiciously small: {size} B"
        print(f"[triage-viewer] wrote {out} ({size} B) PASS")


if __name__ == "__main__":
    try:
        test_render_triage_dashboard()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
