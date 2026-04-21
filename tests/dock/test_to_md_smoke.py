"""Layer 1.3 CSV → markdown table smoke."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.to_md import render_markdown  # noqa: E402


def test_render_markdown_small_table():
    records = [
        dict(ok="True", rank="1", name="compoundA",
             triage="follow_up", dG_kcalmol="-9.2",
             Kd_human="50 nM", strain_band="good",
             pocket_ok="True", ro5_pass="True", QED="0.72",
             triage_reason="ΔG -9.20, pose trustworthy"),
        dict(ok="True", rank="2", name="compoundB",
             triage="review", dG_kcalmol="-8.1",
             Kd_human="1.2 µM", strain_band="suspicious",
             pocket_ok="True", ro5_pass="True", QED="0.55",
             triage_reason="suspicious pose strain"),
        dict(ok="True", rank="3", name="compoundC",
             triage="drop", dG_kcalmol="-5.5",
             Kd_human="95 µM", strain_band="good",
             pocket_ok="True", ro5_pass="True", QED="0.41",
             triage_reason="ΔG > -6"),
    ]
    md = render_markdown(records, top=3)
    # Header, separator, 3 data rows = 5 lines in the table block.
    assert "| rank |" in md
    assert "| 1 | compoundA | follow_up |" in md
    assert "| 2 | compoundB | review |" in md
    assert "| 3 | compoundC | drop |" in md
    # Next-steps block
    assert "## Next steps" in md
    assert "Send **1** follow_up" in md
    assert "Re-examine **1** review" in md
    print("  markdown render ok")


def test_render_markdown_empty_hits():
    records = [
        dict(ok="True", rank="1", name="x", triage="deprioritise",
             dG_kcalmol="-7.0"),
    ]
    md = render_markdown(records, top=3)
    assert "No compounds qualify" in md
    print("  empty-hits warning ok")


if __name__ == "__main__":
    tests = [test_render_markdown_small_table,
             test_render_markdown_empty_hits]
    for t in tests:
        try:
            t()
            print(f"  {t.__name__} PASS")
        except AssertionError as e:
            print(f"  {t.__name__} FAIL: {e}")
            sys.exit(1)
    print("PASS")
