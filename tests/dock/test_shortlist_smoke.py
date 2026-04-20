"""Layer 1.3 shortlist CSV filter smoke."""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.shortlist import filter_csv  # noqa: E402


def test_filter_csv_keeps_followup_and_review():
    rows = [
        dict(rank="1", name="A", triage="follow_up",
             dG_kcalmol="-9.2"),
        dict(rank="2", name="B", triage="review",
             dG_kcalmol="-8.6"),
        dict(rank="3", name="C", triage="deprioritise",
             dG_kcalmol="-7.2"),
        dict(rank="4", name="D", triage="drop",
             dG_kcalmol="-5.4"),
        dict(rank="5", name="E", triage="",
             dG_kcalmol="-6.0"),
    ]
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        in_p = td / "in.csv"
        out_p = td / "out.csv"
        with in_p.open("w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

        kept = filter_csv(in_p, out_p)
        assert kept == 2, f"expected 2 kept; got {kept}"
        out_rows = list(csv.DictReader(out_p.open()))
        assert [r["name"] for r in out_rows] == ["A", "B"]
        assert [r["triage"] for r in out_rows] == \
            ["follow_up", "review"]
        print(f"  kept {kept} rows, order preserved PASS")


def test_filter_csv_custom_verdicts():
    rows = [
        dict(rank="1", name="A", triage="follow_up"),
        dict(rank="2", name="B", triage="drop"),
        dict(rank="3", name="C", triage="deprioritise"),
    ]
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        in_p = td / "in.csv"
        out_p = td / "out.csv"
        with in_p.open("w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

        kept = filter_csv(in_p, out_p,
                           verdicts=("follow_up", "deprioritise"))
        assert kept == 2
        names = [r["name"] for r in csv.DictReader(out_p.open())]
        assert names == ["A", "C"]
        print(f"  custom verdicts kept {names} PASS")


if __name__ == "__main__":
    try:
        test_filter_csv_keeps_followup_and_review()
        test_filter_csv_custom_verdicts()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
