"""Convert a batch-output CSV to a markdown table for lab-notebook
paste (Benchling / eLabFTW / Obsidian / any markdown host).

Wet-lab rationale: biologists keep notes in markdown-capable
tools and don't want to format a ΔG table by hand. This emits
a GitHub-flavoured markdown table of the top-N hits plus the
Next-steps verdict counts — paste-ready.

Usage:
    cellsim to-md run.csv --out report.md --top 10
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Optional


def _row(rec: dict) -> list[str]:
    def _pick(key: str, default: str = "—") -> str:
        v = rec.get(key)
        return str(v) if v not in (None, "", "None") else default
    return [
        _pick("rank"),
        _pick("name"),
        _pick("triage"),
        _pick("dG_kcalmol"),
        _pick("Kd_human"),
        _pick("strain_band"),
        _pick("pocket_ok"),
        _pick("ro5_pass"),
        _pick("QED"),
        _pick("triage_reason"),
    ]


def render_markdown(
    records: list[dict],
    *,
    top: Optional[int] = 10,
    title: str = "CellSim screen — triage summary",
) -> str:
    header = [
        "rank", "name", "triage", "ΔG (kcal/mol)", "K_d",
        "strain", "pocket", "Ro5", "QED", "reason"]
    sep = [":--" for _ in header]

    # Verdict counts (from all ok records, not just the top N).
    from collections import Counter
    verdicts = Counter(
        r.get("triage") for r in records if r.get("ok") == "True")

    subset = records[:top] if top else records

    lines = [
        f"# {title}",
        "",
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(sep) + " |",
    ]
    for rec in subset:
        lines.append("| " + " | ".join(_row(rec)) + " |")

    # Next-steps block.
    n_f = verdicts.get("follow_up", 0)
    n_r = verdicts.get("review", 0)
    n_d = verdicts.get("deprioritise", 0)
    n_x = verdicts.get("drop", 0)
    lines.extend([
        "",
        "## Next steps",
        "",
    ])
    if n_f:
        lines.append(
            f"- Send **{n_f}** follow_up compound(s) to wet lab.")
    if n_r:
        lines.append(
            f"- Re-examine **{n_r}** review compound(s) — "
            "one flag needs a chemist's eye before committing.")
    if n_d:
        lines.append(
            f"- {n_d} deprioritise compound(s): rescore with FEP "
            "before triage if scaffold-important.")
    if n_x:
        lines.append(
            f"- {n_x} drop compound(s): too weak / non-physical "
            "pose / known-bad ADMET.")
    if not (n_f or n_r):
        lines.append(
            "- ⚠ No compounds qualify for wet-lab follow-up in "
            "this batch. Check the target-class reliability table "
            "(TUTORIAL §8) before blaming the library.")
    return "\n".join(lines) + "\n"


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path,
                    help="batch output CSV")
    ap.add_argument("--out", type=Path,
                    help="markdown path (default: print to stdout)")
    ap.add_argument("--top", type=int, default=10,
                    help="emit the top-N rows in the table "
                         "(default 10; 0 = all)")
    ap.add_argument("--title", default="CellSim screen — "
                                         "triage summary")
    args = ap.parse_args(argv)

    if not args.csv.exists():
        print(f"input CSV not found: {args.csv}", file=sys.stderr)
        return 2
    with args.csv.open(encoding="utf-8-sig") as fp:
        records = list(csv.DictReader(fp))

    top = args.top if args.top > 0 else None
    md = render_markdown(records, top=top, title=args.title)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(md, encoding="utf-8")
        print(f"wrote {args.out}")
    else:
        print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
