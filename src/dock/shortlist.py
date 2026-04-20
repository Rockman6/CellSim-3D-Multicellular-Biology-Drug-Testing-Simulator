"""Extract a hand-to-wet-lab shortlist from a batch CSV.

Filters an existing batch-output CSV down to the rows whose
`triage` verdict is `follow_up` or `review` (i.e. the compounds
a chemist should actually look at), preserving the original
column order so the shortlist opens the same way in Excel /
pandas.

Useful when:
  - A docking run was done days ago without --shortlist-csv and
    you want to regenerate the handoff file.
  - You want to merge several batch CSVs into one combined
    shortlist (cat full_*.csv | grep 'follow_up\|review' — this
    module handles the header/row plumbing properly).

Non-AI. Pure CSV filter.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Optional


DEFAULT_VERDICTS = ("follow_up", "review")


def filter_csv(
    in_path: Path,
    out_path: Path,
    *,
    verdicts: tuple[str, ...] = DEFAULT_VERDICTS,
) -> int:
    """Filter `in_path` to `out_path`. Returns number of rows kept."""
    kept = 0
    with in_path.open(encoding="utf-8-sig") as fi:
        reader = csv.DictReader(fi)
        fieldnames = reader.fieldnames or []
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Write with UTF-8 BOM so Excel on Windows auto-detects the
        # encoding (same reasoning as src/dock/batch.py _write_csv).
        with out_path.open("w", newline="",
                            encoding="utf-8-sig") as fo:
            writer = csv.DictWriter(
                fo, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for row in reader:
                if row.get("triage") in verdicts:
                    writer.writerow(row)
                    kept += 1
    return kept


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path,
                    help="input batch CSV")
    ap.add_argument("--out", type=Path, required=True,
                    help="output shortlist CSV")
    ap.add_argument("--verdicts", default=",".join(DEFAULT_VERDICTS),
                    help="comma-separated triage verdicts to keep "
                         "(default: follow_up,review)")
    args = ap.parse_args(argv)

    verdicts = tuple(v.strip() for v in args.verdicts.split(",")
                      if v.strip())
    if not args.csv.exists():
        print(f"input CSV not found: {args.csv}", file=sys.stderr)
        return 2
    kept = filter_csv(args.csv, args.out, verdicts=verdicts)
    print(f"wrote {kept} rows → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
