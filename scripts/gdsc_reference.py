#!/usr/bin/env python3
"""Extract the GDSC2 reference IC50s for the curated lines and drugs.

GDSC2 fitted dose-response is a large public CSV (not vendored; see
``docs/PLAN.md`` Phase 1). Download it once:

    mkdir -p data/gdsc && cd data/gdsc
    curl -LO https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/\\
release-8.4/GDSC2_fitted_dose_response_24Jul22.csv

Then:

    python scripts/gdsc_reference.py --out benchmarks/cell/gdsc_reference.csv

The emitted CSV is the *target* of the Phase-1 validation: one row per
(drug, cell line) with the published IC50 in µM, the assay duration and
the screen's maximum concentration, which bounds how far an IC50 above
the top dose can be trusted.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.library import CELL_LINES, DRUGS  # noqa: E402

DEFAULT_CSV = REPO_ROOT / "data" / "gdsc" / "GDSC2_fitted_dose_response_24Jul22.csv"
# GDSC drug names differ in case from ours.
DRUG_ALIASES = {"cisplatin": "Cisplatin", "doxorubicin": "Doxorubicin",
                "paclitaxel": "Paclitaxel"}
# GDSC2 exposes cells to drug for 72 h before the CellTiter-Glo readout.
GDSC_ASSAY_HOURS = 72.0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV,
                    help=f"GDSC2 fitted dose-response CSV (default: {DEFAULT_CSV})")
    ap.add_argument("--out", type=Path, default=None, help="write a CSV here")
    a = ap.parse_args(argv)

    if not a.csv.exists():
        print(f"error: {a.csv} not found. See this script's docstring for the "
              f"one-line download.", file=sys.stderr)
        return 2

    want_drugs = {DRUG_ALIASES[d]: d for d in DRUGS}
    want_lines = set(CELL_LINES)
    rows: list[dict] = []
    with a.csv.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["DRUG_NAME"] in want_drugs and row["CELL_LINE_NAME"] in want_lines:
                rows.append({
                    "drug": want_drugs[row["DRUG_NAME"]],
                    "cell_line": row["CELL_LINE_NAME"],
                    "p53_functional": CELL_LINES[row["CELL_LINE_NAME"]].p53_functional,
                    "ic50_uM": round(math.exp(float(row["LN_IC50"])), 4),
                    "auc": round(float(row["AUC"]), 4),
                    "max_conc_uM": float(row["MAX_CONC"]),
                    "assay_hours": GDSC_ASSAY_HOURS,
                    "dataset": row["DATASET"],
                    "drug_id": row["DRUG_ID"],
                })
    rows.sort(key=lambda r: (r["drug"], r["cell_line"]))

    hdr = f"{'drug':<12}{'line':<12}{'p53':<6}{'IC50 µM':>10}{'AUC':>8}{'max µM':>9}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        flag = "wt" if r["p53_functional"] else "mut"
        over = " >top" if r["ic50_uM"] > r["max_conc_uM"] else ""
        print(f"{r['drug']:<12}{r['cell_line']:<12}{flag:<6}{r['ic50_uM']:>10.3g}"
              f"{r['auc']:>8.3f}{r['max_conc_uM']:>9.3g}{over}")
    print(f"\n{len(rows)} pairs. 'IC50 > top dose' means the screen never reached "
          f"50 % kill; treat as a lower bound.")

    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        with a.out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
