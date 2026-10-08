#!/usr/bin/env python3
"""Extract GDSC reference IC50s for the curated lines and drugs.

GDSC fitted dose-response tables are large public CSVs (not vendored; see
``docs/PLAN.md`` Phase 1). Download them once (GDSC1 holds doxorubicin):

    mkdir -p data/gdsc && cd data/gdsc
    curl -LO https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/\\
release-8.4/GDSC2_fitted_dose_response_24Jul22.csv
    curl -LO https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/\\
release-8.4/GDSC1_fitted_dose_response_24Jul22.csv

Then:

    python scripts/gdsc_reference.py --out benchmarks/cell/gdsc_reference.csv

The emitted CSV is the *target* of the GDSC validation: one row per
(drug, cell line) with the published IC50 in µM, the assay duration and
the screen's maximum concentration, which bounds how far an IC50 above
the top dose can be trusted.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.library import CELL_LINES, DRUGS  # noqa: E402

DEFAULT_CSVS = [
    REPO_ROOT / "data" / "gdsc" / "GDSC2_fitted_dose_response_24Jul22.csv",
    REPO_ROOT / "data" / "gdsc" / "GDSC1_fitted_dose_response_24Jul22.csv",
]
# GDSC drug names differ from ours in case and, for nutlin, form.
DRUG_ALIASES = {"cisplatin": "Cisplatin", "doxorubicin": "Doxorubicin",
                "paclitaxel": "Paclitaxel", "etoposide": "Etoposide", "sn-38": "SN-38",
                "gemcitabine": "Gemcitabine", "5-fluorouracil": "5-Fluorouracil",
                "docetaxel": "Docetaxel", "vinorelbine": "Vinorelbine",
                "nutlin-3a": "Nutlin-3a (-)", "palbociclib": "Palbociclib"}
# GDSC2 exposes cells to drug for 72 h before the CellTiter-Glo readout.
GDSC_ASSAY_HOURS = 72.0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, nargs="+", default=DEFAULT_CSVS,
                    help="GDSC fitted dose-response CSV(s); default: GDSC2 and GDSC1 "
                         "release 8.4 under data/gdsc/ (missing files are skipped)")
    ap.add_argument("--out", type=Path, default=None, help="write a CSV here")
    a = ap.parse_args(argv)

    present = [p for p in a.csv if p.exists()]
    if not present:
        print(f"error: none of {[str(p) for p in a.csv]} found. See this script's "
              f"docstring for the one-line download.", file=sys.stderr)
        return 2

    want_drugs = {DRUG_ALIASES[d]: d for d in DRUGS}
    # Match on each line's GDSC name. Exact matching on our own names
    # silently dropped HCT116 (GDSC: HCT-116) until October 2026.
    want_lines = {(cl.gdsc_name or name): name for name, cl in CELL_LINES.items()}
    rows: list[dict] = []
    for path in present:
        with path.open(newline="") as f:
            for row in csv.DictReader(f):
                if row["DRUG_NAME"] in want_drugs and row["CELL_LINE_NAME"] in want_lines:
                    ic50 = math.exp(float(row["LN_IC50"]))
                    top = float(row["MAX_CONC"])
                    ours = want_lines[row["CELL_LINE_NAME"]]
                    rows.append({
                        "drug": want_drugs[row["DRUG_NAME"]],
                        "cell_line": ours,
                        "p53_functional": CELL_LINES[ours].p53_functional,
                        "ic50_uM": round(ic50, 4),
                        "auc": round(float(row["AUC"]), 4),
                        "max_conc_uM": top,
                        "ic50_in_tested_range": ic50 <= top,
                        "assay_hours": GDSC_ASSAY_HOURS,
                        "dataset": row["DATASET"],
                        "drug_id": row["DRUG_ID"],
                    })
    rows.sort(key=lambda r: (r["drug"], r["cell_line"], r["dataset"], r["drug_id"]))

    hdr = (f"{'drug':<12}{'line':<12}{'p53':<5}{'set':<7}{'id':>6}"
           f"{'IC50 µM':>10}{'AUC':>7}{'top µM':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        flag = "wt" if r["p53_functional"] else "mut"
        over = "  >top" if not r["ic50_in_tested_range"] else ""
        print(f"{r['drug']:<12}{r['cell_line']:<12}{flag:<5}{r['dataset']:<7}{r['drug_id']:>6}"
              f"{r['ic50_uM']:>10.3g}{r['auc']:>7.3f}{r['max_conc_uM']:>8.3g}{over}")
    print(f"\n{len(rows)} screens. '>top' means the screen never reached 50 % kill; "
          f"that IC50 is an extrapolation, use it as a lower bound only.")

    # Replicate agreement: the same (drug, line) screened more than once.
    groups: dict[tuple, list[float]] = {}
    for r in rows:
        if r["ic50_in_tested_range"]:
            groups.setdefault((r["drug"], r["cell_line"]), []).append(r["ic50_uM"])
    spreads = [(k, max(v) / min(v)) for k, v in groups.items() if len(v) > 1]
    if spreads:
        print("\nReplicate spread (max/min IC50 across screens, in-range only):")
        for (drug, line), s in sorted(spreads):
            print(f"  {drug:<12}{line:<12}{s:>6.1f}x")
        med = statistics.median(s for _, s in spreads)
        print(f"  median {med:.1f}x — no model should be gated tighter than the "
              f"data agrees with itself.")

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
