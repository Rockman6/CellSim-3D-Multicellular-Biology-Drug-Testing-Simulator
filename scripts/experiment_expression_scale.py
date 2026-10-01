#!/usr/bin/env python3
"""Does per-line expression predict drug sensitivity, at a sample size that can tell?

The Phase-1 diagnosis (`experiment_leverage.py`) says the engine needs a
per-line channel with leverage, and the obvious candidates are repair
capacity and efflux. At ten cell lines nothing could be validated:
correcting for having searched sixteen genes, no correlation beat chance.

This repeats the test on every GDSC line that DepMap has expression for
— 156 to 1204 lines depending on the drug, rather than 6 to 10 — so a
real effect of the size the engine would need becomes visible.

Genes are **pre-specified by mechanism**, declared below before any
correlation is computed, so the result is not the best of a search:

* efflux: ABCB1 (P-gp) for doxorubicin and paclitaxel, which are its
  substrates; ABCC1, ABCG2.
* cisplatin handling: SLC31A1 (CTR1 uptake), ERCC1/ERCC2/XRCC1 (repair
  of Pt-DNA adducts), GSTP1 (conjugation).
* apoptotic set-point: BCL2, BCL2L1, MCL1 against BAX, BAK1, BBC3.
* drug targets: TOP2A (doxorubicin), TUBB3 (class III beta-tubulin, a
  documented taxane-resistance marker).

Data (not redistributed; both are large public files):

    mkdir -p data/gdsc data/depmap
    # GDSC fitted dose response, release 8.4
    curl -L -o data/gdsc/GDSC1_fitted_dose_response_24Jul22.csv \\
      https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/GDSC1_fitted_dose_response_24Jul22.csv
    curl -L -o data/gdsc/GDSC2_fitted_dose_response_24Jul22.csv \\
      https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/GDSC2_fitted_dose_response_24Jul22.csv
    # DepMap 24Q4: Model.csv and the protein-coding expression matrix
    curl -L -o data/depmap/Model.csv https://ndownloader.figshare.com/files/51065297
    curl -L -o data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv \\
      https://ndownloader.figshare.com/files/51065489

Usage:
    python scripts/experiment_expression_scale.py
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

DEPMAP = REPO_ROOT / "data" / "depmap"
GDSC = REPO_ROOT / "data" / "gdsc"
OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "expression_scale_experiment.json"

# drug -> (GDSC release file, pre-specified genes with the direction expected)
# direction +1 = higher expression should mean MORE resistant (higher IC50)
DRUGS = {
    "Cisplatin": ("GDSC2", {
        "SLC31A1": -1,   # CTR1 uptake -> more drug in -> more sensitive
        "ERCC1": +1, "ERCC2": +1, "XRCC1": +1,   # repair -> resistant
        "GSTP1": +1,     # conjugation/detox -> resistant
        "BCL2L1": +1, "MCL1": +1, "BCL2": +1,    # anti-apoptotic -> resistant
        "BAX": -1, "BAK1": -1, "BBC3": -1,       # pro-apoptotic -> sensitive
    }),
    "Doxorubicin": ("GDSC1", {
        "ABCB1": +1, "ABCC1": +1, "ABCG2": +1,   # efflux -> resistant
        "TOP2A": -1,                              # more target -> more sensitive
        "BCL2L1": +1, "MCL1": +1, "BCL2": +1,
        "BAX": -1, "BAK1": -1, "BBC3": -1,
    }),
    "Paclitaxel": ("GDSC2", {
        "ABCB1": +1, "ABCC1": +1, "ABCG2": +1,
        "TUBB3": +1,                              # class III beta-tubulin -> taxane resistance
        "BCL2L1": +1, "MCL1": +1, "BCL2": +1,
        "BAX": -1, "BAK1": -1, "BBC3": -1,
    }),
}


def _strip(n: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", (n or "").upper())


def spearman(a: list[float], b: list[float]) -> float:
    def ranks(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):                       # average ties
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = ranks(a), ranks(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
    return num / den if den else float("nan")


def p_from_rho(rho: float, n: int) -> float:
    """Two-sided p via the t approximation, adequate for n > 30."""
    if n < 5 or not math.isfinite(rho) or abs(rho) >= 1:
        return float("nan")
    t = abs(rho) * math.sqrt((n - 2) / (1 - rho * rho))
    df = n - 2
    # survival function of |t| via the incomplete beta, series-free approximation
    x = df / (df + t * t)
    # regularised incomplete beta I_x(df/2, 1/2) == 2*P(T>|t|)
    a, b = df / 2.0, 0.5
    # continued fraction (Lentz) for I_x(a,b)
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1 - x) - lbeta) / a
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 300):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        d = 1e-30 if abs(d) < 1e-30 else d
        d = 1.0 / d
        c = 1.0 + num / c
        c = 1e-30 if abs(c) < 1e-30 else c
        f *= c * d
        if abs(1.0 - c * d) < 1e-10:
            break
    return min(1.0, front * (f - 1.0))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args(argv)

    for f in (DEPMAP / "Model.csv", DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"):
        if not f.exists():
            print(f"error: {f} missing — see this script's docstring.", file=sys.stderr)
            return 2

    name2mid = {}
    for r in csv.DictReader((DEPMAP / "Model.csv").open(newline="")):
        s = _strip(r.get("StrippedCellLineName"))
        if s:
            name2mid[s] = r["ModelID"]

    genes = sorted({g for _, gs in DRUGS.values() for g in gs})
    print(f"pre-specified genes ({len(genes)}): {', '.join(genes)}")

    # Stream the expression matrix once, keeping only the needed columns.
    expr: dict[str, dict[str, float]] = {}
    with (DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv").open(newline="") as f:
        rdr = csv.reader(f)
        hdr = next(rdr)
        col = {h.split(" (")[0]: i for i, h in enumerate(hdr) if h.split(" (")[0] in genes}
        missing = set(genes) - set(col)
        if missing:
            print(f"  not in matrix, dropped: {sorted(missing)}")
        for parts in rdr:
            expr[parts[0]] = {g: float(parts[i]) for g, i in col.items() if parts[i]}
    print(f"expression loaded for {len(expr)} models\n")

    report: dict = {"genes": genes, "drugs": {}}
    for drug, (release, spec) in DRUGS.items():
        src = GDSC / f"{release}_fitted_dose_response_24Jul22.csv"
        if not src.exists():
            print(f"{drug}: {src.name} missing, skipped")
            continue
        ic: dict[str, list[float]] = {}
        for r in csv.DictReader(src.open(newline="")):
            if r["DRUG_NAME"] != drug:
                continue
            v = math.exp(float(r["LN_IC50"]))
            if v > float(r["MAX_CONC"]):      # never reached 50 % kill
                continue
            mid = name2mid.get(_strip(r["CELL_LINE_NAME"]))
            if mid and mid in expr:
                ic.setdefault(mid, []).append(v)
        lines = sorted(ic)
        y = [math.log10(math.exp(sum(map(math.log, ic[m])) / len(ic[m]))) for m in lines]
        print(f"== {drug} ({release}): {len(lines)} cell lines with in-range IC50 + expression")
        if len(lines) < 20:
            print("   too few lines, skipped\n")
            continue
        n_tests = len(spec)
        rows = []
        for g, direction in sorted(spec.items()):
            pairs = [(expr[m][g], yy) for m, yy in zip(lines, y) if g in expr[m]]
            if len(pairs) < 20:
                continue
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            rho = spearman(xs, ys)
            p = p_from_rho(rho, len(pairs))
            p_adj = min(1.0, p * n_tests)          # Bonferroni over this drug's panel
            agrees = (rho > 0) == (direction > 0)
            rows.append({"gene": g, "expected_sign": direction, "n": len(pairs),
                         "rho": rho, "p": p, "p_bonferroni": p_adj,
                         "direction_as_expected": agrees,
                         "significant": p_adj < 0.05 and agrees})
        rows.sort(key=lambda r: -abs(r["rho"]))
        print(f"   {'gene':<9}{'exp':>5}{'rho':>8}{'p(Bonf)':>11}  verdict")
        for r in rows:
            mark = ("SIGNIFICANT, expected direction" if r["significant"]
                    else "wrong direction" if not r["direction_as_expected"]
                    else "n.s.")
            print(f"   {r['gene']:<9}{'+' if r['expected_sign']>0 else '-':>5}"
                  f"{r['rho']:>8.3f}{r['p_bonferroni']:>11.2e}  {mark}")
        hits = [r["gene"] for r in rows if r["significant"]]
        print(f"   -> {len(hits)}/{len(rows)} survive Bonferroni in the expected direction"
              f"{': ' + ', '.join(hits) if hits else ''}\n")
        report["drugs"][drug] = {"release": release, "n_lines": len(lines), "genes": rows,
                                 "significant": hits}

    if not a.no_write and report["drugs"]:
        OUT_JSON.write_text(json.dumps(report, indent=2, default=float) + "\n")
        print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
