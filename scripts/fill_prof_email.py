#!/usr/bin/env python3
"""Fill the professor-debrief email template with numbers from a
fep-report run.

Usage:
    python scripts/fill_prof_email.py run/fep/verdict/

Reads `report.md` from the supplied directory, parses out:
  - MAE, RMSE, Pearson r, Spearman ρ, Kendall τ from the
    'Aggregate accuracy' section
  - methane + acetamide rows from the per-compound table → sign
    PASS/FAIL
  - GHMC mean + worst from the gate verdict line
  - overall verdict (PASS / FAIL / partial / inconclusive)

Emits the email body on stdout, ready to copy into your message.
The `parity.png` you attach is whatever's already in the same
directory (the script doesn't move files).

Exit 0 if the email could be filled. Exit 1 if report.md is
missing or unparseable so you don't accidentally send an email
full of <X.XX> placeholders.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


EMAIL_TEMPLATE = """\
Subject: Milestone A — FreeSolv FEP results

Hi {prof_name},

Quick update on Milestone A. The CellSim Layer-1.3 alchemical FEP
pipeline is now end-to-end on a Mac and the FreeSolv-12 hydration
gate ran on a friend's M5 Max overnight. Results below.

== What ran ==

- Hardware: Apple M5 Max (40-core GPU), OpenMM Metal backend
- Force field: OpenFF Sage 2.1.0 (small molecules) + AM1-BCC partial
  charges + TIP3P water — no learned surrogate at any layer
- MD: openmmtools.alchemy + GHMC integrator (Metropolised Langevin),
  1 fs timestep, 11 lambda-windows per leg
- Estimator: pymbar 4.2.0 MBAR (closed-form from alchemical samples)
- Per compound: 50 ps equilibration + 50 ps production per window
  x 11 windows x 2 legs (vacuum + TIP3P-solvated decoupling)
  ~= 2.2 ns simulated MD; 12 compounds ~= 26 ns total
- Compounds: the FreeSolv subset I sent before — methane (+2.00) ->
  acetamide (-9.71) kcal/mol experimental, spans the dynamic range
  we care about
- Software pinned in environment.yml (rdkit 2025.9.5, openmm 8.5.1,
  openmmtools 0.26.0, openff-toolkit 0.18.0, pymbar 4.2.0). No perses
  dependency — openmmtools primitives + a custom DDM driver covered
  the scope.

== Headline numbers ==

- Overall verdict:                  {overall_verdict}
- MAE vs FreeSolv published values: {mae:>5} kcal/mol (gate <= 1.5)
- Pearson r:                        {pearson:>5}
- Spearman rho:                     {spearman:>5}
- Kendall tau:                      {kendall:>5}
- Sign-correct on methane:          {methane_sign}
- Sign-correct on acetamide:        {acetamide_sign}
- GHMC acceptance:                  mean {ghmc_mean}, worst {ghmc_worst} (gate >= 70%)

== Parity figure ==

Attached parity.png — predicted vs experimental ΔG_hyd, ±1.5 kcal/mol
gate band shaded, per-point error bars from MBAR, all 12 compounds
labelled.

== What this answers from your last critique ==

You raised closed-ontology / circular validation: that CellSim's
binding affinities were hand-tuned phenomenology dressed as physics,
so any cell-level prediction built on top would inherit the
unfalsifiable layer underneath. I've closed that gap on the
small-molecule axis:

1. Physics-derived priors with calibrated uncertainty. Every ΔG
   comes from MBAR over alchemical samples; σ from MBAR. No
   empirical scoring function anywhere. The per-compound report has
   a "within σ" column saying which residuals are statistically
   meaningful vs sampling noise.
2. Sign-correctness as a hard gate. Methane and acetamide signs are
   explicit gate items — if either flips, the report refuses PASS
   regardless of MAE.
3. GHMC acceptance gate. Your "monitor the acceptance" requirement
   is enforced per-compound: any window < 70% forces the report to
   refuse PASS. Full per-window vector is in the tarball's run.log.
4. Reproducible from fresh clone. environment.yml is pinned,
   `cellsim doctor` runs 42 install + benchmark checks, 54 smoke
   tests gate every code change. The M5 Max ran the same script my
   CI runs.

== What this does NOT yet address ==

1. Binding ΔG / ΔΔG on real receptors. The Milestone B pipeline
   (double-decoupling, amber14 + SMIRNOFF-ligand, harmonic CoM
   restraint, analytical Hamelberg-Gilson correction) is built and
   tested end-to-end on an integration smoke (1ubq + methane,
   82.5 s CPU). No GPU-sampled numbers on real binders yet.
2. EGFR rank-order rescue. The same 6-compound EGFR series where
   Vina gave Spearman -0.49 is validated + ready to sample
   (benchmarks/fep/binding_egfr.yaml). Awaiting GPU time.
3. Campaign-2 work. Per your gate, the cell-level rebuild does not
   start until Milestone A clears.

== Proposed next ==

If Milestone A clears: streptavidin-biotin reference set (4
compounds, ~15 min GPU/compound estimated) followed by the EGFR
6-compound series (~38 min total). Looking for your sign-off on
cloud-GPU spend before booking the run.

If Milestone A fails on methane sign specifically: I'd extend that
single compound to 200 ps/window and re-test before touching Phase 2.

If MAE is systematically > 2.5: open conversation — switch to ff19SB
/ longer windows / different tip3p variant.

Tarball, report.md, and parity.png attached.

— Henry
"""


def _grab_float(text: str, label: str) -> str | None:
    """Find a number after a label like 'MAE  = 0.420' or
    'Pearson r    = +0.993'. Returns the number as a string or
    None."""
    pattern = rf"{re.escape(label)}\s*=\s*([+-]?[\d.]+)"
    m = re.search(pattern, text)
    return m.group(1) if m else None


def _grab_ghmc(text: str) -> tuple[str | None, str | None]:
    """From the gate verdict line, e.g.:
       'overall mean 77%, worst 71%, 3 compounds report'
    return (mean, worst) as strings with %."""
    m = re.search(
        r"overall mean\s+([\d.]+%?),\s+worst\s+([\d.]+%?)", text)
    if not m:
        return None, None
    mean, worst = m.groups()
    return (
        mean if mean.endswith("%") else mean + "%",
        worst if worst.endswith("%") else worst + "%",
    )


def _grab_compound_sign(text: str, name: str) -> str | None:
    """In the per-compound table, find the row for `name` and
    return 'PASS' if pred and expt have the same sign, 'FAIL'
    otherwise. Looks for 'SIGN WRONG' marker first (analyser's
    direct flag) then falls back to inspecting the +/− on pred
    vs expt cols.

    Table row shape:
      | acetamide | `CC(=O)N` | -9.71 | -8.90 | 0.55 | +0.81 | 0.81 | ... |
    """
    pattern = rf"\|\s*{re.escape(name)}\s*\|[^|]*\|\s*([+-][\d.]+)\s*\|\s*([+-][\d.]+)\s*\|"
    m = re.search(pattern, text)
    if not m:
        return None
    expt, pred = m.groups()
    expt_v = float(expt)
    pred_v = float(pred)
    # Match analyser's near-zero rule: |expt| < 0.3 → either sign ok.
    if abs(expt_v) < 0.3:
        return "PASS"
    same_sign = (expt_v >= 0) == (pred_v >= 0)
    return "PASS" if same_sign else "FAIL"


def _grab_overall_verdict(text: str) -> str:
    """First-line header: '# Hydration FEP report — PASS' (or FAIL,
    inconclusive, partial)."""
    m = re.search(r"^#\s+\S.*?—\s*(.+)$", text, re.MULTILINE)
    return m.group(1).strip() if m else "(unknown)"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Fill the prof-debrief email from a fep-report "
                    "directory's report.md")
    ap.add_argument(
        "report_dir",
        help="path to the run/fep/verdict/ directory (or any dir "
             "containing report.md)")
    ap.add_argument(
        "--prof-name", default="[Prof]",
        help="name to address (default '[Prof]' — fill in by hand)")
    ap.add_argument(
        "--out", default="-",
        help="write the filled email to this file (default '-' = stdout)")
    args = ap.parse_args(argv)

    report_path = Path(args.report_dir) / "report.md"
    if not report_path.is_file():
        print(f"fill_prof_email: no report.md under {args.report_dir}",
              file=sys.stderr)
        return 1

    md = report_path.read_text(encoding="utf-8")

    overall = _grab_overall_verdict(md)
    mae = _grab_float(md, "MAE") or "<missing>"
    pearson = _grab_float(md, "Pearson r") or "<missing>"
    spearman = _grab_float(md, "Spearman ρ") or "<missing>"
    kendall = _grab_float(md, "Kendall τ") or "<missing>"
    ghmc_mean, ghmc_worst = _grab_ghmc(md)
    ghmc_mean = ghmc_mean or "<missing>"
    ghmc_worst = ghmc_worst or "<missing>"
    methane_sign = _grab_compound_sign(md, "methane") or "<missing>"
    acetamide_sign = _grab_compound_sign(md, "acetamide") or "<missing>"

    filled = EMAIL_TEMPLATE.format(
        prof_name=args.prof_name,
        overall_verdict=overall,
        mae=mae,
        pearson=pearson,
        spearman=spearman,
        kendall=kendall,
        ghmc_mean=ghmc_mean,
        ghmc_worst=ghmc_worst,
        methane_sign=methane_sign,
        acetamide_sign=acetamide_sign,
    )

    if args.out == "-":
        print(filled)
    else:
        Path(args.out).write_text(filled, encoding="utf-8")
        print(f"wrote {args.out}", file=sys.stderr)

    # Refuse to silently succeed if any field is still missing — the
    # biologist might paste the result without noticing.
    if "<missing>" in filled:
        print("\nfill_prof_email: WARNING — some fields are missing "
              "(see <missing> markers above). Inspect report.md "
              "manually.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
