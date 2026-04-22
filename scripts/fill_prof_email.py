#!/usr/bin/env python3
"""Fill the professor-debrief email template with numbers from a
fep-report run — hydration (Milestone A / FreeSolv) OR binding
(Milestone B / streptavidin / EGFR).

Usage:
    python scripts/fill_prof_email.py run/fep/verdict/

The script auto-detects whether the report is hydration or binding
by reading the first-line markdown title (# Hydration FEP report ...
vs # Binding FEP report ...) and selects the matching email
template.

From the report.md it parses:
  Hydration-specific:
    - methane + acetamide sign-correctness (Milestone A critical pair)
  Binding-specific:
    - per-compound sign-correctness (every ΔG_bind must be < 0)
  Both:
    - MAE, RMSE, Pearson r, Spearman ρ, Kendall τ
    - GHMC mean + worst-window
    - overall verdict (PASS / FAIL / partial / inconclusive)
    - n/m compound count

Emits the filled email body on stdout (or --out file). The
parity.png + table.csv you attach are whatever's in the same dir.

Exit 0 if all fields filled; 2 if any <missing> remains so you
don't accidentally send placeholder markers to the prof.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


HYDRATION_TEMPLATE = """\
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
- Compounds completed:              {n_ok}/{n_total}
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
   `cellsim doctor` runs 42 install + benchmark checks, 54+ smoke
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


BINDING_TEMPLATE = """\
Subject: Milestone B — {subject_tag} FEP results

Hi {prof_name},

Milestone B Phase-2 ran on {bench_name}. Results below.

== What ran ==

- Hardware: {hardware_hint} (OpenMM {platform_hint} backend)
- Force field: AMBER14 (protein + tip3pfb water) + OpenFF Sage 2.1.0
  + AM1-BCC charges (ligand) via openmmforcefields'
  SMIRNOFFTemplateGenerator — no learned surrogate at any layer
- Method: double-decoupling (DDM) absolute binding ΔG
    ΔG_bind = −(ΔG_decouple_complex − ΔG_decouple_solvent)
              + ΔG_restraint_correction (Hamelberg-Gilson analytical)
- MD: openmmtools.alchemy + GHMC integrator, 1 fs timestep,
  11 lambda-windows per leg × 2 legs (complex + solvent)
- Estimator: pymbar 4.2.0 MBAR
- Per compound: 50 ps equilibration + 50 ps production per window
  ~= 2.2 ns simulated MD per compound
- Compounds: {n_total_compounds} ({bench_name} series)

== Headline numbers ==

- Overall verdict:                     {overall_verdict}
- Compounds completed:                 {n_ok}/{n_total}
- MAE vs published ΔG_bind:            {mae:>5} kcal/mol (gate <= 2.0)
- Pearson r:                           {pearson:>5}
- Spearman rho:                        {spearman:>5}
- Kendall tau:                         {kendall:>5}  (rank-correlation gate: >= 0.6)
- All compounds predicted as binders:  {binding_sign}
- GHMC acceptance:                     mean {ghmc_mean}, worst {ghmc_worst} (gate >= 70%)

== Parity figure ==

Attached parity.png — predicted vs experimental ΔG_bind, ±2.0 kcal/mol
gate band shaded, per-point error bars from MBAR. Compound labels
adjacent to each point.

== What this answers ==

This is the test that retires your 'physics-FEP vs Vina on kinases'
critique. The same {bench_name} chemical series where Vina's
empirical scoring gave Spearman −0.49 (anti-correlated with
experiment) — this FEP run reports Spearman {spearman} and Kendall
{kendall} against the same published reference ΔG values.

The absolute-ΔG MAE carries a Cheng-Prusoff offset (we used ΔG = RT
ln(IC50) since papers rarely report the [ATP]_Km pairs needed for the
correction). That offset is constant across the series, so the rank-
correlation (Kendall τ) is the load-bearing metric; absolute MAE is
informative but offset-limited.

== What this does NOT yet address ==

1. Protein-specific force-field transferability. ff14SB covers
   standard amino acids; non-standard residues / cofactors / metal
   sites need separate parametrisation.
2. Slow conformational changes. 50 ps/window doesn't capture
   rearrangements > 100 ns timescale (rare-event activation loops,
   large domain motions).
3. Campaign-2 cell-level numbers. Per your gate, those start after
   Milestone B clears.

== Proposed next ==

{proposed_next}

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
    """Hydration-specific: check sign of methane/acetamide rows."""
    pattern = rf"\|\s*{re.escape(name)}\s*\|[^|]*\|\s*([+-][\d.]+)\s*\|\s*([+-][\d.]+)\s*\|"
    m = re.search(pattern, text)
    if not m:
        return None
    expt, pred = m.groups()
    expt_v = float(expt)
    pred_v = float(pred)
    if abs(expt_v) < 0.3:
        return "PASS"
    same_sign = (expt_v >= 0) == (pred_v >= 0)
    return "PASS" if same_sign else "FAIL"


def _grab_binding_all_negative(text: str) -> str | None:
    """Binding-specific: every predicted ΔG_bind must be < 0.
    Scrape the per-compound table and check the 'pred' column
    for all rows. Scaffolded/failed rows are ignored — they get
    'scaffolded' or 'FAIL' in the pred cell and don't count."""
    rows = re.findall(
        r"\|\s*(\S+)\s*\|\s*`[^`]*`\s*\|\s*"
        r"[+-]?[\d.]+\s*\|\s*([+-][\d.]+)\s*\|",
        text)
    if not rows:
        return None
    non_binders = [n for n, pred in rows if float(pred) >= 0]
    if non_binders:
        return f"FAIL ({', '.join(non_binders)} predicted non-binder)"
    return "PASS"


def _grab_overall_verdict(text: str) -> str:
    """First-line header: '# Hydration FEP report — PASS' (or FAIL,
    inconclusive, partial)."""
    m = re.search(r"^#\s+\S.*?—\s*(.+)$", text, re.MULTILINE)
    return m.group(1).strip() if m else "(unknown)"


def _detect_yaml_kind(text: str) -> str:
    """From the markdown title '# Hydration FEP report' or
    '# Binding FEP report'. Default 'hydration'."""
    m = re.search(r"^#\s+(Hydration|Binding)\s+FEP report",
                  text, re.MULTILINE)
    if m:
        return m.group(1).lower()
    return "hydration"


def _grab_counts(text: str) -> tuple[str, str]:
    """From the 'compounds: N ok / M total' line, return
    (n_ok, n_total) as strings."""
    m = re.search(
        r"compounds:\s+(\d+)\s+ok\s*/\s*(\d+)\s+total", text)
    if m:
        return m.group(1), m.group(2)
    return "?", "?"


def _infer_bench_name(text: str) -> str:
    """For binding: try to pick a human-readable bench name from
    rows present. Fallback: generic 'binding'."""
    rows = [m.group(1) for m in re.finditer(
        r"^\|\s*(\S+)\s*\|\s*`[^`]*`\s*\|", text, re.MULTILINE)]
    row_names = {r.lower() for r in rows}
    # Streptavidin markers
    if {"biotin", "desthiobiotin"} & row_names:
        return "streptavidin"
    if {"erlotinib", "gefitinib", "ag-1478", "lapatinib"} & row_names:
        return "EGFR kinase"
    return "binding"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Fill the prof-debrief email from a fep-report "
                    "directory's report.md. Auto-detects hydration "
                    "vs binding from the report title.")
    ap.add_argument(
        "report_dir",
        help="path to the run/fep/verdict/ directory (or any dir "
             "containing report.md)")
    ap.add_argument(
        "--prof-name", default="[Prof]",
        help="name to address (default '[Prof]' — fill in by hand)")
    ap.add_argument(
        "--hardware", default="rented GPU",
        help="binding-template hardware hint (default 'rented GPU'; "
             "override with 'Apple M5 Max (40-core GPU)' etc.)")
    ap.add_argument(
        "--platform", default="CUDA",
        help="binding-template platform hint (CUDA | Metal | OpenCL)")
    ap.add_argument(
        "--next-step", default=None,
        help="binding-template 'Proposed next' paragraph override")
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
    kind = _detect_yaml_kind(md)

    overall = _grab_overall_verdict(md)
    n_ok, n_total = _grab_counts(md)
    mae = _grab_float(md, "MAE") or "<missing>"
    pearson = _grab_float(md, "Pearson r") or "<missing>"
    spearman = _grab_float(md, "Spearman ρ") or "<missing>"
    kendall = _grab_float(md, "Kendall τ") or "<missing>"
    ghmc_mean, ghmc_worst = _grab_ghmc(md)
    ghmc_mean = ghmc_mean or "<missing>"
    ghmc_worst = ghmc_worst or "<missing>"

    if kind == "binding":
        bench_name = _infer_bench_name(md)
        # Subject line reads cleaner with a named bench vs the
        # fallback. 'streptavidin' → 'Subject: Milestone B —
        # streptavidin binding FEP', 'binding' → 'Subject: Milestone
        # B — binding FEP'.
        if bench_name == "binding":
            subject_tag = "binding"
        else:
            subject_tag = f"{bench_name} binding"
        binding_sign = _grab_binding_all_negative(md) or "<missing>"
        proposed = args.next_step or (
            "If Kendall τ >= 0.6 AND all compounds predicted as "
            "binders: the EGFR rank-order rescue claim is evidenced "
            "— happy to discuss Campaign-2 sequencing next.\n\n"
            "If τ < 0.6 or any non-binder prediction: we'd identify "
            "the problem compound(s) and extend sampling / inspect "
            "pose before trusting downstream.")
        filled = BINDING_TEMPLATE.format(
            prof_name=args.prof_name,
            bench_name=bench_name,
            subject_tag=subject_tag,
            hardware_hint=args.hardware,
            platform_hint=args.platform,
            n_total_compounds=n_total,
            n_ok=n_ok,
            n_total=n_total,
            overall_verdict=overall,
            mae=mae,
            pearson=pearson,
            spearman=spearman,
            kendall=kendall,
            binding_sign=binding_sign,
            ghmc_mean=ghmc_mean,
            ghmc_worst=ghmc_worst,
            proposed_next=proposed,
        )
    else:
        # Hydration template (Milestone A)
        methane_sign = _grab_compound_sign(md, "methane") or "<missing>"
        acetamide_sign = _grab_compound_sign(md, "acetamide") or "<missing>"
        filled = HYDRATION_TEMPLATE.format(
            prof_name=args.prof_name,
            overall_verdict=overall,
            n_ok=n_ok,
            n_total=n_total,
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

    if "<missing>" in filled:
        print("\nfill_prof_email: WARNING — some fields are missing "
              "(see <missing> markers above). Inspect report.md "
              "manually.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
