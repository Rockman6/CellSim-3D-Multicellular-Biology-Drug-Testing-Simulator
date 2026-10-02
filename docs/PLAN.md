# CellSim plan: from two half-products to one finished one

*Written October 2026 for everyone working on the repository. Supersedes
`docs/archive/campaign1/GOAL.md`, `MISSION.md` and `ROADMAP.md`.*

## Why the project needed a plan

Until October 2026 the repository held two different products stacked on
top of each other: the 2026 C++/Metal 3D cell simulator the project is
named for, and a Python "non-AI drug-discovery triage" stack built after
the April restart. Each stopped at "the pipeline runs" rather than "the
answer is right", and the original five-campaign roadmap was sized for a
large team over four to seven years. One person executed about five per
cent of each of seven layers in parallel. The result read as a half
product because it was two quarter products.

What was verified before writing this plan (details in
`VALIDATION.md` and the archive):

- Docking pose numbers were computed on scrambled molecules for months.
  The corrected 15-cocrystal blind set gives 87 % top-3 pose recovery.
- The hydration FEP "pass" (MAE 1.42 kcal/mol on 10 of 12) hides a
  size-dependent bias: the solvated box runs NVT with no barostat.
  Binding FEP never produced a number on a real binder.
- The legacy simulator builds headless and passes its 8 benchmarks at
  two seeds, but cells never execute death under cisplatin and the p53
  period is 3.0 h against 5.5 h in the literature.
- `cellsim/cell` is 20 sound textbook PK/PD modules with no shared state and
  demo inputs typed by hand.
- 204 000 words of documentation for 25 000 lines of Python.

## Definition of done: v2.0

An interactive dish. A user picks a cell line and a drug (with a measured
or computed affinity and its ADME properties), sets a dose schedule,
watches the colony respond, and gets a dose-response curve with a
calibrated error bar that has been checked against public cell-line data
(GDSC2 or NCI-60). The molecular layer is an optional input plugin.
Everything that does not serve that is archived.

The mechanistic core stays physics-first and interpretable. The old
"non-AI at every layer" rule no longer applies to *inputs*: binding
affinities come from ChEMBL or BindingDB with provenance first, from
docking with its reliability flag second. The non-circular validation the
project was asked for comes from external cell-line data, not from
re-implementing free-energy codes.

## Phases

### Phase 0, done October 2026

- PRs #10, #11, #12 merged; their 29 test files gated in CI.
- CI green; the red badge since May was a noise-boundary gate on
  scrambled poses.
- One README, one identity, GitHub metadata and license detection fixed.
- Campaign-era docs archived under `docs/archive/campaign1/`; handoff
  tooling, run logs and retired-UI screenshots removed; empty scaffolds
  deleted.

### Phase 1, the validated slice — CLOSED October 2026, with a finding

No UI work until this exists.

1. *Done.* `cellsim/cell/engine.py`: one state array, one RK4
   integrator, one `simulate`. A population is weighted representative
   cells, so a whole dose-response is one run.
2. *Done.* Ported from `OLD/`: the CDK/cyclin cycle, ATM → p53 ⇄ MDM2
   (3-variable delay oscillator) → p21 and PUMA → BAX → MOMP. The p53
   period is now 4.9 h (published 5.5, prototype 3.0) and caspase-3
   execution actually kills cells. Ten phenotype gates pass
   (`tests/cell/test_engine_smoke.py`).
3. Drug input: PK/disposition from `cellsim/cell` (permeation, pH trapping,
   efflux, binding sink) feeding occupancy at the target, then fate.
   Affinity and ADME from a small curated table with provenance
   (ChEMBL / BindingDB / literature), docking as a flagged fallback.
4. Validation harness against GDSC fitted dose-response.

   *Revised October 2026 after pulling the data* (`scripts/gdsc_reference.py`
   reads GDSC1 and GDSC2 release 8.4; 30 screens extracted to
   `benchmarks/cell/gdsc_reference.csv`). The original plan named
   cisplatin, doxorubicin and paclitaxel with a "within 3×" gate. What
   the public data actually supports:

   - **The gate cannot be tighter than the data's own replicates.** GDSC1
     screened doxorubicin twice per line, and the two in-range IC50s of
     the same line disagree by 1.5× to 9.3× (median 5.6×). So the
     quantitative gate is: predicted IC50 inside the span of the
     replicate screens, widened to at least 3× either side of the
     geometric mean when there is only one screen; plus the correct
     ordering of the most and least sensitive line.
   - **Paclitaxel (GDSC2) and doxorubicin (GDSC1) are the quantitative
     targets.** Paclitaxel 0.011–0.089 µM, one screen per line, all in
     range. Doxorubicin 0.010–0.37 µM, nine of ten screens in range.
   - **Cisplatin cannot be a quantitative target.** Only 2 of its 15
     screens reached 50 % kill below their top dose (A549 9.8 µM and
     HeLa 6.9 µM, GDSC1 at 10 µM top). The rest are extrapolations up to
     122 µM. Gate it on direction and rank only: weak, high-µM, never
     sub-µM.
   - **Camptothecin is not a usable substitute**: four of five lines sit
     above its 0.1 µM top dose in GDSC2.
   - **The p53 split is not clean** and must not be assumed. Cisplatin is
     less potent in p53-mutant HT-29 than in wild-type A549, as expected,
     but p53-mutant MDA-MB-231 is more sensitive than wild-type MCF7.
     Report the p53 effect as a measured outcome, never as a gate.

   Fit the one potency gain per drug on the p53 wild-type lines; predict
   the mutant lines with nothing re-tuned. Publish the table including
   every miss.

   *First cycle done* (`scripts/validate_gdsc.py`, results in
   `docs/VALIDATION.md`). Nine of ten held-out IC50s land inside the
   GDSC span, but a constant-IC50 null does equally well there, and on
   log error the engine beats the null only for paclitaxel. So the
   potency scale is right and cross-line discrimination is not yet.
   **The Phase-1 exit gate is therefore: beat the null's log10 RMSE on
   held-out lines for at least two of the three drugs.**

   *Status, October 2026: NOT met, and now measured properly.* On ten
   lines (five TP53 wild-type, five mutant, fitted on two) the engine is
   closer than a constant-IC50 null on 10 held-out lines and further on
   11, sign-test p = 1.0. An earlier five-line run appeared to pass, but
   the gate then compared two RMSE numbers and the winning margin was
   0.45 vs 0.46 — noise. The gate is now a paired per-line sign test.

   What IS established: each drug's potency scale, from one fitted
   constant, with 22 of 25 held-out predictions inside GDSC's replicate
   span against the null's 18.

   What the measurement says is missing, now pinned down: the engine has
   no channel through which lines can differ. Sweeping each per-line
   field across its full observed range moves the predicted IC50 by
   1.0x for doubling time, 1.1x for the Bcl-2 reserve, and 7-9x only
   for DNA repair rate (and not at all for paclitaxel, whose IC50 is
   set by a drug-intrinsic arrest threshold). Lines differ by 8-51x in
   reality. So the engine can produce at most two answers per drug, and
   tying a constant is arithmetic, not bad luck.

   Also established: DepMap 24Q4 expression for all ten lines and
   sixteen candidate genes correlates with sensitivity no better than
   chance once corrected for the search (family-wise p 0.25-0.86), and
   DepMap's CRISPR growth rate does not match curated doubling times
   (rho -0.15). Building a mapping now would fit noise.

   At ten lines no mapping could be validated. Repeating the test on
   every GDSC line DepMap has expression for (156-683 per drug), with
   genes pre-specified by mechanism and Bonferroni-corrected, real
   signal does appear: Bcl-xL (BCL2L1) and P-glycoprotein (ABCB1) for
   doxorubicin and paclitaxel, plus class III beta-tubulin (TUBB3) for
   paclitaxel. Cisplatin gives nothing, and notably the repair genes do
   not predict it (ERCC1 rho -0.03).

   Putting that beside the leverage measurement explains the whole
   failure. The engine's only strong channel is DNA repair, where the
   data has no signal. The strongest real signal is the apoptotic
   set-point, where the engine HAS the field (CellLine.bcl2_level) but
   its MOMP switch is so sharp that an 8-fold change moves the IC50 by
   1.1x, against the ~2x the data implies. The second strongest is
   efflux, which the engine does not model at all.

   Both channels were then built (apoptotic buffer, efflux) and
   neither moved the gate, which led to the measurement that settles
   it: fitting ALL eleven markers to IC50 by least squares across every
   GDSC line DepMap covers gives best-case R2 of 0.09 / 0.15 / 0.20.
   80-92 % of line-to-line variation is not in these markers, so a
   mechanistic model restricted to them cannot pass the paired test
   however faithfully it is built.

   DECISION: per-line IC50 prediction is not a goal this engine can
   reach with available inputs, and is dropped as the headline metric.
   Phase 2 re-aims at within-line dynamics -- schedule dependence,
   wash-out recovery, resistance under selection, combination timing --
   where the mechanism is load-bearing and no per-line marker is
   required. See docs/VALIDATION.md for the numbers.

5. *Done.* Headless state stream: `cellsim/cell/stream.py` runs one
   population under a piecewise-constant dosing schedule (wash-outs
   included) and emits per-tick JSON Lines, schema
   `cellsim.cell.stream/v1` (`cellsim cell-stream`). Any UI is a
   renderer of this stream.
6. A Colab notebook that reproduces the validation figure from a clean
   environment.
7. Packaging (done in Phase 0): `pyproject.toml`, the `cellsim` console
   script, pytest collection; the `src` package is now `cellsim`.

### Phase 2, within-line dynamics and the dish (now current)

**Re-aimed October 2026.** Phase 1 measured that per-line IC50
prediction is unreachable from canonical markers (R² 0.09–0.20 across
156–683 lines). Phase 2 therefore targets the questions mechanism
answers well and statistics answer badly — all of them *within* one
line, so no per-line marker is needed:

- **Schedule dependence.** Does a long low exposure differ from a short
  high one at matched AUC? It should, and differently per drug: a
  tubulin binder kills only cells that reach mitosis during exposure,
  so it is exposure-time limited in a way a DNA-damage agent is not.
- **Wash-out and recovery.** Already shown: a 6 h pulse regrows ×4.0
  where the same dose held continuously eradicates the colony.
- **Resistance under selection**, and whether intermittent dosing
  delays relapse relative to continuous dosing at matched exposure.
- **Combination timing**, where a cytostatic given first antagonises a
  phase-specific partner.

These are validated against published *phenotypes* rather than a
ten-line IC50 table, which is the kind of evidence that does not run
into the ceiling above.

- Space: lattice or off-lattice agents, diffusion of oxygen, glucose
  and drug, a vessel source, contact inhibition. Reuse
  `cellsim/cell/tissue.py` and `agents.py`.
- Validate on spheroid growth curves from the literature and on the
  Cell Tracking Challenge HeLa counts already in the repository.
- New UI, web first (three.js or WebGPU) so it runs anywhere, reading
  the engine's state stream; per-cell pathway readouts, dose schedules,
  uncertainty bands, CSV export.
- Three outside users: a course, an iGEM team, a lab.

### Phase 3, credibility (target October 2027)

- Ten drugs × five lines; combination schedules; resistance evolution.
- Uncertainty calibrated per target class, with coverage reported.
- JOSS paper for the software, a short validation preprint.
- conda-forge package, documentation site, hosted demo.

## Keep / kill

**Keep.** Novak-Tyson and p53 biology from `OLD/` (as source material),
the six headless validators as regression tests, `cellsim/cell`,
`cellsim/bridge`, ADMET descriptors, the cache, docking as a plugin with its
reliability table, the docking blind set and calibration bundles.

**Kill or freeze.** Hand-rolled binding FEP as a product path (frozen as
experimental; adopt OpenFE later only if GPU time appears), the Martini
scaffold (deleted), DFT rescoring and the MD layer as user-facing
features, the Metal renderer, the campaign framing, and all per-person
handoff tooling.

## Success metrics

- CI green on every merge; `pip install` under five minutes.
- A published validation table against GDSC with named misses.
- Three outside users by month six, one paper by month twelve.

## Decisions taken

- One product: the validated cell and drug-response simulator.
- New UI is web-based and reads a headless state stream.
- Affinities come from measured data by default; docking is a flagged
  fallback; in-house FEP is not required for validation.
