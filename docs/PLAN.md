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
   constant, with 20 of 25 held-out predictions inside GDSC's replicate
   span against the null's 18 (22 of 25 before the buffer and efflux
   channels).

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
   Phase 3 confirmed this on ten drugs (32 wins to 41 losses over 91
   held-out predictions) and found the exception it does NOT cover:
   where the drug's target is itself part of the modelled mechanism —
   an MDM2 inhibitor acting through p53 — the engine calls all eleven
   lines correctly against a constant's five, HeLa included, whose
   wild-type TP53 would mislead a genotype-only model.
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

*Status, October 2026.* All four within-line targets are done and were
re-measured over seeds before anything was built on them (one claim
retracted, one qualified; `docs/VALIDATION.md`, "Phase 2: audited").

- *Done.* Space: `cellsim/cell/dish.py` puts the engine's cells on a
  lattice (monolayer or spheroid) with oxygen and drug diffusing from
  the medium, contact inhibition, hypoxic quiescence and necrosis. Built
  on the engine rather than on `agents.py`, whose closed-form fate rates
  are not the validated biology. Glucose and a vessel source are not in
  yet.
- *Done, with misses on record.* Validated on the Cell Tracking Challenge
  HeLa movie (colony size within 3-4 % at 46 h; the single-cell
  cycle-time split is a miss) and calibrated on DLD-1 spheroids (growth
  speed fixes the one spatial constant; the necrotic core is too large,
  a miss).
- *Done.* Web UI: `web/viewer/` renders both stream schemas in a
  browser — 3-D cells, colour by state, cut-open spheroids, time
  scrubbing, charts, CSV export. Uncertainty bands wait for Phase 3's
  calibrated intervals.
- Three outside users: a course, an iGEM team, a lab. Needs people, not
  code; the hosted demo (Phase 3) is what makes it possible.

### Phase 3, credibility (target October 2027)

- Ten drugs x five lines; combination schedules; resistance evolution.
- Uncertainty calibrated per target class, with coverage reported.
- JOSS paper for the software, a short validation preprint.
- conda-forge package, documentation site, hosted demo.

*Status, October 2026 — most of it done, honestly itemised.*

| Item | State |
|---|---|
| Ten drugs x ten lines | **Done.** 91 held-out predictions. Phase-1 conclusion confirmed; the MDM2 exception found (11/11 lines, null 5/11). |
| Combination schedules | **Done** in Phase 2, re-measured over seeds, and exposed through `cellsim.api.combination`. |
| Resistance evolution | **Done.** Heritable change at division, with the homogeneous-start control that makes the claim testable. |
| Uncertainty calibrated, coverage reported | **Done.** Leave-one-drug-out; within 15 points of nominal at 68/90/95 %. |
| JOSS paper | **Parked, not in progress.** A draft sits in `paper/` from when this plan listed it; the goal is a product, not a publication. |
| Validation preprint | **Dropped.** `docs/VALIDATION.md` serves the purpose that matters — making the tool trustworthy to someone using it. |
| conda-forge | **Recipe prepared** (`packaging/conda-forge/`), not submitted: the sdist must be on PyPI first, and submission is the owner's call. |
| Documentation site | **Served by** `web/index.html` plus `docs/cell_tutorial.md`; no separate site generator, which would be machinery for its own sake at this size. |
| Hosted demo | **Workflow ready** (`.github/workflows/pages.yml`), which simulates the demo runs at deploy time. Needs Pages enabled in repository settings once; its steps were verified locally against a pip-installed package. |

The two mechanisms added for the panel were S-phase-restricted damage
(TOP1 poisons and antimetabolites) and MDM2 inhibition, with the HPV E6
route that explains why HeLa's wild-type p53 does not respond.

## The goal is a product

*Set by Henry, October 2026: "we are not going to write any paper yet,
our goal is a product."*

So the measure of a piece of work is whether a cell biologist would
notice it, not whether it would make a figure. Academic deliverables
inherited from the original 2026 roadmap — the JOSS paper, the
validation preprint — are parked. `docs/VALIDATION.md` stays, because a
tool nobody can check is a tool nobody should trust, but it exists to
serve the user rather than a referee.

What that ordering implies, shortest path first:

1. **Someone hears about it and sees it work** — the hosted demo, live
   and self-updating. **Done.**
2. **They install it in one line** — `pip install cellsim`. **Done**
   (PyPI, October 2026).
3. **They run their own experiment without writing code** — **THE GAP.**
   The site can only replay runs made in advance. A biologist who does
   not write Python cannot ask it anything, and that is most biologists.
4. **They use it on THEIR cells** — plate import and calibration. Done.
5. **They trust the answer** — measured error bars and published misses.
   Partly done; see the honesty table below.

### Stage 3: an experiment anyone can run

The target user types no code. They open a link, choose from menus, press
Run, and get an answer with its error bar.

**It runs in the browser, with no server.** Small runs take 1.5 s
natively and Pyodide is a few times slower, so ten seconds for a
dose-response is realistic — fast enough to feel interactive. A server
would be faster and would also mean hosting bills, an attack surface and
something to maintain; a page that is just files on GitHub Pages will
still work in five years with nobody tending it.

**The design is borrowed; the code is not.** PhysiCell Studio is the
closest thing to this — a GUI over an agent-based cell simulator — and
its layout is worth learning from: tabs that separate *what the cells
are* from *what you do to them* from *what came out*. But it is **GPL
v3**, and copying its code would drag CellSim from MIT to GPL v3,
deciding our licence by way of a UI choice and shutting out commercial
users. So: study the interaction design, write our own code, keep MIT.
Where we do take code it must be MIT, BSD or Apache — Three.js (already
in the viewer) and Plotly are both fine.

Four screens, in the order a person actually thinks:

| Screen | What the user does | Non-scientist wording |
|---|---|---|
| **1. Cells** | pick a line, optionally knock out a gene | "Which cells?" |
| **2. Treatment** | pick a drug, dose, and schedule | "What do you do to them?" |
| **3. Run** | press one button, watch progress | "Run the experiment" |
| **4. Results** | 3-D view, curves, IC50 with its band, CSV | "What happened?" |

Plus **presets** — one click for "does pulsing beat continuous dosing?",
"what if these cells lose p53?", "how deep does the drug reach in a
tumour?" — because a new user does not know what to ask yet, and a blank
form is the commonest reason a tool goes unused.

Rules that keep it honest rather than merely pretty:

* every number comes with its band, and the UI never shows a bare point
  estimate;
* anything the engine is known to get wrong is labelled *in the result*,
  not in a document the user will not read;
* the preset that produced a figure is shareable as a URL, so a result
  can be checked by someone else.

### Stage 3b: inside a cell, then any drug (October 2026)

The Lab's second page shows one cell at a time: every species the engine
tracks, minute by minute, on a signalling map. Watching it is the most
direct check the engine has had — its first week found two wrong time
courses and, chasing them, seven drugs running placeholder constants
(`docs/VALIDATION.md`). The plan from there, in order:

1. **Fix what the view shows.** Done in 1.5.1: the library carries its
   fitted constants; p53 kills only with damage signalling, fitted to
   measured nutlin outcomes; notes on screen say what real cells do.
2. **Any drug by name or structure.** PubChem for the structure,
   ChEMBL's curated mechanism for the target, RDKit (its browser build)
   for drawing and for molecules nobody has made. A drug whose target the
   engine models runs with its potency *borrowed* from a library drug of
   the same class and labelled uncalibrated; any other target is refused
   by name. An impossible structure is rejected with the atom and the
   rule it breaks.
3. **More mechanisms on the existing map.** CDK4/6 inhibitors (cyclin D),
   BCL-2 inhibitors, ATM/ATR/CHK1, PARP, mitotic kinases — each a hook
   into a node already there, each calibrated on GDSC and checked on held
   out lines. In GDSC's 542 drugs, about 131 act on pathways the network
   already contains.
4. **Design mode.** Block or boost any molecule on the map, or give a
   drawn molecule a mechanism, and watch the cell. Labelled a
   hypothesis, kept apart from validated results.
5. **Growth signalling** (EGFR → RAS → ERK, PI3K → AKT → mTOR, feeding
   cyclin D), with per-line mutations from DepMap. About another 131 GDSC
   drugs act there. The largest single step.
6. **The cell's anatomy** — nucleus, mitochondria, spindle, membrane —
   with each species drawn where it acts.

What stays out: predicting a new molecule's target or strength from its
structure alone. Docking places poses well (87 % top-3) but its energies
are not reliable enough to set a potency, and nothing else available is.

### Stage 4: finish the honesty work

| Feature | State |
|---|---|
| Drug response, schedules, combinations | validated against GDSC and the published schedule literature |
| Spheroid O₂ and necrosis | calibrated on DLD-1, one documented miss |
| Gene knockouts | validated on four independent right answers |
| Migration | validated against random-walk theory |
| Immune killing | mechanism validated; absolute lysis unchecked |
| **Glucose** | **fails its validation and is labelled so** — fix or remove |
| Calibration to a user's plate | hold-out tested, 80 % band coverage |

### Stage 5: real users

Two or three biologists, their own questions, and a record of what
confused them. This cannot be compressed by writing code faster, and it
is the only stage that can tell us whether any of the rest worked.

## After Phase 3: what makes this useful to a working cell biologist

*Added October 2026. The list below is deliberately concrete; each item
names the data it would be checked against, because an unvalidated
feature is not an improvement.*

**The gap to close first.** Everything validated so far uses OUR cell
lines and OUR fitted constants. A biologist has their own line, their own
plate reader, and their own drug. Nothing in the repository yet turns
their numbers into a calibrated simulation of their cells. That is the
single highest-value piece of work left, and it is Phase 4.

### Phase 4 — the user's own cells

- *Done.* Read plate-reader exports (96/384-well with a plate map) and
  tidy CSV; normalise to controls, fit four-parameter curves with
  bootstrap confidence intervals whose coverage has been measured
  (`cellsim/plate.py`).
- *Done for measured data.* Report GR metrics (GR50, GRmax, GR_AOC;
  Hafner et al. 2016 Nat Methods 13:521), so fast- and slow-growing
  cultures are compared fairly. Still to do: the same metrics out of a
  simulation, so model and experiment sit on one axis. The engine
  already produces the counts they need.
- *Done.* Calibrate: fit the line's doubling time from its controls and
  each drug's single constant from its curve, with uncertainty, and say
  plainly when a curve cannot identify the constant
  (`cellsim/calibrate.py`).
- Carry that uncertainty into every prediction (sample the posterior),
  so schedules and combinations come out as bands, not lines.
- Validation: hold out a concentration range or a time point from the
  user's own plate and predict it; publish coverage of the bands.

### Phase 5 — biology the engine still cannot express

Each item is a mechanism with a validation target, ordered by how often
a cell biologist would hit it:

- **Quiescence at mitotic exit.** Real HeLa splits into fast cyclers and
  cells that do not divide for days (Spencer et al. 2013 Cell 155:369);
  the engine's cells all cycle. This is already a recorded miss against
  the Cell Tracking Challenge movie, and it changes every
  schedule-dependence answer.
- **Drug retention after wash-out.** The engine empties a cell as fast as
  it fills it, so no short exposure can kill. Tissue data say otherwise
  (Kuh et al. 2000 JPET 293:761 for paclitaxel). Validate against the
  published 1 h exposure IC50s the engine currently cannot reproduce.
- **Glucose, lactate and pH** in the dish, which set the necrotic
  boundary as much as oxygen does. Validate against Freyer & Sutherland's
  EMT6/Ro spheroids, where viable rim thickness was measured at four
  glucose and two oxygen levels.
- **Senescence** as a fate distinct from death, since DNA-damage agents
  produce it at the doses used clinically.
- **Resistance by mutation**, not only selection from a pre-existing
  tail: heritable changes at division, validated against the fold-
  resistance range that pulsed and continuous selection actually produce
  (2-8x; McDermott et al. 2014 Front Oncol 4:40).

### Housekeeping that keeps the record honest

Found the hard way, three times in one day, and worth fixing before it
bites a fourth: a change lands in one place while a DERIVED ARTEFACT
elsewhere goes on asserting the old world.

* `benchmarks/cell/*.json` are outputs, not assertions, so nothing fails
  when they go stale. Changing HeLa's parameters silently invalidated
  `dish_validation.json`, which kept reporting figures measured under the
  old ones.
* `scripts/validate_gdsc.py` carries its own copy of the IC50 estimator,
  so a change to `engine.ic50` left the largest published table untouched
  while appearing to succeed.
* A table of C50 figures went on quoting four significant figures after
  the text around it said they were approximate.

The fix is one target that regenerates every benchmark artefact and fails
if the committed file disagrees with what the code now produces —
`make validate`, or a CI job on a schedule, since the full set takes
hours. Until it exists, every change to a cell line or an estimator means
re-running `validate_gdsc.py`, `validate_dish.py`, `validate_plate_fit.py`
and `validate_calibration.py` by hand, and the duplicate estimator should
be unified away.

### Phase 6 — the parts that make it a tool others use

- A notebook-first API (`cellsim.api`) that returns data frames, because
  biologists live in notebooks and R, not in a CLI.
- Hosted demo and documentation site; conda-forge and PyPI releases with
  a pinned environment.
- Worked tutorials that answer real questions end to end: "should I pulse
  or hold", "when will my culture regrow", "does my combination
  antagonise".
- An experiment-design mode: given a question and a plate budget, propose
  the concentrations and time points that would best distinguish the
  competing answers.

### What stays out of scope

Ranking cell lines by sensitivity from markers (measured as unreachable,
see above); clinical dosing advice; and anything that claims to replace
the experiment rather than help design it.

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
