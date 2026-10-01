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

### Phase 1, the validated slice (target January 2027)

No UI work until this exists.

1. `cellsim/cell/engine.py`: one `CellState` dataclass, one integrator
   (SciPy LSODA or a fixed-step RK4 with sub-stepping), one `step(dt)`.
   Modules plug into it instead of being free functions.
2. Port from `OLD/`: Novak-Tyson cycle, ATM → p53 ⇄ MDM2 (3-variable
   delay oscillator) → p21 and PUMA → BAX → MOMP, with real units and
   the published rate constants already cited there. Re-tune the p53
   period to 5.5 h. Fix the XIAP/Smac balance so caspase-3 execution
   actually kills the cell.
3. Drug input: PK/disposition from `cellsim/cell` (permeation, pH trapping,
   efflux, binding sink) feeding occupancy at the target, then fate.
   Affinity and ADME from a small curated table with provenance
   (ChEMBL / BindingDB / literature), docking as a flagged fallback.
4. Validation harness against GDSC2 fitted dose-response.

   *Revised October 2026 after pulling the data* (`scripts/gdsc_reference.py`,
   extract in `benchmarks/cell/gdsc_reference.csv`). The original plan named
   cisplatin, doxorubicin and paclitaxel. What the public data actually
   supports:

   - **Paclitaxel is the quantitative target.** IC50 0.011–0.089 µM across
     our five lines, all inside the tested range, AUC 0.72–0.97. Gate:
     predicted IC50 within 3× on each line, and the correct ordering of
     the most and least sensitive line.
   - **Cisplatin cannot be a quantitative target.** In GDSC2 it is tested
     to 8 µM, and every one of our lines has a fitted IC50 *above* the top
     dose (20–122 µM) with AUC 0.92–0.99, i.e. the screen never reached
     50 % kill. Those IC50s are extrapolations. Gate it only on direction
     and rank: cisplatin must kill far less than paclitaxel at its tested
     range, and the model must not produce a sub-µM cisplatin IC50.
   - **Doxorubicin is not in GDSC2** for these lines; it is in GDSC1,
     which is a different assay generation. Either use GDSC1 explicitly
     and say so, or substitute **camptothecin** (in GDSC2, also a
     replication-coupled DNA-damage agent, so the same S-phase mechanism).
   - **The p53 split is not clean in this data** and must not be assumed:
     cisplatin IC50 is higher in p53-mutant HT-29 (122 µM) than in
     wild-type A549 (20 µM), as expected, but p53-mutant MDA-MB-231
     (43 µM) is *more* sensitive than wild-type MCF7 (100 µM). Report the
     p53 effect as a measured outcome, never as a gate that assumes it.

   Publish the table including every miss.
5. Headless API: `engine.run(schedule) -> per-tick state records`
   (JSON or Arrow). Any UI is a renderer of this stream.
6. A Colab notebook that reproduces the validation figure from a clean
   environment.
7. Packaging (done in Phase 0): `pyproject.toml`, the `cellsim` console
   script, pytest collection; the `src` package is now `cellsim`.

### Phase 2, the dish (target April 2027)

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
