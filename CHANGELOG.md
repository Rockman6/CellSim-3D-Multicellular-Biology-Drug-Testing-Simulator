# CellSim — Changelog

## Unreleased — Phase 3: ten drugs, calibrated intervals, evolved resistance

### Added
- Ten-drug panel covering every mechanism the engine can express, each
  with one fitted constant and GDSC screens on our lines: + etoposide,
  sn-38, gemcitabine, 5-fluorouracil, docetaxel, vinorelbine, nutlin-3a.
  Two new mechanisms for it: 's_phase' (damage only where replication
  forks are moving) and 'mdm2' (an inhibitor occupying MDM2's p53
  pocket), plus CellLine.p53_mdm2_independent_deg_per_h for the HPV E6
  route.
- `scripts/calibrate_cell_uncertainty.py`: prediction intervals scored
  leave-one-drug-out, so no width is tested on the errors that built it.
- `DishParams.mutation_rate` / `mutation_effect_sd`: heritable change at
  division, so resistance can arise during treatment rather than only be
  selected. `scripts/experiment_evolution.py` and its homogeneous-start
  control.
- `cellsim/api.py`: five tidy-table functions for notebook use
  (curve, exposure, washout, combination, spheroid) plus lines/drugs.
- `docs/cell_tutorial.md`: one real question answered end to end, every
  command run and its actual output pasted in.
- `web/index.html` landing page; `.github/workflows/pages.yml`, which
  simulates the demo runs at deploy time instead of committing them;
  `packaging/conda-forge/meta.yaml` (prepared, not submitted).
- `tests/test_pip_install.py`: builds a bare venv, pip-installs the
  project and runs the simulator there, asserting RDKit is absent — the
  gate behind "the cell simulator needs only NumPy and SciPy".

### Results
- The Phase-1 conclusion holds at ten drugs: engine closer than a
  constant-IC50 null on 34 held-out lines, further on 39, p = 0.64.
- THE EXCEPTION: nutlin-3a, whose target is part of the mechanism. All
  eleven lines called correctly (nine held out) against the null's five,
  including HeLa, resistant despite wild-type TP53 because HPV E6
  degrades p53 independently of MDM2.
- Intervals calibrated within 15 points of nominal at 68/90/95 %.
- Resistance evolves from an identical population (1.85x) where
  selection alone can produce nothing, and a hard pulsed schedule
  eradicates before it gets started.

### Fixed
- `tests/cell/test_gdsc_targets_smoke.py` asserted the reference covered
  exactly three named drugs, so it would have passed while validating
  nothing as the panel grew. It now checks against the library, and
  immediately found DLD-1 and SW480 match no GDSC screens — both
  genuinely absent from release 8.4, now recorded rather than silently
  substituted with a sister line.
- `scripts/gdsc_reference.py` matched cell-line names exactly, silently
  dropping HCT116 (GDSC: HCT-116). New `CellLine.gdsc_name`.
- The IC50 estimator used a single 4.6-fold-step grid, which moved a
  fitted line's own prediction by up to 15 % between runs; it is now
  two-pass.

## Unreleased — Phase 2: the dish

### Added
- `cellsim/cell/dish.py` and `cellsim dish`: the engine's cells on a
  lattice, one per site, each with a full engine state. Monolayer (contact
  inhibition) or spheroid (radial oxygen and drug fields, hypoxic
  quiescence, necrosis, apoptotic-body clearance), with daughters
  inheriting their founder's traits. Field solvers checked against closed
  forms. Stream schema `cellsim.dish.stream/v1`.
- `web/viewer/`: a static web viewer for both stream schemas (three.js):
  3-D cells coloured by phase, p53, caspase-3, damage, oxygen, growth
  signal or lineage; cut-open spheroids; time scrubbing; charts; cell
  inspector; CSV export. Runs no biology. Demo runs in `web/demo/`.
- `scripts/ctc_reference.py`: counts and cell-cycle times from the Cell
  Tracking Challenge HeLa ground truth (Kaplan-Meier for censoring);
  `scripts/validate_dish.py`: the dish against that movie and against
  DLD-1 spheroids (Grimes et al. 2014).
- Engine: a per-cell G1 growth signal (`aux["gs"]`) for contact
  inhibition and hypoxia, and optional cycle-time variability
  (`Params.cycle_cv`); with neither used, results are bit-identical.
- DLD-1 in the library (TP53 S241F, high ABCB1), and
  `Drug.accumulation_ratio` for tissue penetration.

### Fixed
- HeLa's doubling time was 20 h, read off the CTC movie's complete
  cycles, which a 46 h window biases short. Now 31 h (Cellosaurus 1.3 d;
  the movie's own counts double every 27-31 h). HeLa's held-out GDSC
  predictions moved 3-7 % and still pass.

### Known misses (documented, open)
- Spheroid necrotic cores are larger than Grimes's anoxic-core relation
  allows (implied diffusion limit below the measured 233 um).
- The HeLa movie's cells split into fast cyclers and non-dividers
  (proliferation/quiescence at mitotic exit); the engine's do not.

## Unreleased — Phase 2 audit: refit, then re-measure everything

### Fixed
- The library's fitted potency constants predated the apoptotic buffer
  and the efflux channel, though the docstring said a re-run reproduced
  them within 1 %. Refitted under the current defaults: cisplatin
  k_damage 0.000951 -> 0.001181, doxorubicin 0.00794 -> 0.01749,
  paclitaxel partition 0.153 -> 0.212 (`scripts/validate_gdsc.py`).

### Changed
- Experiments state doses as multiples of each drug's own engine IC50,
  computed at run time, so refitting the library cannot change what a
  design means. Tests pin the potency their doses were designed at
  (exactly a dose rescaling for cisplatin and paclitaxel).
- `scripts/experiment_schedule.py` now measures the iso-effect curve
  C50(T) for exposures of 1-72 h and compares it with Ozawa 1989,
  Liebmann 1993 and Georgiadis 1997, including one documented miss
  (24 h paclitaxel exposure is far too effective in the engine).
- `scripts/experiment_combination.py` reports five seeds and two dose
  designs, mean +- sd.

### Corrected (see docs/VALIDATION.md, "Phase 2: audited")
- RETRACTED: "pulsing selects a pharmacokinetic escape, not an apoptotic
  one". It rested on one surviving lineage. Measured exactly on a trait
  grid, the direction is the reverse for all three drugs: continuous
  exposure selects on drug accumulation, pulsing shifts selection
  toward the apoptotic reserve.
- QUALIFIED: the platinum/taxane ordering effect. Absent with both drugs
  at their IC50; ~12 % at 72 h with cisplatin at 1.25x IC50, fading to
  ~3 % by 120-144 h. The earlier "0.93 -> 0.98 -> 0.99" was one seed.
- REPLACED: the "242x / 17x / 7.6x" exposure-time ranking, which
  depended on the dose chosen. The iso-effect curve keeps the ranking.

### Added
- `cellsim.cell.engine.ic50` (two-pass IC50), `exposure_h` on
  `dose_response`/`ic50` (pulse then wash-out), and `traits=` on
  `simulate` to set each cell's heterogeneity instead of drawing it.
- `cellsim/cell/selection.py`: survival landscapes over trait space and
  exact expected selection under the population's trait distribution.

## Unreleased — Phase 1, first validation cycle

### Added
- `cellsim/cell/engine.py`, the single-cell drug-response engine: the
  CDK/cyclin cycle, ATM → p53 ⇄ MDM2 → p21 / PUMA → BAX → MOMP → caspase,
  ported from `OLD/` with cited constants. Cells now actually die (the
  prototype's caspase-3 was clamped at 0.04) and the p53 pulse period is
  4.9 h (published 5.5 h; prototype 3.0 h).
- Cell-to-cell variability: log-normal (σ 0.3) multipliers on each cell's
  Bcl-2 reserve and drug uptake (Spencer 2009). Widens the 85→15 %
  viability window from 1.2× to about 2×; still far steeper than a Hill
  slope of 1, documented as open.
- Dosing schedules and a UI contract: `cellsim/cell/stream.py` emits
  per-tick JSON Lines (schema `cellsim.cell.stream/v1`) for a population
  under a piecewise-constant schedule, wash-outs included.
  `cellsim cell-stream`; also `cellsim cell-sim` for a dose-response.
- `scripts/gdsc_reference.py` (30 GDSC1/GDSC2 screens, replicate spread)
  and `scripts/validate_gdsc.py` (`cellsim validate-gdsc`): fit one
  constant per drug on A549 + MCF7, predict every line, compare with a
  constant-IC50 null.

### Added (continued)
- `Params.p73_gain`: the ATM/c-Abl -> p73 -> PUMA route by which
  TP53-mutant cells still die of DNA damage (Gong 1999, Agami 1999).
  Default 0 (p53-only). `scripts/experiment_p73.py` tests it
  leave-one-mutant-out: the gain is chosen on one mutant line and
  scored on the other, and it improves out-of-sample prediction in
  both folds while beating the constant-IC50 null in both.

### Fixed — streptavidin reference data (closes #14)
- `benchmarks/dock/streptavidin_calibration.yaml` listed desthiobiotin at
  K_d 5e-5 M, ~6 orders of magnitude too weak and inconsistent with the
  source cited on the entry itself. Corrected to ~1e-11 M / -15.0
  kcal/mol (Hirsch 2002 Anal Biochem 308:343; bracketed by Green 1990
  Methods Enzymol 184:51, which puts desthiobiotin 10^2-10^4 fold weaker
  than biotin).
- Consequences followed through, not hidden: measured Spearman on this
  set falls +0.80 -> +0.40 and MAE rises 4.98 -> 6.66 kcal/mol. The
  claim that docking ranks biotin-site binders usably does not hold, so
  the class moved out of the tutorial's "reliable for ranking" table
  into "pose filter only", next to kinases.
- `tests/uq/test_calibration_smoke.py` asserted Spearman >= 0.8, which
  passed only because of the bad number. It now pins the saturation
  itself: four compounds spanning 12.0 kcal/mol of experimental affinity
  produce predictions spanning 0.33 kcal/mol. That is reproducible and
  is what the reliability table's `do_not_trust_absolute` verdict rests
  on.
- Regenerated from corrected data: `benchmarks/dock/uq_coverage.json`,
  the `ultra_tight_binder` row of `benchmarks/dock/reliability_table.yaml`.

### Validation result — exit gate NOT met, and now measured properly
- Expanded to **ten** cell lines (five TP53 wild-type, five mutant),
  fitted on two, so eight are held out.
- **Established:** each drug's potency scale from one fitted constant.
  22/25 held-out predictions inside GDSC's replicate span vs the
  constant-IC50 null's 18/25; A549 cisplatin 9.0 µM vs 9.8 measured.
- **Not established:** line-to-line discrimination. Paired per-line
  sign test on held-out lines — engine closer on 10, null closer on 11,
  p = 1.0.
- **Corrects an earlier entry in this file.** A five-line run appeared
  to meet the gate; the gate then compared two RMSE numbers and the
  margin that passed it (0.45 vs 0.46) was inside the noise. The gate
  is now a paired sign test that cannot pass on a difference that small.
- The p53-independent route still earns its place: it is a large
  improvement over the p53-only model on four of five mutants (which
  that model misses by up to 12×), but it beats the null on only two and
  hurts T47D, the one genuinely resistant mutant. Recorded in
  `docs/VALIDATION.md`.

### Diagnosed — why the engine cannot tell cell lines apart
- Swept every per-line field across its full observed range and measured
  how far the predicted IC50 moves: doubling time 1.0x, Bcl-2 reserve
  1.1x, DNA repair rate 7-9x (and 1.0x for paclitaxel). Real lines
  differ by 8-51x. Doubling time is the engine's main line-specific
  input and has no leverage at all, so with TP53 status binary the
  engine can give at most two answers per drug — tying a constant is
  arithmetic, not misfortune.
- Paclitaxel has no line-specific channel whatsoever: its IC50 is set by
  the tubulin-occupancy arrest threshold, a property of the drug. It
  needs a per-line efflux channel the engine does not model.
- Pulled DepMap 24Q4 expression for all 10 lines and 16 candidate genes
  (efflux, BCL2 family, repair, target). After correcting for having
  searched 16 genes, by permutation, nothing correlates with sensitivity
  better than chance (family-wise p 0.86 / 0.25 / 0.54). The top hit for
  both cisplatin and doxorubicin is TUBB3, which has no mechanism for
  either. DepMap's CRISPR growth rate also fails as a scalable
  doubling-time substitute (rho -0.15 against curated values).
- Conclusion recorded rather than coded around: building an expression
  mapping now would fit noise, the same error behind the retracted gate
  above. The blocking constraint is the number of cell lines.

### Measured — what per-line data actually predicts, and why the engine misses it
- Repeated the expression test on every GDSC line DepMap covers
  (156-683 per drug instead of 6-10), genes pre-specified by mechanism,
  Bonferroni-corrected per drug. Real signal appears:
  doxorubicin — BCL2L1 rho +0.19 (p 5e-06), ABCB1 rho +0.17 (p 9e-05);
  paclitaxel  — BCL2L1 rho +0.31 (p 9e-10), ABCB1 rho +0.23 (p 1e-05),
                TUBB3 rho +0.22 (p 3e-05);
  cisplatin   — nothing, and the repair genes in particular do not
                predict it (ERCC1 rho -0.03), consistent with ERCC1's
                record of failing to replicate as a platinum biomarker.
- Side by side with the leverage sweep this explains the tie exactly:
  the engine's one strong channel (repair, 7-9x) is where the data has
  no signal; the strongest real signal (apoptotic set-point, implying
  ~2x) maps to a field the engine already has but throttles to 1.1x
  because its MOMP switch is too sharp; and efflux, the second
  strongest, is not modelled at all.
- Ceiling recorded before any build: the real effects are ~1.5-2.3x
  each, composing to perhaps 3-5x against the 8-51x lines actually
  span. Wiring them in should beat a constant but will not close the
  gap.
- `scripts/experiment_expression_scale.py` added; its p-value routine
  was checked against scipy to 4 decimal places before use.

### Added — Phase 2: real two-drug combinations
- The engine carries one intracellular concentration per drug
  (`state_width`, `N_BASE`), so a combination is a genuine two-drug run.
  `simulate(line, [drugA, drugB], [c_a, c_b], ...)`; a single drug or
  None behaves exactly as before and the whole suite is unchanged.
- Antagonism reproduced and gated: cisplatin arrests the cycle via
  p53 -> p21 and so removes the cells paclitaxel needs in mitosis. Every
  arrangement kills 1.3-1.5x LESS than independent action predicts.
  Nothing fitted; both halves were already in the engine.
- The apparent sequence effect is an ARTEFACT, which is the more useful
  result. Platinum-first looks better at a 72 h readout (ratio 0.93) and
  the advantage dissolves as the readout lengthens (0.98 at 96 h, 0.99
  at 120 h), because cisplatin acts slowly and placing it second
  truncates its effect inside a fixed window. So the engine does NOT
  reproduce a biological platinum/taxane sequence dependence; it exposes
  the timing confound that fixed-endpoint assays fall into.
- `tests/cell/test_combination_timing_smoke.py` (3 gates, 22 s) pins the
  antagonism and the artefact, and deliberately does not assert a
  sequence effect.

### Added — Phase 2: selection under treatment
- `scripts/experiment_resistance.py` +
  `tests/cell/test_resistance_selection_smoke.py` (5 gates, 32 s).
- Selection emerges without being modelled: survivors of a week of
  cisplatin are enriched for high anti-apoptotic reserve and low drug
  accumulation, while untreated controls stay at the population mean.
- Dose intensity decides how hard AND along which axis. At matched total
  exposure, 3x the concentration for a third of the time kills 13x more
  and leaves a more resistant remnant (combined score 2.88 vs 1.87) --
  but selects almost entirely on drug accumulation (0.33 vs 0.61), not
  on the apoptotic reserve, because no reserve survives a triple dose.
- Consistent with the schedule result built independently: pulsing helps
  cisplatin and doxorubicin most and paclitaxel least.
- Limits recorded in the test and docstring: heterogeneity is drawn once,
  so this is selection from a pre-existing tail rather than evolution of
  new resistance, and lineages do not compete for space, so
  adaptive-therapy questions are out of reach here by construction
  (`cellsim/cell/agents.py` has the lattice for those).

### Added — Phase 2 opens: schedule dependence
- `scripts/experiment_schedule.py` + `tests/cell/test_schedule_smoke.py`
  (4 gates, 38 s). Entirely within-line comparisons, so the per-line
  prediction ceiling does not apply.
- Holding paclitaxel's AUC fixed and varying the split gives an INTERIOR
  optimum at 36 h (surviving fraction 1.00 / 0.48 / 0.28 / 0.69 / 1.00
  across 3 / 24 / 36 / 48 / 72 h). Both extremes fail for different
  reasons: too brief and few cells reach mitosis; too dilute and tubulin
  occupancy never crosses the arrest threshold. Reproduces across seeds
  and is a prediction of the mechanism, not a fit.
- Exposure-time dependence ranks paclitaxel 242x, doxorubicin 17x,
  cisplatin 7.6x (survival at 3 h over 72 h at fixed concentration) --
  the published ordering, and the reason AUC is the wrong exposure
  metric for taxanes (Gianni 1995, Huizing 1993).

### Added — per-line efflux channel, and the ceiling it exposed
- Saturable P-glycoprotein efflux on intracellular drug, for substrate
  drugs only (doxorubicin, paclitaxel; not cisplatin), scaled by
  `CellLine.efflux_level`. That level is DERIVED, not fitted:
  2^(ABCB1 log2(TPM+1) - panel median) from DepMap 24Q4, a fixed formula
  with no free parameter, so held-out lines stay out of sample.
  Leverage 1.9x / 2.0x, matching the 1.5-2.1x the data implies, and
  exactly 1.0x for cisplatin. `Params.efflux_vmax_uM_per_h = 0` disables.
- It did NOT improve prediction: paired count 11:10 -> 10:11, a one-line
  change at p = 1.0. HeLa carries ~10x the median ABCB1 yet sits
  mid-pack for doxorubicin, so the engine now calls it too resistant.

### Measured — the ceiling on per-line prediction from canonical markers
- Fitting ALL eleven markers to IC50 by least squares (far more freedom
  than a mechanistic model gets) across every GDSC line DepMap covers:
  best-case R2 = 0.09 (cisplatin, n=156), 0.15 (doxorubicin, n=683),
  0.20 (paclitaxel, n=427). So 80-92 % of line-to-line variation is not
  captured by these markers at all.
- Consequence, recorded rather than worked around: a mechanistic model
  restricted to these inputs cannot pass the paired sign test however
  faithfully it is built, and building it more faithfully (the efflux
  channel) did not help. rho ~ 0.2 is a population trend, not a
  per-line predictor.
- Phase-2 direction revised: aim at WITHIN-LINE dynamics (schedule
  dependence, wash-out recovery, resistance under selection, combination
  timing), where the engine's mechanism is load-bearing and no per-line
  marker is needed. Per-line IC50 is dropped as the headline metric.

### Changed — Bcl-2 reserve is now a stoichiometric buffer, not a rate penalty
- The anti-apoptotic reserve divided the Bax ACTIVATION RATE, so a
  high-Bcl-xL line merely died later: a slower climb still crosses a
  fixed MOMP threshold within 72 h. That is why an 8-fold change in
  `CellLine.bcl2_level` moved the IC50 by only 1.1x. It now acts as a
  buffer -- only Bax in excess of the reserve can form pores -- which is
  the mitochondrial-priming picture (Certo 2006 Cancer Cell 9:351).
  `Params.bcl2_buffer = 0` restores the old behaviour.
- Effect, as predicted by the ceiling recorded above: doxorubicin
  leverage 1.2x -> 1.7x (data implies 1.8x) and its held-out RMSE
  0.45 -> 0.40 against the null's 0.46; the paired sign test flipped
  from 10:11 to 11:10 but stays at p = 1.0, so the gate is still NOT
  met. A 1.7x channel cannot resolve an 8-51x spread.
- Cisplatin got worse (0.54 -> 0.65), consistent with it being the drug
  with no Bcl-xL signal.
- Side benefit: TP53-mutant lines now keep PARTIAL resistance (survival
  0.08 vs wild-type 0.00 under sustained damage) instead of being killed
  outright, which was previously flagged as stronger than the
  literature supports. Pinned by a rewritten test.

### Changed
- `cellsim/cell/library.py`: placeholder potencies replaced by the fitted
  values (refitting reproduces them within 0.3 %). Paclitaxel's fitted
  constant is now its effective intracellular accumulation, because the
  mitotic-arrest threshold, not the death rate, sets its 72 h IC50.

### Fixed
- `cellsim cell-sim` crashed on a wrong result-field name.

## v1.4 — 2026-10-01 — Phase 0: one identity, green CI, July fixes merged

The project is refocused on one product: a validated, interactive
cell and drug-response simulator. See `docs/PLAN.md`.

### Merged
- PR #11 — Meeko atom-order scramble fix. Docked poses are now
  reconstructed via Meeko's reverse conversion; every earlier
  pose-recovery number was computed on scrambled molecules
  (biotin/1STP reported 2.02 Å, real value 0.56–0.69 Å). New
  15-cocrystal blind set, RCSB-fetched: top-1 73 %, top-3 87 %,
  PoseBusters 100 %.
- PR #12 — CYP3A4 SoM ranks C–H bonds only; `src/cell` (20 PK/PD
  modules), `src/uq/reliability.py` target-class accuracy table,
  Sobol and coverage artifacts.
- PR #10 — all 13 BUG_AUDIT fixes in the FEP binding path
  (standard-state sign, leg force-field mismatch, restraint
  geometry, seed plumbing, stereo).

### CI
- The 29 test files added by those PRs are now gated (54 CI steps).
- Root cause of the red badge since 2026-05-30: the 1STP re-dock
  gate sat on the 2.0 Å noise boundary because of the scramble bug.
  Fixed by PR #11; the gate now reads 0.56 Å.

### Repository
- README rewritten around the single product identity; GitHub
  description and topics corrected (no longer "C++20/Metal").
- `LICENSE` is pure MIT again; attributions moved to
  `THIRD_PARTY_NOTICES.md` so GitHub detects the license.
- Campaign-era planning docs (GOAL, MISSION, ROADMAP, BENCHMARKS,
  campaign scope/status/close-out, professor debrief, friend
  handoff) moved to `docs/archive/campaign1/` with the FreeSolv
  pilot-3 report. `TUTORIAL.md` → `docs/docking_tutorial.md`.
  OLD-specific references (`UNITS.md`, LaTeX reference) moved under
  `OLD/`.
- Removed: prof-email / csv_tldr / finalize_run handoff tooling and
  their tests, committed `run/` logs, the retired Metal-UI
  screenshots, and the empty `cg`, `core`, `render`, `viewer`
  scaffolds (`core` was a byte-identical copy of `OLD/src/core`).
- 18 prototype-era scripts (C++ calibration, validation, window
  capture, molecule generation) moved from `scripts/` to
  `OLD/scripts/` with their paths fixed; they had pointed at the
  pre-restart layout since April. `OLD/README.md` added with the
  verified headless build.

### Packaging
- The `src` package is renamed `cellsim` (every import, `python -m`
  invocation and path reference rewritten; `git mv` keeps history).
- `pyproject.toml` added: `pip install -e .` inside the conda env
  registers the package and a `cellsim` console command
  (`cellsim/cli.py`, a Python twin of `scripts/cellsim`). `pytest`
  collects `tests/` directly.

### Known, documented limits carried forward
- Hydration FEP runs NVT on a packmol box with no barostat; the
  FreeSolv-12 "PASS" (MAE 1.42 on 10/12) hides a size-dependent
  bias (toluene −3.7 kcal/mol). Binding FEP has never produced a
  number on a real binder. Both are experimental, not product paths.
- SoM is advisory (2/3; blind to N-dealkylation).
- `OLD/` biology: cells never execute death under cisplatin
  (caspase-3 clamped), p53 period 3.0 h vs 5.5 h literature.

## v1.3 — 2026-05-23 — Milestone A passes (FreeSolv hydration FEP gate)

**The chemistry-axis half of the professor's "closed ontology /
circular validation" critique is closed.** Hydration ΔG_hyd on
the FreeSolv-12 benchmark now reproduces published experiment to
MAE 1.42 kcal/mol (gate ≤ 1.5), with the full alchemical pipeline
running end-to-end on Apple-silicon CPU and producing a tarball
that goes through `cellsim fep-report` → `csv_tldr.py` →
`fill_prof_email.py` for prof-ready handoff. Tagged as
`milestone-a-pilot-3` on `origin`.

### Milestone A pass numbers (10/10 completed compounds)

- MAE = 1.42 kcal/mol (gate ≤ 1.5): **PASS**
- Pearson r = +0.913
- Spearman ρ = +0.903
- Kendall τ = +0.778
- GHMC acceptance: 99% mean, 99% worst (gate ≥ 70%): **PASS**
- Sign-critical (methane positive, hydrophobe direction): **PASS**
- Wall: 7.5 h on Henry's MacBook Pro M-series CPU
- Compounds completed: 10/12 (acetic_acid + acetamide need
  bumped sampling per translator hint — `--n-windows 22
  --n-production-steps 75000` queued)

### Diagnosis trail (BENCHMARKS.md § "Milestone A post-mortem")

Six hypotheses tested over four investigation passes; five
rejected (padding / sampler / softcore / PME / FF-bug); one
landed (`b89dd51` split-schedule decoupling, electrostatics
first then sterics, prevents water penetration into the softcore
region at intermediate λ — yank/perses standard). One wrong fix
shipped and reverted (`c461053` sign-flip looked right under
smoke sampling, was wrong under production sampling, reverted in
`8ab37a1` once today's pilot-3 production data proved the
original Hummer-Szabo direction). Lesson: never use under-
sampled FEP output to assert physics — the methane sign test
now inspects source code directly instead of relying on a
sampling-noise-biased value.

### Added (Milestone A infrastructure)

- `_split_lambda_schedule` in `src/fep/sampling.py` — canonical
  absolute-hydration FEP schedule (decouple elec first, sterics
  second). Replaces the legacy coupled `λ_elec = λ_sterics` that
  let waters penetrate the softcore well and inflated F(decoupled)
  by 4-6 kcal/mol.
- `_biologist_reason_for_mbar_error` — translator that converts
  pymbar internal exceptions ("column sum W_nk = 0") into
  biologist-actionable guidance ("adjacent λ-windows don't overlap;
  try --n-windows 22 --n-production-steps 75000"). 10/10 regression.
- `--resume` + atomic incremental CSV write on `cellsim
  fep-freesolv` (and parallel pattern on binding bench from
  earlier). Multi-day CPU runs survive lid-close / power loss
  without losing completed compounds.
- `scripts/csv_tldr.py` — one-line Slack/standup summariser of a
  bench CSV, exit-code-gated, verdict logic aligned with
  `cellsim fep-report`. 10/10 regression.
- `estimate_sampling_wall_hours` + `format_wall_estimate_block`
  in binding — auto-prints projected CPU / Metal / CUDA wall time
  after any scaffold-only run, with a "> 48 h CPU not viable" flag.
- `--hardware` and `--platform` flags now flow into the hydration
  template in `scripts/fill_prof_email.py` (was only binding before).
- Provenance fallback in `cellsim fep-report` — scans `env.log`
  and `doctor.log` for `git commit:` when `run.log` doesn't have
  it, so legacy / pre-fix tarballs are still traceable.
- `docs/milestone_b_run.md` — friend handoff for the M5 Max
  binding runs (streptavidin + EGFR).

### Fixed

- Hydration sign convention restored to Hummer-Szabo `dG_hyd =
  -ddg = vac - solv` after a smoke-test-driven mis-fix shipped
  the wrong sign (8ab37a1).
- GHMC integrator timestep + friction tuned per prof's checklist
  (acceptance > 99% under split schedule).
- env.log / run.log header-tee fix so `git commit:` reaches the
  tarball provenance.

### Tested (smoke gates added)

| Suite | Cases |
|---|---|
| fep-freesolv --resume + incremental CSV | 7/7 |
| MBAR-error translator | 10/10 |
| wall-time estimator | 10/10 |
| csv_tldr verdict alignment | 10/10 |
| fep-report provenance fallback | 7/7 |
| fill_prof_email hardware override | 10/10 |
| hydration_dg (post-revert composition pin) | 3/3 |

All 15 FEP suites + 107+ FEP test cases green on commit `8643604`.

## v1.2 — 2026-04-19

Chemistry-first bioagent foundation + rigorous apoptosis + faster sim.

### Added

**Multi-threshold apoptosis engine**
- Every cell owns an `ApoptosisEngine` (Bcl-2 / Bax / Bak / MOMP / caspase-3/8/9 / Smac / IAP cascade) driven by 11 input triggers (p53 from DNA damage, ROS, ATP collapse, hypoxia, survival-factor loss, anoikis, crowding, replicative, Fas ligand, mito dysfunction, drug pressure).
- `ApoPhase` state machine: ALIVE → PRIMED → MOMP → EXECUTION → FRAGMENTATION → BODIES → CLEARED with literature-calibrated dwell times (Ye 2017, Spencer 2009, Silva 2010).
- Apoptotic bodies spawned per dying cell (5-15 bodies, 0.1-0.5 wu radius) with cytosol / membrane / receptor ledgers; drift with fluid current; decompose over 3 bio-h; secondary necrosis at 6 bio-h via nanopore-driven osmotic swelling (Chen 2016).
- Mass-conserved release to `MediumField` through the closed-system `exchange()` API — tracked per species (AA, water, ions, Ca²⁺, pyruvate, lactate, glucose).
- Efferocytosis: live cells within 4 µm absorb bodies into lysosomal pools; 85 % recycle efficiency into biomass (Ravichandran 2015, Appelqvist 2013); extracellular proteases released from death sites boost neighbor substrate availability.
- DAMP signalling: lysed bodies emit Ca²⁺ + stress bumps to neighbors within 2 wu (Bratton 2010).

**Chemistry-first emergent drug system (MOA-free)**
- Replaces the old MOA-labelled `DRUG_LIBRARY`. The simulator no longer hardcodes what a drug *does* — it computes binding affinities and effects from the drug's chemical structure.
- [src/simulation/biochem/ChemicalEntity.h](src/simulation/biochem/ChemicalEntity.h) — universal record (kind, structure, partial charges, LJ, Morgan FP, pharmacophores, binding affinities, reactive tags) for drugs / metabolites / organelles / viruses / bacteria.
- [src/simulation/biochem/Bioagent.h](src/simulation/biochem/Bioagent.h) / [.cpp](src/simulation/biochem/Bioagent.cpp) — registry with `loadFromDisk()` that reads [data/bioagents/drugs.csv](data/bioagents/drugs.csv) + per-drug Tier-0 JSON caches.
- [src/simulation/biochem/BindingMatcher.h](src/simulation/biochem/BindingMatcher.h) — descriptor-based drug × target Kd (Hill-like physicochemical fit, calibrated 0 → 100 mM, 1 → 1 nM).
- [src/simulation/biochem/TargetLibrary.h](src/simulation/biochem/TargetLibrary.h) / [.cpp](src/simulation/biochem/TargetLibrary.cpp) — 13 built-in druggable targets (CDK1, CDK2, Rb, TopII, tubulin, Bcl-2, DNA intercalation, HDAC, EGFR, proteasome, ribosome, Pol II, Complex I) with per-target function modulators that perturb *existing* `SimCell` fields — no new outcome code paths.
- `SimCell.applyDrug()` + `updateDrugPK()`: Lipinski-derived permeability, per-cell target binding (k_on × free × free_target − k_off × bound), occupancy-driven modulators.
- [scripts/chem_worker.py](scripts/chem_worker.py) — background RDKit worker computing Tier-0 descriptors + Morgan fingerprints + 3D SDF for any SMILES. Five seed compounds pre-computed.

**Drug Lab UI**
- Right-side ImGui panel: compound dropdown (from registry), live descriptor readout, log-scale concentration slider, `APPLY UNIFORMLY` button, top-affinities table.
- Drug molecules render as bright yellow ball-and-stick particles via the existing `MediumChemical` swarm; drift with curl-noise current, home on cells, bind at membrane, stream inside — reusing the `MC_FREE → ATTRACTED → BINDING → TRANSPORT` state machine that already handles metabolite uptake.

**Export Everything button**
- Native folder-picker dialog → timestamped export bundle with 6 files: `manifest.txt` (run summary + mass-balance drift), `population.csv` (time-series), `cells.csv` (per-cell snapshot), `medium_field.csv` (grid × 12 species), `apo_bodies.csv`, `medium_chems.csv`. `Reveal in Finder` button.

**Validation & tooling**
- [scripts/compare_against.py](scripts/compare_against.py) — multi-schema comparator (CTC, viability, growth, single-cell).
- [scripts/validate_unseen.sh](scripts/validate_unseen.sh) — runs headless against every reference dataset.
- [scripts/scale_test.sh](scripts/scale_test.sh) — (init cells × bio-h) sweep.
- [scripts/calibrate_rate.sh](scripts/calibrate_rate.sh) — rate-parameter grid search.
- [scripts/derive_ctc_counts.sh](scripts/derive_ctc_counts.sh) — converts raw CTC `man_track.txt` into cell-count CSVs.
- New DIC-C2DH-HeLa reference curves (2 sequences, 14 bio-h each).

### Changed

**Calibration**
- `SLOW_DT_SCALE`: 0.052 → 0.055, `MEDIUM_DT_SCALE`: 0.190 → 0.201, `MECH_P21_COUPLING`: 0.018 → 0.002.
- Cell-count mean-relative error against CTC Fluo-N2DL-HeLa seq02: **6.4 % → 5.1 %**.
- Calibration metric reproducible across 5 independent runs.

**Apoptosis engine integration cadence**
- Rewrote `ApoptosisEngine::step()` as `stepOnce()` + outer loop with `DT_MAX = 0.08` bio-s sub-steps. Previous hard-clamp `dt = fminf(dt, 0.01)` left the cascade effectively frozen over a physics tick.

**Drug subsystem rewrite**
- Deleted the hardcoded `Drug` struct / `DRUG_LIBRARY` / `MOA_*` enum / `updateDrugResponse()` / `MS_DRUG` species — everything the old PK/PD model bolted onto the simulator.
- `MS_COUNT` reduced from 13 → 12; `MS_WATER` shifted to index 11.

**Rendering**
- Dying cells: shrinkage (up to 35 %), warm red-brown tint, PS-flip "eat-me" warmth, alpha ghost fade, bleb-driven `furrowDepth` oscillation — all via existing cell-shell shader, no shader changes.
- Apoptotic bodies: warm-tinted translucent spheres; nuclear fragments deep-red, mito fragments green-tinged (residual cyt-c).
- [MoleculeCache](src/molgen/MoleculeCache.h) now also scans `assets/drugs/*.sdf` at startup so user-added drugs render automatically without edits to the hardcoded library list.
- `CellInstance` extended with `apoProgress`, `blebPhase`, `nuclearFragmentation`.

**Telemetry**
- Removed drug-related columns (`viability_pct`, `avg_drug_damage`, `med_drug_uM`); cleaner schema.

### Fixed

- Avogadro scaling in drug target-binding math (was 6.02e17, corrected to 2.41e6 = N_A × V × 10⁻⁶).
- Closed-system mass invariant now holds through apoptosis + body decomposition + efferocytosis cycle (`MediumField::checkBalance()` stays < 1e-3 drift).
- `TRACKING_LOSS_PROB_PER_BIOSEC` scaffold kept at rate=0 by default, gated on rate > 0 so dead code doesn't perturb RNG state.
- "Halved cells" / floor-clipping fixed: `c.position.y = FLOOR_Y + c.radius * c.size * 0.85f` at init.

### Tested

| Scenario | Metric | Value |
|---|---|---|
| CTC Fluo-N2DL seq02 (calibration) | mean \|rel err\| | 5.1 % |
| CTC seq02 reproducibility | 5 runs | all 5.1 % |
| Scale: 1200 cells × 48 bio-h | wall time | 65 s |
| Scale: 800 cells × 120 bio-h | deaths emerged | 687 (chronic-crowding apoptosis) |
| Drug matrix (125→186 warmup + 48 h drug, every compound) | proliferation inhibition | 75 % (257 → 64 divisions) |
| Mass balance after 10+ apoptotic events | max drift | < 1e-3 |

## Earlier releases

See `git log` for commits prior to v1.2.
