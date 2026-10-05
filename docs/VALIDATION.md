# Validation: the numbers that are actually true today

Every row names its reproducer. Every caveat is part of the number.
The long form, including the post-mortems that produced these numbers,
is in [`archive/campaign1/BENCHMARKS.md`](archive/campaign1/BENCHMARKS.md).
Last full re-check: 2026-10-01 on an Apple-silicon laptop in the pinned
`cellsim` environment.

## Docking: blind pose recovery (Layer 1.3 / 1.7)

15 diverse cocrystals; structures and ligand SMILES fetched from RCSB at
run time, so the benchmark is not circular on bundled inputs.

| Metric | Result | Gate |
|---|:-:|:-:|
| PoseBusters physical validity | 100 % (15/15) | ≥ 95 % |
| Pose recovery, best of top-3 ≤ 2 Å | 87 % (13/15) | ≥ 75 % |
| Pose recovery, top-1 ≤ 2 Å | 73 % (11/15) | ≥ 75 %, just under |

Reading: Vina's sampling finds the right pose; its scoring ranks it
first less reliably. A UFF-strain re-rank does not recover the misses.
Earlier numbers (the "2/3 mini-bench") were computed on molecules
scrambled by a Meeko atom-order bug; biotin/1STP reported 2.02 Å then,
and reads 0.56–0.69 Å now.

Reproduce: `python scripts/run_blind_dock_bench.py benchmarks/pdbbind/blind_set.yaml`
(results committed in `benchmarks/pdbbind/blind_set_results.csv`).

## Corrected: the streptavidin reference set (October 2026)

`benchmarks/dock/streptavidin_calibration.yaml` listed desthiobiotin at
**K_d = 5×10⁻⁵ M** (ΔG −5.9 kcal/mol). That was wrong by about six
orders of magnitude and inconsistent with the source cited on the entry
itself: Hirsch 2002 (Anal Biochem 308:343) is the standard reference for
**~10⁻¹¹ M** (ΔG −15.0), and Green 1990 (Methods Enzymol 184:51) puts
desthiobiotin 10²–10⁴ fold weaker than biotin, which brackets it. The
value is what makes desthiobiotin useful as a reversible affinity tag:
tight enough to capture, loose enough to elute with free biotin.

**Corrected, and the consequences followed through rather than hidden.**

| | Before (wrong K_d) | After (corrected) |
|---|:-:|:-:|
| Spearman ρ | +0.80 | **+0.40** |
| MAE | 4.98 kcal/mol | **6.66 kcal/mol** |

The claim this set used to support — that docking ranks biotin-site
binders usably — does not hold. The class has been moved out of the
"reliable for ranking" table in the tutorial into "pose filter only",
alongside kinases.

**What replaced it is a stronger test.** `tests/uq/test_calibration_smoke.py`
asserted Spearman ≥ 0.8, which only passed because of the bad number.
It now pins the saturation itself: across four compounds spanning
**12.0 kcal/mol** of experimental affinity, the predictions span
**0.33 kcal/mol**, all between −7.1 and −7.4. That is Vina's tight-binder
saturation measured directly, it is reproducible, and it is the property
the reliability table's `do_not_trust_absolute` verdict rests on.

Downstream artifacts regenerated from the corrected data:
`benchmarks/dock/uq_coverage.json` and the `ultra_tight_binder` row of
`benchmarks/dock/reliability_table.yaml`.

## Docking: accuracy per target class

Measured absolute error of docking ΔG against experiment, so a
prediction can carry an accuracy estimate, not only seed scatter.
Monte-Carlo over seeds understated biotin's real error about 40-fold.

| Target class | n | MAE (kcal/mol) | Verdict |
|---|:-:|:-:|---|
| Trypsin-like S1 pocket (3PTB) | 6 | 0.91 | absolute ΔG usable |
| Kinase ATP site (1M17, EGFR) | 6 | 2.16 | rank-order only, and even that is weak (Spearman −0.49) |
| Ultra-tight binders (1STP) | 4 | 4.98 | absolute ΔG meaningless; Vina saturates |

n < 19 per class, so these are typical errors seen so far, not
guaranteed conformal bounds. Pooled split-conformal across classes
gives ± 11.6 kcal/mol, which is honest and useless; calibrate per class.

Reproduce: `python scripts/run_uq_coverage.py`; table in
`benchmarks/dock/reliability_table.yaml`.

## CYP3A4 site of metabolism (Layer 1.4)

| Path | Literature set (aspirin, midazolam, diazepam) |
|---|:-:|
| xTB C–H BDE ranking | 1/3 |
| + heme-accessibility re-rank (default) | 2/3 |

Systematic blind spot: N-dealkylation sites (diazepam and many
tertiary-amine substrates). The tool is advisory and is labelled so.

Reproduce: `cellsim som-validate benchmarks/quantum/cyp3a4_som_validation.yaml`.

## Alchemical FEP (experimental)

FreeSolv-12 hydration, 11 windows × 50 ps per leg, CPU, tag
`milestone-a-pilot-3` (report archived in
`archive/campaign1/freesolv_pilot3/`):

- 10 of 12 compounds converged (acetic acid and acetamide failed MBAR
  overlap on the hand-rolled sampler; replica exchange converges them).
- MAE 1.42 kcal/mol, Pearson 0.91, but only 1 of 10 residuals inside the
  MBAR error bar, and residuals grow with solute size: methane −0.2,
  ethane −1.6, propane −1.6, benzene −1.9, pyridine −3.1, toluene −3.7.
- Cause identified 2026-10-01: the solvated system is a packmol box run
  in NVT; nothing in `cellsim/fep` applies a barostat or equilibrates the
  density. This is a systematic error, not sampling noise.
- Binding ΔG (double decoupling) has never produced a sampled number on a
  real binder. PR #10 fixed 13 audited defects, including a ~10 kcal/mol
  sign error in the standard-state correction.

Status: kept as reference code with its smoke gates; not a product path.

## Legacy simulator (`OLD/`)

Builds headless with `cmake -S OLD -B build-old && cmake --build build-old`
(no vcpkg needed for the validators). Re-run 2026-10-01:

| Validator | Seed 0x1 | Seed 0xDEADBEEF |
|---|:-:|:-:|
| `cellsim_blind_benchmarks` (BM01–BM08) | 8/8 | 8/8 |
| `cellsim_p53_pulsing` (gates a–h) | PASS | — |
| `cellsim_cisplatin_arrest` | PASS | — |
| `cellsim_cisplatin_apoptosis` | PASS | — |

What those passes do and do not mean: the p53 axis responds, oscillates
(7 peaks / 24 h, period 3.03 h against 5.5 h published), arrests the
cycle (100 % G2 at 24 h, stronger than the literature's "dominant"), and
builds the PUMA/BAX/MOMP commitment signature. Cells do **not** die:
caspase-3 stays at 0.04 under saturated MOMP, so death counts are 0 under
cisplatin in every run. Units are "sim units" in several subsystems
(`OLD/UNITS.md`). Treat `OLD/` as a reference implementation of the
pathway topology and the citations, not as a validated predictor.

## Phase-1 cell engine (`cellsim/cell/engine.py`)

Twelve gates in `tests/cell/test_engine_smoke.py`, all passing, each
asserting a phenotype rather than "it ran" (the two not tabulated below
cover cell-to-cell variability), plus five in
`tests/cell/test_stream_smoke.py` for the UI stream, including that a
6 h pulse washed out regrows (×4.0 at 72 h) where the same dose held
continuously kills the colony:

| Gate | Result |
|---|---|
| Untreated colony doubles at the line's quoted rate | A549 grows ×10.3 in 72 h vs ×9.7 expected for a 22 h doubling |
| Asynchronous phase fractions are physiological | G1 0.67 / S 0.18 / G2 0.13 / M 0.03 |
| p53 pulses under sustained damage, Purvis-like period | 4.9 h (published 5.5 h; the C++ prototype gave 3.0 h) |
| No pulsing without damage (specificity) | basal p53 stays below 0.12 |
| Damage *executes* death, not just commitment | caspase-3 reaches 1.0, cells are removed |
| p53-mutant line resists damage-induced death | HT-29 deaths ≪ A549 deaths at the same damage |
| Dose-response is monotone with a finite IC50 | holds |
| Fixed seed is reproducible | exact |
| State stays finite and bounded | holds |
| Cycle scale tracks doubling time | HCT116 (18 h) > MCF7 (29 h) |

The death-execution gate is the one that matters most: it is the defect
that made the C++ prototype unusable (MOMP saturated while caspase-3 sat
clamped at 0.04, so no cell ever died). The port fixes it by putting
Smac in excess over XIAP, as measured (Rehm 2006).

**Curve shape is not yet realistic.** Each simulated cell carries
log-normal variability (σ 0.3) in its Bcl-2 reserve and drug uptake,
the documented source of fractional killing (Spencer 2009). That widens
the dose window over which viability falls from 85 % to 15 % from 1.2×
(identical cells) to about 2×, measured on a 10-point-per-decade grid.
A Hill slope of 1 would give about 32×. So the simulated curves are far
steeper than typical measured ones, and only IC50s are compared below.

## Phase-1 GDSC validation (`scripts/validate_gdsc.py`)

Protocol: each drug's single `fit_target` constant is fitted by
bisection on the clean TP53 wild-type lines A549 and MCF7 (geometric
mean of the per-line fits) and applied unchanged to every line. HeLa and
both TP53-mutant lines are held out. Results in
`benchmarks/cell/gdsc_validation_results.csv` and `…_summary.json`.

Ten cell lines (five TP53 wild-type, five mutant), fitted on A549 and
MCF7 only, so eight lines are held out. The table is the first ten-line
run, before the apoptotic buffer and the efflux channel were built (see
below); the current engine's figures follow it.

| Drug | Fitted constant | Engine in span | Null in span | log10 RMSE engine (held-out) | log10 RMSE null |
|---|---|:-:|:-:|:-:|:-:|
| cisplatin | `k_damage_per_uM_h` = 0.000954 | 9/10 | 7/10 | 0.54 (n=5) | 0.47 |
| doxorubicin | `k_damage_per_uM_h` = 0.00794 | 8/10 | 8/10 | 0.45 (n=8) | 0.46 |
| paclitaxel | `partition` = 0.153 | 10/10 | 8/10 | 0.31 (n=8) | 0.35 |
| **held-out** | | **22/25** | 18/25 | | |

With the current engine (buffer and efflux on, library refitted
October 2026 to cisplatin 0.001181, doxorubicin 0.01749, paclitaxel
partition 0.212): held-out in span **20/25** against the null's 18/25;
held-out log10 RMSE 0.65 / 0.43 / 0.32 against the null's
0.47 / 0.46 / 0.35; paired test 10 : 11, p = 1.0
(`benchmarks/cell/gdsc_validation_summary.json`).

### Exit gate: NOT met

*Re-measured at ten drugs in Phase 3 (below): 34 wins to 39 losses over
91 held-out predictions, p = 0.64. The conclusion holds, and the larger
panel also found the one exception — a drug whose target is part of the
modelled mechanism, where the engine does call lines correctly.*

Paired per-line comparison against the null on held-out lines: the
engine is closer on **10**, the null is closer on **11**, two-sided sign
test **p = 1.0**. The engine does not yet tell cell lines apart better
than predicting one number for all of them.

**This corrects an earlier claim in this file.** The gate first used to
read "beat the null's held-out RMSE for 2 of 3 drugs", and on a
five-line set it passed. That was wrong twice over: the set was too
small, and the margin that passed it was 0.45 against 0.46, far inside
the noise on eight points. The gate is now the paired sign test, which
cannot pass on a difference that small, and the whole set is ten lines.

**What is solid.** The engine puts each drug's potency on the right
scale from one fitted constant: 22 of 25 held-out predictions landed
inside GDSC's replicate span against the null's 18 in this run (20 of 25
with the current engine), and A549 cisplatin predicts 9.0 µM against
9.8 measured. Paclitaxel is the best case,
10/10 in span and the lowest error of the three.

**What is missing, now measured.** Line-to-line variation. The only
line-specific inputs are doubling time and TP53 status, and the
reference data shows TP53 status only predicts sensitivity for one of
the three drugs (table above). Doubling times themselves disagree
between sources by up to 1.6× for these lines. So the engine has almost
no information with which to separate lines, and it shows.

**The p53-independent death route: helps, but not enough.** ATM/c-Abl →
p73 → PUMA is the documented route by which TP53-mutant cells still die
of DNA damage (Gong 1999 Nature 399:806; Agami 1999 Nature 399:809).
Added as `Params.p73_gain`, default 0.1. Tested leave-one-out across all
five mutant lines (`scripts/experiment_p73.py`): the gain is chosen on
the other four and scored on the one left out, so the test line never
influences it. Held-out log10 error per mutant:

| Held out | p53-only | With route | Null |
|---|:-:|:-:|:-:|
| HT-29 | 0.388 | **0.332** | 0.410 |
| MDA-MB-231 | 0.524 | **0.325** | 0.405 |
| MDA-MB-468 | 0.980 | 0.337 | **0.321** |
| MIA-PaCa-2 | 1.100 | 0.590 | **0.422** |
| T47D | **0.215** | 0.716 | 0.739 |

The route is a large improvement over the p53-only model on four of five
mutants, where that model is badly wrong (off by up to 12×). But it
beats the null on only two, and it actively **hurts T47D**, the one
mutant that is genuinely drug-resistant: the route kills it when the
data says it survives.

That is the whole problem in one row. Treating p53 loss as "these cells
still die" fixes the lines that are sensitive and breaks the line that
is resistant, because what separates them is not p53 at all. The earlier
two-mutant version of this table looked like a clean win and was simply
too small.

### Why it cannot discriminate lines: measured, not guessed

The sign-test tie is not a calibration problem. The engine has no
channel through which cell lines can differ enough to matter. Measured
by sweeping each per-line field across its full observed range and
recording how far the predicted IC50 moves:

| Per-line field | swept over | cisplatin | doxorubicin | paclitaxel |
|---|---|:-:|:-:|:-:|
| doubling time | 18–38 h | 1.0× | 1.1× | 1.0× |
| Bcl-2 reserve | 0.5–4× | 1.1× | 1.2× | 1.0× |
| DNA repair rate | t½ 32 h–2 h | **9.1×** | **7.4×** | 1.0× |
| *observed between lines* | | *10.4×* | *51.3×* | *8.6×* |

Three things follow.

1. **Doubling time, the engine's main line-specific input, has no
   leverage at all.** Sweeping the entire range of the ten lines moves
   the IC50 by 1.0×. Since TP53 status is the only other per-line input
   and it is binary, the engine can produce at most two distinct
   answers per drug. Tying with a single constant is the expected
   result, not a surprise.
2. **DNA repair rate is the one channel that works**, and it nearly
   spans what cisplatin needs (9.1× against 10.4×). This is where
   per-line information should enter for DNA-damage drugs — ERCC1 and
   the rest of the NER machinery, not proliferation rate.
3. **Paclitaxel has no line-specific channel whatsoever.** Every field
   gives 1.0×, because its IC50 is set by the tubulin-occupancy arrest
   threshold, which is a property of the drug and not of the line. To
   make it line-specific the engine needs per-line efflux (ABCB1 acting
   on intracellular drug), which it does not currently model.

### What per-line data actually predicts, at a sample size that can tell

At ten lines nothing could be validated. Repeating the test on every
GDSC line DepMap has expression for — **156 to 683 lines** rather than
6 to 10 — with genes **pre-specified by mechanism** before any
correlation was computed, and Bonferroni-corrected within each drug
(`scripts/experiment_expression_scale.py`):

| Drug | n lines | Survives correction, in the expected direction |
|---|:-:|---|
| doxorubicin | 683 | **BCL2L1** ρ +0.19 (p 5e-06), **ABCB1** ρ +0.17 (p 9e-05) |
| paclitaxel | 427 | **BCL2L1** ρ +0.31 (p 9e-10), **ABCB1** ρ +0.23 (p 1e-05), **TUBB3** ρ +0.22 (p 3e-05) |
| cisplatin | 156 | **nothing** |

So real, mechanistically sensible signal does exist once the panel is
large enough: the anti-apoptotic set-point (Bcl-xL), P-glycoprotein
efflux, and class III β-tubulin for taxanes, which is the textbook
taxane-resistance marker. Cisplatin gives nothing at all — and notably
**the repair genes do not predict cisplatin sensitivity** (ERCC1
ρ −0.03, ERCC2 +0.03, XRCC1 −0.06), consistent with ERCC1's long record
of failing to replicate as a platinum biomarker.

### The mismatch that explains everything

Put the two measurements side by side. Translating each significant
correlation into the IC50 change implied across its 10th-to-90th
percentile expression range:

| Channel | What the data implies | What the engine delivers |
|---|:-:|:-:|
| anti-apoptotic set-point (Bcl-xL) | 1.8–2.3× | **1.1×** |
| efflux (ABCB1) | 1.5–2.1× | **not modelled** |
| class III β-tubulin | 1.9× | not modelled |
| DNA repair | **no signal** (cisplatin) | 7–9× |
| proliferation rate | not tested | 1.0× |

The engine has its leverage exactly where the data has no signal, and
has little or no leverage exactly where the data does. That is the
whole failure, stated in one table:

1. **Repair rate** is the engine's only strong channel, and repair
   genes do not predict cisplatin sensitivity.
2. **The apoptotic set-point** is the strongest real signal and the
   engine already has the field for it (`CellLine.bcl2_level`), but its
   MOMP switch is so sharp that changing that field 8-fold moves the
   IC50 by 1.1×. The channel exists and is throttled.
3. **Efflux** is the second strongest real signal and the engine does
   not model it at all.

**A ceiling worth knowing before building.** These effects are modest:
roughly 1.5–2.3× each. Even used perfectly and independently they
compose to maybe 3–5×, against the 8–51× by which these lines actually
differ. So wiring them in should make the engine beat a constant, but
it will not close the gap — most line-to-line variation in drug
response is not captured by the canonical markers. That is a property
of the biology, not of this implementation.

### Step 1 built: the apoptotic set-point is no longer throttled

The first item was acted on. The Bcl-2 reserve was modelled as a **rate
penalty** on Bax activation, which is why it barely mattered: a slower
climb still crosses a fixed threshold inside 72 h, so the reserve
changed *when* a cell died, not *whether*. It is now a **stoichiometric
buffer** — only Bax in excess of the reserve can form pores — which is
the mitochondrial-priming picture (Certo et al. 2006 Cancer Cell 9:351).

| | Before | After |
|---|:-:|:-:|
| IC50 span over `bcl2_level` 0.5–4×, doxorubicin | 1.2× | **1.7×** (data implies 1.8×) |
| doxorubicin held-out log10 RMSE vs null | 0.45 / 0.46 | **0.40** / 0.46 |
| paclitaxel held-out RMSE vs null | 0.31 / 0.35 | 0.30 / 0.35 |
| cisplatin held-out RMSE vs null | 0.54 / 0.47 | **0.65** / 0.47 (worse) |
| paired sign test (engine : null) | 10 : 11 | **11 : 10**, p = 1.0 |

**Predicted and confirmed.** The ceiling above said wiring these in
"should make the engine beat a constant, but will not close the gap".
That is exactly what happened: the paired count flipped from losing to
winning, and stayed far from significance. A 1.7× channel cannot resolve
an 8–51× spread.

It also improved the model's honesty elsewhere. Treating the reserve as
a buffer restored **partial** resistance in TP53-mutant lines — survival
0.08 against wild-type 0.00 under sustained damage, where the previous
version killed mutants outright. That was flagged above as stronger than
the literature supports, and is now pinned by a test.

Cisplatin got worse, which fits: it is the drug with no Bcl-xL signal,
so the added mechanism only adds noise there.

### Step 2 built: a per-line efflux channel

ABCB1 was the remaining signal the engine could not express, and the
only possible route to line-specific behaviour for paclitaxel, whose
IC50 is set by a tubulin-occupancy threshold belonging to the drug. A
saturable P-glycoprotein pump now acts on intracellular drug for
substrate drugs only, scaled by the line's ABCB1.

`CellLine.efflux_level` is **derived, not fitted**: 2^(ABCB1 log2(TPM+1)
− panel median) from DepMap 24Q4, a fixed formula with no free
parameter, so held-out lines stay out of sample. At the pump strength
the measured ABCB1 range implies, leverage is 1.9× for doxorubicin and
2.0× for paclitaxel — matching the data's 1.5–2.1× — and exactly 1.0×
for cisplatin, which is not a substrate.

**It did not improve prediction.** The paired count moved from 11 : 10
to 10 : 11, a one-line change at p = 1.0 — noise, not degradation, but
certainly not the improvement the population statistics suggested. The
high-ABCB1 lines are the tell: HeLa carries ~10× the panel's median
ABCB1 yet sits mid-pack for doxorubicin sensitivity, so the engine now
predicts it too resistant.

### The ceiling, measured: these markers cannot support per-line prediction

That failure is not an implementation flaw, and the point generalises.
Taking all eleven candidate markers together and fitting them to IC50
by least squares — far more freedom than any mechanistic model gets —
across every GDSC line DepMap covers:

| Drug | n lines | best-case R² | spread it must explain |
|---|:-:|:-:|:-:|
| cisplatin | 156 | **0.09** | 4× |
| doxorubicin | 683 | **0.15** | 10× |
| paclitaxel | 427 | **0.20** | 7× |

**80–92 % of why one cell line is more sensitive than another is not
captured by the canonical markers at all.** A correlation of ρ ≈ 0.2 is
real at n = 683 and nearly useless for an individual line: it is a
population trend, not a predictor. So a mechanistic model restricted to
these inputs cannot pass the paired test, however faithfully it is
built — and building it more faithfully, as the efflux channel was,
does not help.

This is the honest answer to the question Phase 1 was posed to settle,
and it redirects the project rather than ending it. **Predicting which
cell line is more sensitive is not a goal this engine can reach with
available inputs.** What it does well is the other thing: reasoning
about one system over time — dose schedules and wash-out, resistance
emerging under selection, the timing of combinations, the difference
between arrest and kill. Those depend on mechanism the engine has and
on no per-line marker at all.

**Still to do, in order.** Revised by the measurement above:

1. **Re-aim Phase 2 at within-line dynamics**, which is where the
   engine's mechanism is load-bearing and where the validation data is
   not a ten-line IC50 table: schedule dependence, wash-out recovery,
   resistance under selection, combination timing.
2. **Keep the per-line channels** — apoptotic set-point and efflux are
   correct biology and are needed for drug-level reasoning — but stop
   treating per-line IC50 as the headline metric.
3. **If per-line prediction is revisited**, it needs inputs beyond
   canonical markers (pathway state, baseline apoptotic priming by
   BH3 profiling, or a much larger panel), not more mechanism.
4. A shallower death response (curve shape above).

Re-run `cellsim validate-gdsc` after each; the gate is the paired sign
test, not the RMSE table.

## Phase 2: audited before building on it (October 2026)

The library's three fitted potency constants predated the apoptotic
buffer and the efflux channel, although its docstring said re-running
the fit reproduced them within 1 %. Refitting under the current defaults
moved them 1.2–2.2× (cisplatin k_damage 0.000951 → 0.001181,
doxorubicin 0.00794 → 0.01749, paclitaxel partition 0.153 → 0.212).
That silently moved every Phase-2 dose relative to its drug's potency,
and re-running the experiments showed which headlines were properties of
the mechanism and which were properties of one dose or one random seed.

Two changes make that class of error structural rather than lucky.
Experiments now state doses as multiples of each drug's own engine IC50,
computed at run time (`cellsim.cell.engine.ic50`), so a refit cannot
change what a design means. Tests pin the potency their doses were
designed at; cisplatin and paclitaxel enter the engine only as potency ×
concentration, so for them this is exactly a rescaling of the dose.

| Claim as first published | Status | What re-measurement shows |
|---|---|---|
| Pulsing selects a pharmacokinetic escape (low uptake), not an apoptotic one | **Retracted: backwards** | It rested on ONE surviving lineage for cisplatin and two for doxorubicin. Measured exactly on a trait grid, continuous exposure selects on accumulation and pulsing shifts selection toward the apoptotic reserve, for all three drugs. |
| The sequence effect is a readout artefact (0.93 → 0.98 → 0.99) | **Qualified** | A single-seed reading; one 128-cell run carries ±0.02–0.03 on that ratio. Over five seeds there is no ordering effect with both drugs at their IC50, and ~12 % at 72 h with cisplatin at 1.25× its IC50, shrinking to ~3 % by 120–144 h. |
| Exposure-time ranking 242× / 17× / 7.6× (3 h ÷ 72 h survival) | **Replaced** | The ratio depends on the dose chosen, and once no drug acts within 3 h it reduces to 1/S(72 h). The iso-effect curve below replaces it; the ranking holds. |
| Antagonism 1.32–1.50× | Confirmed | 1.29–1.35 (±0.08 over seeds) with both drugs at their IC50, rising to 1.6–1.9 with more cisplatin. |

## Phase 2: schedule dependence (within-line, no marker needed)

Everything here is one cell line compared against itself under
different schedules, so the per-line ceiling above does not apply.

**The measurement is the iso-effect curve**: for each exposure time T,
the concentration C50(T) that halves the 72 h cell count when the drug
is washed out at T (A549; `scripts/experiment_schedule.py`). If killing
were governed by AUC, C50 × T would be constant and log C50 against
log T would have slope −1.

| Drug (72 h IC50) | 1 h | 3 h | 6 h | 12 h | 24 h | 72 h | slope, 1–6 h | shortest exposure reaching half kill |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| cisplatin (10.7 µM) | 326 | 120 | 59.1 | 30.6 | 16.3 | 10.4 | **−0.95** | 1 h |
| doxorubicin (0.063 µM) | 2.06 | 0.68 | 0.34 | 0.17 | 0.089 | 0.060 | **−1.01** | 1 h |
| paclitaxel (0.030 µM) | none | none | none | none | 0.037 | 0.030 | — | **24 h** |

(C50 in µM; "none" = no concentration reaches half kill.)

**Cisplatin and doxorubicin are AUC-governed over short exposures**:
C50 × T stays at 326–367 µM·h for cisplatin from 1 to 12 h and at
2.0–2.1 µM·h for doxorubicin from 1 to 24 h, then rises as repair
outpaces slow accumulation. Ozawa et al. measured exactly this C × T
law for cisplatin (WiDr cells; Cancer Res 1989 49:3823) and time
dependence instead for phase-specific agents.

**Paclitaxel is time-governed**: no exposure shorter than 24 h halves
the colony at ANY concentration, because only cells that attempt
mitosis while the drug is present are affected. At a fixed AUC both
extremes therefore fail, for different reasons — too brief and few
cells reach mitosis, too dilute and tubulin occupancy never crosses the
arrest threshold — and the best split is interior (C50 × T is 0.88 µM·h
at 24 h against 2.1 at 72 h). That is the documented clinical finding
that taxane efficacy follows time above a threshold concentration
rather than AUC (Gianni 1995 J Clin Oncol 13:180; Huizing 1993 J Clin
Oncol 11:2127).

Against published exposure-duration studies, matched by assay:

| Study | Published | Engine | Verdict |
|---|---|---|---|
| Ozawa 1989, cisplatin, WiDr | C × T governs (slope −1) | slope −0.95 over 1–6 h | agrees |
| Liebmann 1993 (Br J Cancer 68:1104), 8 lines, clonogenic | 24 h kill saturates above ~50 nM | 0.46 / 0.43 / 0.40 / 0.49 surviving at 0.05 / 0.2 / 1 / 5 µM | agrees |
| Liebmann 1993 | 24 → 72 h exposure raises cytotoxicity 5–200× | C50 0.037 → 0.030 µM (1.2×) | **miss** |
| Georgiadis 1997 (Clin Cancer Res 3:449), 14 NSCLC lines, MTT, read at 120 h: 3 / 24 / 120 h exposure | IC50 > 32 / 9.4 / 0.027 µM | none / 0.036 / 0.029 µM (A549) | 3 h agrees; 120 h agrees but is inherited from the GDSC fit; **24 h misses by 260×** |

**The miss, stated plainly.** Both studies say a 24 h paclitaxel
exposure is far weaker than a 72–120 h one; the engine says it is
nearly as good. The engine kills almost every cell that reaches
mitosis inside the 24 h window. The likely causes are that its cells
cycle almost in lockstep (real cycle times vary by 20–30 % and include
slow-cycling cells that a 24 h window misses) and that cells slipping
out of mitotic arrest re-enter the cycle instead of arresting as
tetraploids. Drug retention after washout cannot be the cause, because
it would make short exposures stronger, not weaker. Open; recorded so
the schedule result is not over-read.

Gated by `tests/cell/test_schedule_smoke.py` (5 gates, 38 s): the C × T
slope for cisplatin, the 6 h paclitaxel pulse that cannot reach half
kill, the sub-threshold thin split, the interior optimum, and
continuous-beats-washout.

## Phase 2: selection, and what the schedule selects *for*

Each representative cell is a lineage carrying its own anti-apoptotic
reserve (Bcl-2 multiplier) and drug accumulation (uptake multiplier),
and a surviving lineage keeps both as it doubles, so selection emerges
without being modelled.

**Measured exactly, not sampled.** `cellsim.cell.selection` runs a
13 × 13 grid of trait values (six cells per point, random cycle phases)
through each schedule, records survival as a function of trait, and
integrates that map against the population's log-normal trait
distribution. That gives the survivors' expected means however rare
survival is; a logistic fit gives the boundary's steepness along each
axis, and |b_uptake / b_reserve| says which trait decides survival.

Two weeks, A549, matched total exposure: continuous at the drug's own
72 h IC50 versus 3× IC50 for 24 h then 48 h drug-free
(`scripts/experiment_resistance.py`):

| Drug | Arm | Expected survival | Survivors' reserve | Survivors' uptake | \|b_u / b_r\| |
|---|---|:-:|:-:|:-:|:-:|
| cisplatin | continuous | 0.235 | 1.01 | 0.64 | 5.2 |
| | pulsed 1 d / 2 d | 0.0016 | **1.88** | 0.50 | **1.35** |
| doxorubicin | continuous | 0.070 | 1.21 | 0.56 | 1.8 |
| | pulsed 1 d / 2 d | 0.0014 | **2.45** | 0.76 | **0.91** |
| paclitaxel | continuous | 0.476 | 1.00 | 0.75 | uptake only |
| | pulsed 1 d / 2 d | 0.45 | 1.14 | 0.98 | **0.39** |

(1.00 = no selection.)

**Continuous exposure selects on drug accumulation.** A sustained dose
sets a damage steady state, and whether it stays under threshold depends
on how much drug a cell takes up; for cisplatin the reserve barely moves
(1.01). **Pulsing shifts selection toward the apoptotic reserve.** A
brief high dose pushes nearly every cell over threshold for a while, and
what decides survival is whether the reserve can absorb the transient —
it can buffer a pulse but not a siege.

**Pulsing kills far more of a concentration-driven drug**: 143× fewer
survivors for cisplatin and 49× for doxorubicin at matched exposure,
but 1.06× for paclitaxel, which wastes a compressed exposure because it
needs cells to transit mitosis while it is present — the same mechanism
as the schedule result, measured independently.

**A prediction an experimental lab can test.** Pulsed and continuous
selection are the two standard ways resistant lines are derived
(McDermott et al. 2014 Front Oncol 4:40). The engine predicts that the
continuously derived line resists mainly by accumulating less drug
(e.g. raised efflux), and the pulse-derived one mainly by a larger
apoptotic reserve (lower mitochondrial priming on BH3 profiling).

**What this cannot show.** Heterogeneity is drawn once, so this is
selection from a pre-existing tail, not evolution of new resistance;
lineages grow without a carrying capacity, so they never compete, and
adaptive-therapy questions are out of reach until the spatial dish.

Gated by `tests/cell/test_resistance_selection_smoke.py` (5 gates,
18 s), on a 7 × 7 trait grid.

## Phase 2: combinations, and a timing confound

The engine carries one intracellular concentration per drug, so a
combination is a real two-drug run. Every arm gives each drug for 36 h
at the same concentration in the same line; only the timing differs.
Five seeds × 256 cells, mean ± sd over seeds
(`scripts/experiment_combination.py`):

| | Both at 1× IC50 | Cisplatin 1.25× IC50 |
|---|:-:|:-:|
| together ÷ independent expectation | 1.29 ± 0.08 | 1.64 ± 0.07 |
| cisplatin → paclitaxel ÷ independent | 1.34 ± 0.08 | 1.68 ± 0.08 |
| paclitaxel → cisplatin ÷ independent | 1.35 ± 0.06 | 1.92 ± 0.17 |

**Antagonism is robust.** Every arrangement kills less than independent
action predicts, and more so with more cisplatin, because the platinum
arrests the cycle through p53 → p21 and removes the cells paclitaxel
needs in mitosis — the textbook interaction of a cytostatic with a
phase-specific agent. Nothing was fitted to produce it.

**The ordering effect is small, dose-dependent, and mostly timing.**
Cisplatin-first survival divided by paclitaxel-first survival:

| Readout | 72 h | 96 h | 120 h | 144 h |
|---|:-:|:-:|:-:|:-:|
| both at 1× IC50 | 0.99 ± 0.02 | 1.00 ± 0.02 | 1.00 ± 0.02 | 1.01 ± 0.02 |
| cisplatin 1.25× IC50 | **0.88 ± 0.06** | 0.93 ± 0.05 | 0.96 ± 0.03 | 0.97 ± 0.03 |

With both drugs at their IC50 there is no ordering effect at all. With
more cisplatin the platinum-first order looks ~12 % better at 72 h, and
most of that disappears once the slow drug given second has had time to
act. So the engine does not reproduce a biological platinum/taxane
sequence dependence; what it shows is how a fixed-endpoint assay can
confound "the order matters" with "the second drug had less time".

Gated by `tests/cell/test_combination_timing_smoke.py` (4 gates, 41 s):
two-drug state, antagonism, no ordering effect at matched IC50 (3-seed
mean), and the fade of the above-IC50 ordering effect with readout.

## Phase 2: the dish — cells in space (`cellsim/cell/dish.py`)

The well-mixed engine says how many cells survive; the dish says which
ones and where. Each lattice site holds at most one cell and each cell
carries a full engine state row, so the validated biology runs per cell
unchanged. Space adds contact inhibition (a dividing cell must push its
neighbours to free space within a set reach, or its G1 growth signal
drops to zero), diffusing oxygen and drug consumed by the cells they
pass, hypoxic quiescence, necrosis, and clones that inherit their
founder's traits. Two geometries: a monolayer under stirred medium and a
3-D spheroid whose fields are solved radially, which is exact for a
symmetric aggregate and is how the spheroid literature models it.

Every environmental constant is an input with a source: oxygen
diffusivity, medium pressure and consumption from DLD-1 spheroids
(Grimes et al. 2014 J R Soc Interface 11:20131124: 2e-9 m²/s, 100 mmHg,
22.1 mmHg/s in packed tissue); oxygen-dependent proliferation from
PhysiCell's documented defaults (Ghaffarizadeh et al. 2018); apoptotic
clearance 8.6 h (Macklin et al. 2012); small-molecule interstitial
diffusivity 1e-6 cm²/s (Nugent & Jain 1984). How much drug a cell
sequesters, which sets how slowly drug penetrates tissue, is known only
to an order of magnitude and is flagged as such per drug.

**Two engine changes it needed, both inert by default.** A per-cell G1
growth signal: lowering cyclin D synthesis alone does not arrest this
cycle, whose Rb/E2F/cyclin E feedback restarts itself, so the signal
scales progression through G1 only — cells past the restriction point
finish their cycle and their daughters wait (gs = 0 leaves exactly the
56 % of cells that were past G1 dividing once and nothing dividing
twice). And optional cycle-time variability (`Params.cycle_cv`). With
neither in use, every earlier result is bit-identical.

**Numerics against closed forms** (`tests/cell/test_dish_smoke.py`). The
oxygen solver reproduces the zero-order sphere profile to within 1 % of
the medium value and Grimes's 233 µm diffusion limit; with first-order
uptake the slab solver reproduces `cellsim.cell.tissue`'s analytic
profile to 0.5 %; the transient drug solver relaxes to the steady one.

**A sparse monolayer: the Cell Tracking Challenge HeLa movie.** Both
46 h sequences re-run as closed fields of the same size, seeded with the
same counts (five seeds each; `scripts/validate_dish.py`, targets from
the curated tracks by `scripts/ctc_reference.py`):

| | sequence 01 | sequence 02 |
|---|:-:|:-:|
| cells at 46 h, measured | 43 → **137** | 125 → **363** |
| engine (cycle CV 0.25) | **140** (trajectory error 12 %) | **377** (15 %) |
| complete cycles, measured | 20.2 h, CV 0.22 | 18.5 h, CV 0.33 |
| engine | 26.4 h, CV 0.23 | 27.3 h, CV 0.22 |
| still undivided at 30 h, measured | 0.51 | 0.71 |
| engine | 0.41 | 0.45 |

The colony grows right; the single-cell picture does not. The movie's
cells are a mix of fast cyclers (complete cycles ~19 h) and a large
fraction that does not divide within the window — the proliferation /
quiescence decision at mitotic exit (Spencer et al. 2013 Cell 155:369),
which the engine does not model. Open. With every cell on one clock
(cycle CV 0) no cell divides within 30 h of birth at all, so measured
variability is clearly needed.

This also caught a provenance error. The library gave HeLa a 20 h
doubling time, read off this movie's COMPLETE cycles, which a 46 h
window biases short. Cellosaurus gives 1.3 days (PubMed 29156801; DSMZ
~48 h) and the movie's own counts double every 27–31 h, so it is now
31 h. HeLa's held-out GDSC predictions moved 3–7 % and still pass.

**A spheroid: DLD-1, calibrated, not yet independently validated.**
The dish has one free spatial constant, how far a dividing cell can
push its neighbours (`DishParams.push_sites`). Grimes's DLD-1 spheroids
grow ~15 µm/day in radius after their core turns anoxic (233 µm to
~400 µm between days 6 and 17), and from 3000 cells over 12 days:

| mechanical reach | necrosis starts | radial growth after that | Grimes |
|---|:-:|:-:|:-:|
| **1 site** | day 6 (R 206 µm) | **15.3 µm/day** | **~15 µm/day** |
| 2 sites | day 4 (R 211 µm) | 28.2 µm/day | |
| 3 sites | day 4 (R 245 µm) | 36.3 µm/day | |
| unlimited | day 3 (R 277 µm) | 58.1 µm/day | |

A reach of one site matches, and leaves ~80 % of viable cells
quiescent with proliferation confined to the outer layer, consistent
with Ki-67 rims; solid stress is known to limit spheroid proliferation
this way (Helmlinger et al. 1997 Nat Biotechnol 15:778). Without a
mechanical limit the spheroid grows four times too fast.

**A miss in the necrotic core.** Feeding the model's (radius, necrotic
radius) pairs through Grimes's own anoxic-core relation should return
their 233 ± 22 µm diffusion limit if the model's core tracks anoxia.
With PhysiCell's generic necrosis threshold (5 mmHg) it returns 155 µm;
moving the threshold to near-anoxia, which is how Grimes define the
core, gives 173 µm. The core is still too large. The likely
cause is an overshoot built into its formation: the first cells die
under the profile of a fully consuming sphere, and once the dead core
stops consuming, oxygen reaches deeper than the dead region, so the core
is bigger than the anoxic region it now sits in. A surface-growth
lattice also roughens the surface, so the equivalent radius understates
how far oxygen must travel. Open.

**What the dish does not do yet.** One cell per site, so no compression
or multilayering beyond a fixed reach; necrotic debris never lyses; no
glucose, lactate or pH; drug sequestration is order-of-magnitude; and
spheroid drug response has not been checked against data (multicellular
resistance to a 1 h paclitaxel exposure, Nicholson et al. 1997 Eur J
Cancer 33:1291, cannot be used: the engine keeps no paclitaxel in a cell
after wash-out, and kills nothing with exposures under 24 h).

Gated by `tests/cell/test_dish_smoke.py` (7 gates, 20 s).
## Phase 3: ten drugs, and the one place the engine does tell lines apart

The panel is ten drugs across eleven lines — 91 held-out predictions
against 25 in Phase 1 — chosen so that every mechanism the engine can
express is covered and each has GDSC screens on our lines. Each carries
exactly one fitted constant, fitted on the clean TP53 wild-type lines and
applied unchanged everywhere else (`scripts/validate_gdsc.py`).

| Drug | Mechanism | Fitted | In GDSC span | Sensitive/resistant call | Spearman |
|---|---|---|:-:|:-:|:-:|
| cisplatin | DNA adduct | `k_damage` 0.00118 | 8/11 | 6/11 | −0.03 |
| doxorubicin | TopII | `k_damage` 0.0175 | 9/11 | 11/11 | 0.58 |
| etoposide | TopII | `k_damage` 0.00958 | 6/10 | 10/10 | 0.41 |
| sn-38 | TOP1, S phase | `k_damage` 2.07 | 8/11 | 11/11 | −0.19 |
| gemcitabine | antimetabolite | `k_damage` 1.23 | 10/11 | 10/11 | 0.72 |
| 5-fluorouracil | antimetabolite | `k_damage` 0.00313 | 9/11 | 6/11 | 0.94 |
| paclitaxel | tubulin | `partition` 0.212 | 9/11 | 10/11 | −0.06 |
| docetaxel | tubulin | `partition` 0.274 | 7/11 | 11/11 | 0.03 |
| vinorelbine | tubulin | `partition` 0.302 | 7/11 | 10/11 | 0.44 |
| **nutlin-3a** | **MDM2** | `partition` 0.013 | **11/11** | **11/11** (null 5/11) | −0.30 |

**The Phase-1 conclusion survives the larger panel, and is now much
better measured.** On held-out lines the engine is closer than a
constant-IC50 null on 34 and further on 39 (sign test p = 0.64), and
lands inside GDSC's replicate span 67 times against the null's 68. Ten
drugs and 91 predictions say the same thing three drugs and 25 said:
**the engine reproduces each drug's potency scale and does not rank cell
lines better than a single number.** That is now a result rather than a
small-sample suspicion.

**But the panel also found the exception, and it is the informative
one.** Nutlin-3a is the only drug whose line-to-line differences are
*mechanism* rather than an unmeasured marker: it blocks MDM2, so p53
accumulates without any DNA damage, and whether that kills depends on
whether functional p53 is there to accumulate. The engine gets all
eleven lines right, nine of them held out, where the null gets five:

| Line | TP53 | Engine | GDSC | |
|---|---|---|---|---|
| A549, MCF7, HCT116, HT-1080, U-2-OS | wild-type | sensitive | sensitive | ✓ |
| HT-29, MDA-MB-231, MDA-MB-468, MIA-PaCa-2, T47D | mutant | resistant | resistant | ✓ |
| **HeLa** | **wild-type** | **resistant** | **resistant** | ✓ |

HeLa is the sharp one. Its TP53 is wild-type by sequence, so a model
reading the genotype would call it sensitive and be wrong. HPV18 E6 hands
p53 to the E6AP ligase instead of MDM2, and in cervical cancer cells that
switch is complete (Hengstermann et al. 2001 PNAS 98:1218), so blocking
MDM2 cannot rescue a protein that is no longer being degraded by MDM2.
That route is in the engine as `CellLine.p53_mdm2_independent_deg_per_h`,
and switching it off makes the engine predict HeLa sensitive — the
specificity check is in the test suite.

Nothing about this was fitted to the p53 split: nutlin's one fitted
constant is its partition coefficient, fitted on two wild-type lines.

**So the honest statement is narrower and more useful than either
extreme.** The engine cannot rank lines for cytotoxic drugs, because the
differences there live in markers it cannot see. It *can* call which
lines respond to a drug whose target is part of the modelled mechanism,
including a line whose genotype points the wrong way. That is the case
for building more mechanism rather than more markers.

Gated by `tests/cell/test_mdm2_p53_smoke.py` (4 gates, 8 s).

## Phase 3: how wrong is a predicted IC50?

A prediction without an interval is not usable, and an interval that has
never been scored is not an interval
(`scripts/calibrate_cell_uncertainty.py`). The width comes from the
spread of held-out log10 errors; the test is **leave-one-drug-out**, so
each drug is scored with a width built only from the other nine. Scoring
a width on the errors that built it would read as well calibrated however
badly it generalised.

| Nominal | Engine interval | Engine coverage | Null interval | Null coverage |
|---|:-:|:-:|:-:|:-:|
| 68 % | ±2.6× | 64 % | ±2.7× | 70 % |
| 90 % | ±4.8× | 85 % | ±5.0× | 90 % |
| 95 % | ±6.5× | **95 %** | ±6.8× | 93 % |

Every level lands within 15 points of nominal, so the stated interval
means what it says. Two things to read with it. The intervals are wide —
a factor of five at 90 % — and that is the honest width given that most
line-to-line variation is not in the markers the engine can see; GDSC's
own replicate screens of the same line and drug disagree by a median
5.3×, which is the floor nothing can honestly go below. And the null's
intervals are just as narrow at the same coverage, which is the same
finding as the paired test, now in the uncertainty picture: the engine's
*point* predictions are not better than a constant per drug, so neither
are its error bars.

## Phase 3: resistance that evolves, not only resistance selected

Everything before this could select from the spread a population already
had: heterogeneity was drawn once, so every lineage that survived had
been resistant from the start. `DishParams.mutation_rate` lets a daughter
differ heritably from her mother, and the control that makes the claim
testable is a **homogeneous start** — identical cells, so selection has
nothing to act on (`scripts/experiment_evolution.py`, A549 + cisplatin,
seven days, fold against the untreated arm of the same condition):

| Start | Mutation | Schedule | Survivors | Accumulation | Reserve | IC50 |
|---|---|---|:-:|:-:|:-:|:-:|
| identical | off | continuous | **0** | — | — | — |
| identical | off | pulsed 3× | **0** | — | — | — |
| identical | on | continuous | 230 | **0.61** | 1.03 | **1.85×** |
| identical | on | pulsed 3× | **0** (10 mutations first) | — | — | — |
| spread 0.3 | off | continuous | 2535 | 0.68 | 1.05 | 1.29× |
| spread 0.3 | on | pulsed 3× | 27 | 0.56 | **1.75** | **2.18×** |

Four things follow, and the first is the control:

1. **Identical cells with no heritable change are eradicated.** No
   resistance is possible when there is neither standing variation nor
   mutation, which is what makes the rest of the table mean something.
2. **Resistance arises where nothing was there to select.** From
   identical cells, continuous exposure leaves survivors that accumulate
   40 % less drug and need 1.85× the dose.
3. **A hard intermittent schedule can win the race.** The same total
   exposure given as brief 3× pulses eradicates the identical population
   after only ten mutations have occurred — resistance never gets
   started.
4. **Where standing variation exists, pulsing leaves the most resistant
   remnant**, and selects on the apoptotic reserve (1.75) rather than on
   accumulation. That is what the corrected Phase-2 selection measurement
   found by a completely different method.

Fold-resistance reaches 1.3–2.2×, against the two- to eight-fold that
derived lines actually reach (McDermott et al. 2014 Front Oncol 4:40), so
the assumed mutation rate is if anything conservative. It is a modelling
choice, not a measurement: resistance here arises from many routes at
once whose combined rate per division is not a published constant, and
the script reports what a given rate produces rather than claiming one.

Gated by `tests/cell/test_evolution_smoke.py` (4 gates, 18 s), including
the control that identical cells without mutation cannot drift.

## Phase 4: reading a lab's own plate (`cellsim/plate.py`)

Everything above uses our cell lines and our fitted constants. A lab has
its own line, its own reader and its own drug, and until their numbers
can enter the simulator none of this is a tool for them. This is the
first half of that door: get a plate in, normalised, fitted, with an
error bar that has been checked. (Fitting the *simulator's* constants to
those data is the second half, and it needs this uncertainty, so it comes
after.)

`cellsim plate --data reading.csv --map map.csv` reads a reader's grid
export or a tidy table, subtracts blanks, sets the plate's own controls
to 1, fits a four-parameter logistic in log concentration, and reports
GR metrics when the plate carries a time-zero read.

**Why GR as well as IC50.** A 72 h viability assay confounds potency with
division rate: a line that doubles twice in the assay has more to lose
than one that doubles once, so the same drug looks more potent in the
faster line. GR metrics divide by the control's own growth over the same
window (Hafner et al. 2016 Nat Methods 13:521) and are comparable across
lines and assay lengths. They need one extra read, at time zero.

### The confidence interval was wrong twice before it was right

The only way to know what a confidence interval is worth is to simulate
curves whose IC50 is **known** and count how often the interval contains
it (`scripts/validate_plate_fit.py`, 200 plates per row):

| Method | Coverage at a nominal 95 % |
|---|:-:|
| resample the raw residuals | **80 %** |
| scale residuals by √(n/(n−p)) | 88–92 % |
| bias-corrected and accelerated (BCa) | **88–90 %**, and no bias in the point estimate |

The first gap has a definite cause: residuals from a fitted curve are
smaller than the true errors, because the fit has already absorbed p of
the n degrees of freedom. The rest is what a small-sample bootstrap on a
nonlinear parameter does, and ~90 % is where ten concentrations land.
That number is printed next to every interval rather than left implied.

### What the plate design has to be, for the number to mean anything

This is the part worth giving a biologist before they run the assay.

**Concentrations per curve** (noise sd 0.05):

| Points | Coverage | Interval width |
|:-:|:-:|:-:|
| 5 | **68 %** | 1.85× |
| 7 | 86 % | 1.83× |
| 10 | 88 % | 1.75× |
| 12 | 84 % | 1.68× |

A five-point curve fits four parameters to five numbers. The interval it
returns looks no wider than a ten-point one and contains the truth two
times in three. `fit_4pl` therefore attaches an explicit note below seven
points, and `summary()` quotes 68 % rather than 90 % there.

**Noise** (ten points): coverage holds at 88–90 % up to a residual sd of
0.10, and falls to 84 % at 0.20 — where the interval has widened to 18×
and is no longer telling anyone anything.

**Range**, which turned out to matter more than either (50 plates per
row, since these fits are slow by construction):

| Tested range | Coverage | Interval width | Flagged extrapolated |
|---|:-:|:-:|:-:|
| 0.01–100 µM (brackets the IC50) | 86 % | 1.67× | 0/50 |
| 0.1–30 µM (brackets it) | 80 % | 1.76× | 0/50 |
| 0.001–0.5 µM (clips it from below) | **64 %** | **121×** | 20/50 |
| 10–1000 µM (clips it from above) | **44 %** | **168×** | 20/50 |

A range that sits entirely to one side gives a flat plate with no
midpoint to find, and is **refused** rather than fitted — the earlier
code fitted it anyway, taking minutes per call to extract a confident
IC50 from noise.

The in-between case is the one to understand. A range that merely
*clips* the response still fits, and the interval does not hide the
problem: it blows up to a hundredfold and more, which is the method
saying it does not know. What it does not do is contain the truth —
even that enormous interval misses it between a third and a half of the
time, because the data constrain the curve's shape only on one side and
the extrapolation can go anywhere. So a wide interval here is not a
conservative answer, it is an absent one. Those fits are flagged
`extrapolated`, and the flag, not the width, is what to act on.

### Does a calibrated engine predict what it was NOT calibrated on?

Fitting a constant so the simulated IC50 matches a measured one proves
nothing: the constant was chosen to make that number right. The claim
worth testing is the next one — that the calibrated model then says
something true about a condition it never saw — and
`scripts/validate_calibration.py` measures it with ground truth we
control. Perturb a drug's constant to a value the calibrator is never
told; simulate the noisy 72 h plate those cells would give; hand over
**only** the fitted IC50 and its interval; then ask about a 12 h pulse
washed out long before the readout, and check whether the predicted band
contains the truth.

**It failed the first time, and the failure was the useful part.**
Calibrated bands contained the held-out truth in **3 of 7 trials**. A
band advertised as covering the answer was wrong more often than right,
and nothing upstream was broken: the plate fit's coverage had been
measured, the calibration's propagation had been measured, and each
piece passed its own test. What had never been tested was the
composition.

The cause is that the band counted the measurement's uncertainty and
ignored the model's own.

### The engine's IC50 is a random variable, and more cells do not fix it

Repeating the identical simulation with a different seed moves the
engine's IC50 by 1.1–1.2×. That is not counting noise, and the
measurement that settles it is switching the cell-to-cell heterogeneity
off:

| `Params.het_sigma` | 24 cells/dose | 96 cells/dose |
|---|:-:|:-:|
| 0.30 (default) | **1.187×** | **1.108×** |
| 0 (identical cells) | 1.009× | 1.008× |

With identical cells the engine is deterministic to under a per cent.
With the default heterogeneity the spread is ten to twenty per cent and
**quadrupling the cells barely moves it**, because a population's IC50
is set by which resistant lineages happened to be drawn, not by how many
cells were counted — and a tail converges far more slowly than a mean.

That is a property of the model rather than a defect, and it is also the
right biology: a real colony's survival is dominated by its resistant
minority (Spencer et al. 2009 Nature 459:428, the same observation the
heterogeneity was built from). But it has three consequences that matter
to anyone using these numbers:

1. **A single-seed IC50 is good to about ±20 %**, so differences smaller
   than that between two runs are not differences.
2. **Any calibration from one measurement inherits that**, which is why
   the recovered constant scatters between 0.7× and 1.0× of the truth
   even with no measurement noise at all.
3. **More cells is the wrong remedy; more seeds is the right one.**

`predict_band` therefore spans the calibrated interval **and** seeds,
and the hold-out is scored against the truth's own median over the same
seeds, since no band can honestly be narrower than the quantity it is
predicting.

### Which published claims that ±20 % actually touches

A spread that large in every IC50 is alarming until it is checked
against the claims that rest on it, so here is that check rather than a
reassurance. Five seeds, A549:

| Claim | Across seeds | Verdict |
|---|---|:-:|
| Nutlin calls 11/11 lines correctly (Phase 3 headline) | **0 of 11 lines change call**; sensitive lines give viability 0.00, resistant 0.85–1.11 | safe |
| Paclitaxel cannot halve the colony below 24 h at any dose | half kill reached in **0/5 seeds at 3, 6 and 12 h; 5/5 at 24 h** | safe |
| Cisplatin follows C × T over 1–6 h (slope −1) | slope median **−0.974**, range −1.00 to −0.88; the −0.95 quoted above sits mid-range | safe |
| The C50 values themselves (e.g. 318 µM at 1 h) | spread **1.16–1.35×** between seeds | **quoted too precisely** |

The pattern is consistent and it is worth stating as a working rule for
this engine. **A comparison made within one seed is far more robust than
either number in it.** The seed fixes which lineages were drawn, and that
draw shifts a whole dose-response curve up or down together rather than
tilting it — so a slope, a ratio between two arms, or a threshold with a
wide margin survives, while the absolute concentration does not.

That is why the slope holds to ±0.06 while the C50s it is computed from
move by a third, and it retroactively justifies how the Phase-2 and
Phase-3 experiments were built: every one of them compares arms inside a
run rather than quoting a number from one.

What it does NOT justify is the single C50 figures quoted in the
schedule table above, which carry an unstated ±20 % and should be read
as one significant figure.

### What a calibrated prediction is actually worth

Twenty trials, with the band spanning the calibrated interval and seeds:

| | before the fix | after |
|---|:-:|:-:|
| band contains the held-out truth | 3/7 (43 %) | **16/20 (80 %)** |
| median band width | — | 0.435 surviving fraction |
| the truth's OWN spread across seeds | — | 0.173 |
| constant recovered | 0.79–0.83× | 0.84× median (0.69–0.98) |

Read the second and third rows together. The band is 0.44 wide and the
quantity it predicts moves by 0.17 on its own between runs, so roughly
40 % of the width is irreducible: no honest band here can be much
narrower. What that buys a user is an answer like *"between a quarter
and two thirds of the cells survive this pulse"* — enough to rule out
*nearly all die* and *nearly all live*, not enough to separate 35 % from
45 %. Two schedules whose bands overlap have not been shown to differ,
which is the whole reason `predict_band` returns a range.

**One thing tried and rejected.** If the engine's viability carries
±20 % run-to-run spread, averaging the calibration criterion over seeds
ought to help. It does not: three seeds returned byte-identical
constants in every case tried, at three times the cost, and moved the
offset below by under half a per cent. The reason is that the offset is
definitional, not noise, so averaging noise cannot reach it. The option
is kept in the signature and the default is one seed. What *does* move
with the seed is the answer itself — the same target calibrated under
different seeds spans about 1.14× — and that is what the prediction band
spans.

**A residual bias, measured and left in place.** The fitted constant
lands about 10 % below the truth, one-sidedly, and with the measurement
noise removed entirely it is still 0.88–0.92×. Two definitions of "the
IC50" are in play: calibration targets the concentration where median
viability is one half, while `engine.ic50` interpolates a crossing on a
log grid and on a curve this steep lands ~1.1× low. Targeting the
estimator instead would cost a full concentration series per iteration —
the thing that makes calibration fast enough to use — to chase a bias
that sits well inside the ±20 % the engine varies by anyway. It is
recorded rather than removed.

### Two bugs the plate files caught

Writing real plate files for the tests, rather than passing arrays
directly, found both:

* A grid export pads its unused wells, and `empty` was being classified
  as a **blank**. Forty-two zero-reading filler wells were averaged into
  the background, which silently rescaled every viability on the plate.
  `empty` is now its own role and is ignored.
* Readers write `B07` where plate maps write `B7`. Unnormalised, every
  well fails to match and the error surfaces far from its cause; well
  names are now normalised, and an unmapped well raises instead of being
  dropped.

Gated by `tests/cell/test_plate_smoke.py` (11 gates, 9 s), which builds
plate files on disk in both layouts.

## GDSC reference data (the Phase-1 validation target)

Extracted by `scripts/gdsc_reference.py` into
`benchmarks/cell/gdsc_reference.csv` from the public GDSC1 and GDSC2
fitted dose-response tables, release 8.4: **60 screens across the three
drugs and ten cell lines** (HCT116 and SW480 are in the library but have
no screens for these drugs). The raw files are not vendored; the
script's docstring has the two-line download.

| Drug | Source | Screens | In tested range | Quantitative target? |
|---|---|:-:|:-:|---|
| paclitaxel | GDSC2 | 10 | 10 | **Yes** |
| doxorubicin | GDSC1 | 20 | 19 | **Yes**, with the replicate spread below |
| cisplatin | GDSC1 + GDSC2 | 30 | 10 | **Partly.** Two-thirds of its screens never reached 50 % kill, so those lines give lower bounds only. |

Camptothecin is not a usable stand-in for doxorubicin: in GDSC2 four of
five lines have IC50 above its 0.1 µM top dose.

**Does TP53 status predict sensitivity? Only for one drug.** Geometric
mean IC50 over the lines with an in-range reference (October 2026, ten
lines: five TP53 wild-type, five mutant):

| Drug | p53 wild-type | p53 mutant | mutant / wild-type |
|---|:-:|:-:|:-:|
| doxorubicin | 0.059 µM (n=5) | 0.125 µM (n=5) | **2.1×** more resistant |
| paclitaxel | 0.039 µM (n=5) | 0.037 µM (n=5) | 0.94× — no effect |
| cisplatin | 5.74 µM (n=4) | 2.46 µM (n=2) | 0.43× — mutants *more* sensitive |

Only doxorubicin behaves as the textbook p53 story predicts. Paclitaxel
showing nothing is mechanistically right: it kills through mitotic
arrest, not DNA damage, so p53 should not matter. Cisplatin runs the
other way, on just two mutant lines, both of which happen to be
sensitive. So a model is right to make p53 matter for DNA-damage drugs
and not for tubulin drugs, and is wrong to treat p53 loss as blanket
resistance.

**Replicate spread sets the gate.** GDSC1 screened doxorubicin twice per
line. The two in-range IC50s of the same line differ by:

| Line | Spread |
|---|:-:|
| HeLa | 1.5× |
| A549 | 4.9× |
| HT-29 | 6.3× |
| MCF7 | 9.3× |

Median 5.6×. A model cannot be held to "within 3×" when the reference
data disagrees with itself by more than that, so the Phase-1 gate is the
replicate span (at least 3× either side of the geometric mean where
there is one screen) plus the correct most/least-sensitive ordering.

The p53 split is **not** clean in this data and must not be assumed:
cisplatin is less potent in p53-mutant HT-29 (122 µM) than in wild-type
A549 (20 µM) as expected, but p53-mutant MDA-MB-231 (43 µM) is *more*
sensitive than wild-type MCF7 (100 µM). The p53 effect is a measured
outcome to report, never a gate.

## Cell-level modules (`cellsim/cell`)

22 standalone tests, each asserting an analytic limit (steady state
reduces to each component, transient converges to steady state, mass
conserved, Monte-Carlo agrees with interval arithmetic). No calibration
against cell-line data yet; that is Phase 1 of `PLAN.md`.

## CI

`smoke.yml`: 54 always-on steps. The 51 that existed before the October
wiring were replayed locally on the merged tree: 51 pass, 0 fail, 14 min.
The five newly wired steps (29 test files) pass in seconds each.
