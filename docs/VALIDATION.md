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
MCF7 only, so eight lines are held out.

| Drug | Fitted constant | Engine in span | Null in span | log10 RMSE engine (held-out) | log10 RMSE null |
|---|---|:-:|:-:|:-:|:-:|
| cisplatin | `k_damage_per_uM_h` = 0.000954 | 9/10 | 7/10 | 0.54 (n=5) | 0.47 |
| doxorubicin | `k_damage_per_uM_h` = 0.00794 | 8/10 | 8/10 | 0.45 (n=8) | 0.46 |
| paclitaxel | `partition` = 0.153 | 10/10 | 8/10 | 0.31 (n=8) | 0.35 |
| **held-out** | | **22/25** | 18/25 | | |

### Exit gate: NOT met

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
scale from one fitted constant: 22 of 25 held-out predictions land
inside GDSC's replicate span against the null's 18, and A549 cisplatin
predicts 9.0 µM against 9.8 measured. Paclitaxel is the best case,
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
