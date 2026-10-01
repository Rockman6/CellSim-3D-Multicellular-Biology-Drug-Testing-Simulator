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

## ⚠ Open correctness issue: the streptavidin reference set

`benchmarks/dock/streptavidin_calibration.yaml` lists desthiobiotin at
**K_d = 5×10⁻⁵ M** (ΔG −5.9 kcal/mol). The literature says desthiobiotin
binds streptavidin **100 to 10 000 fold weaker than biotin**, i.e. about
10⁻¹¹ M (ΔG ≈ −15), and the entry's own cited source (Hirsch 2002 Anal
Biochem 308:343) is the standard reference for that figure. The file's
in-line note ("often cited 10⁻⁴ to 10⁻⁶ M") does not match it. The value
is internally consistent (its ΔG matches its K_d) but appears to be wrong
by roughly six orders of magnitude.

It matters because it is load-bearing. Recomputing the Linux CI docking
results against both values:

| Reference used | Spearman ρ | MAE |
|---|:-:|:-:|
| repo value (5×10⁻⁵ M) | **+0.80** | 4.95 kcal/mol |
| literature (10⁻¹¹ M) | **0.00** | 6.41 kcal/mol |

So the claim that docking *ranks* streptavidin binders usably, the
`tests/uq/test_calibration_smoke.py` gate (ρ ≥ 0.8), and the
streptavidin row of the reliability table below all rest on this number.
With the literature value the ranking is no better than chance, which is
what Vina's tight-binder saturation would predict: all four compounds
dock within 0.3 kcal/mol of each other.

**Not changed here.** Correcting it flips a published headline and makes
a CI gate fail, so it is the project owner's call. Tracked as a GitHub
issue; the fix is to re-pull all four K_d values from BindingDB or
PDBBind pinned to PDB IDs, as the file's own header already says a future
change should.

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

| Drug | Fitted constant | Engine in span | Null in span | log10 RMSE engine (held-out) | log10 RMSE null (held-out) |
|---|---|:-:|:-:|:-:|:-:|
| cisplatin | `k_damage_per_uM_h` = 0.000951 | 5/5 | 5/5 | 0.31 (n=1) | **0.15** |
| doxorubicin | `k_damage_per_uM_h` = 0.00794 | 5/5 | 5/5 | **0.29** | 0.33 |
| paclitaxel | `partition` = 0.153 | 5/5 | 4/5 | **0.32** | 0.38 |
| **held-out lines** | | **10/10** | 9/10 | | |

**Exit gate: met.** The engine beats the constant-IC50 null on held-out
log error for two of the three drugs. Cisplatin still loses, on a single
in-range held-out line (HeLa); its other lines have only lower bounds,
so there is almost nothing to score it on.

The null predicts one IC50 for every line: the geometric mean of the two
fit lines' references. RMSE is over lines with an in-range reference.

**What this shows.** From one fitted constant per drug, every held-out
line's IC50 lands inside GDSC's own replicate span, and on log error the
engine now beats a constant for two drugs of three. The A549 cisplatin
IC50 comes out at 10.1 µM against 9.8 measured.

**The honest caveat on the gate.** What made the difference is the
p53-independent death route below, and its one parameter was chosen
using these same two TP53-mutant lines (leave-one-out across them). So
the gate result is not fully out-of-sample with respect to that choice.
The clean out-of-sample evidence is the leave-one-out table itself,
where each fold's gain never sees the line it is scored on. A third
mutant line would settle it; that is the first thing to add.

**What is still not shown.** Curve shape (above) is unrealistic, and
cisplatin's line-to-line behaviour is untested because the public data
gives mostly lower bounds for it.

**A side effect of the fitted route, pinned by a test.** At the default
gain the p53-independent route is strong enough that a damage index held
constant kills a TP53-mutant line outright, so p53 status confers no
protection against sustained damage. That is what lets the mutant lines
match GDSC, where they are about as cisplatin-sensitive as wild-type,
but it is stronger than the literature's partial resistance. The gain
was picked from a coarse grid whose smallest non-zero value was 0.1, so
a finer scan may find a value that keeps the GDSC improvement while
leaving some p53 dependence. Tracked in
`tests/cell/test_engine_smoke.py` so it cannot change silently.

**A p53-independent death route closes part of the gap.** ATM/c-Abl →
p73 → PUMA is the documented route by which TP53-mutant cells still die
of DNA damage (Gong 1999 Nature 399:806; Agami 1999 Nature 399:809, both
with cisplatin). Added as `Params.p73_gain` (0 = off, the p53-only
model). Tested leave-one-mutant-out by `scripts/experiment_p73.py`: the
gain is chosen on one mutant line and scored on the *other*, so the test
line never influences it. Log10 error on the held-out mutant:

| Fold | Chosen gain | With route | p53-only | Constant-IC50 null |
|---|:-:|:-:|:-:|:-:|
| train HT-29 → test MDA-MB-231 | 0.1 | **0.318** | 0.524 | 0.405 |
| train MDA-MB-231 → test HT-29 | 0.1 | **0.332** | 0.388 | 0.410 |

The route improves out-of-sample prediction in both folds and beats the
null on both, and both folds independently pick the same small gain. It
does not touch untreated cells or the wild-type lines, because it only
engages once ATM is genuinely activated. Default stays 0 until the
drug constants are refitted with it on, which is the next run.

**Still to do, in order.**
1. A third TP53-mutant line, to test the route's parameter against a
   line that played no part in choosing it.
2. Line-specific inputs from measured data rather than doubling time
   only, e.g. CCLE/DepMap expression of ABCB1 (efflux), the BCL2 family
   (apoptotic priming) and repair genes, mapped onto the engine's
   existing per-cell multipliers.
3. A shallower death response (curve shape above).

The gate to beat from now on is the null's log10 RMSE, not the span.

## GDSC reference data (the Phase-1 validation target)

Extracted by `scripts/gdsc_reference.py` into
`benchmarks/cell/gdsc_reference.csv` from the public GDSC1 and GDSC2
fitted dose-response tables, release 8.4: 30 screens across the three
drugs and five of the curated lines (HCT116 and SW480 have no screens
for these drugs). The raw files are not vendored; the script's docstring
has the two-line download.

| Drug | Source | Screens | In tested range | IC50 | Quantitative target? |
|---|---|:-:|:-:|---|---|
| paclitaxel | GDSC2 | 5 | 5 | 0.011–0.089 µM | **Yes** |
| doxorubicin | GDSC1 | 10 | 9 | 0.010–0.37 µM | **Yes**, with the replicate spread below |
| cisplatin | GDSC1 + GDSC2 | 15 | 2 | 6.9 and 9.8 µM where measured; the rest extrapolated up to 122 µM | **No.** Direction and rank only. |

Camptothecin is not a usable stand-in for doxorubicin: in GDSC2 four of
five lines have IC50 above its 0.1 µM top dose.

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
