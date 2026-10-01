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

**Still to do, in order.** The measurement above says what is missing:
information that separates one line from another.

1. **Line-specific inputs from measured data.** CCLE/DepMap expression
   of ABCB1 (efflux), the BCL2 family (apoptotic priming) and repair
   genes, mapped onto the per-cell multipliers the engine already has.
   This is the one change that could move the paired test, because it
   is the only one that adds per-line information.
2. **Better doubling times**, or treat them as uncertain inputs. They
   disagree by up to 1.6× between sources and are currently one of only
   two line-specific numbers.
3. A shallower death response (curve shape above).

Re-run `cellsim validate-gdsc` after each; the gate is the paired sign
test, not the RMSE table.

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
