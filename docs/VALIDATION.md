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

Ten gates in `tests/cell/test_engine_smoke.py`, all passing, each
asserting a phenotype rather than "it ran":

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

Not yet calibrated: the per-drug potency gains in `cellsim/cell/library.py`
are placeholders, so absolute IC50s are not meaningful yet. Fitting them
against the reference table below is the next step.

## GDSC2 reference data (the Phase-1 validation target)

Extracted by `scripts/gdsc_reference.py` into
`benchmarks/cell/gdsc_reference.csv` from the public GDSC2 fitted
dose-response release 8.4. The raw file is not vendored; the script
prints the one-line download.

| Drug | Lines | Published IC50 range | Usable as a quantitative target? |
|---|:-:|---|---|
| paclitaxel | 5 | 0.011–0.089 µM, all within the tested range | **Yes** |
| cisplatin | 5 | 20–122 µM, all **above** the 8 µM top dose, AUC 0.92–0.99 | **No** — the screen never reached 50 % kill, so those IC50s are extrapolations. Gate on direction and rank only. |
| doxorubicin | 0 | absent from GDSC2 | Use GDSC1 explicitly, or substitute camptothecin (in GDSC2, same replication-coupled damage mechanism) |

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

**Known CI cost problem.** `tests/fep/test_binding_scaffold_smoke.py`
takes **35 minutes on a cold cache** on an Apple-silicon laptop (37 s
warm). PR #10 expanded it to 14 tests, several of which build full
solvated protein-ligand complexes (streptavidin, ubiquitin) with
AM1-BCC charges. It is scaffold-only (`sample=False`, no MD), so the
cost is parametrisation and solvation, not sampling. On a 2-core
GitHub runner with no warm cache this will dominate the job and may
exceed sensible limits. Fix before relying on CI: cache the built
systems as fixtures, or move the heavy builders to `workflow_dispatch`
and keep the cheap assertions always-on.

`smoke.yml`: 54 always-on steps. The 51 that existed before the October
wiring were replayed locally on the merged tree: 51 pass, 0 fail, 14 min.
The five newly wired steps (29 test files) pass in seconds each.
