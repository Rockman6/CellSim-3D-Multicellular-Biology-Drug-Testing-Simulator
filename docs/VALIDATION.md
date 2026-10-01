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
  in NVT; nothing in `src/fep` applies a barostat or equilibrates the
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

## Cell-level modules (`src/cell`)

22 standalone tests, each asserting an analytic limit (steady state
reduces to each component, transient converges to steady state, mass
conserved, Monte-Carlo agrees with interval arithmetic). No calibration
against cell-line data yet; that is Phase 1 of `PLAN.md`.

## CI

`smoke.yml`: 54 always-on steps. The 51 that existed before the October
wiring were replayed locally on the merged tree: 51 pass, 0 fail, 14 min.
The five newly wired steps (29 test files) pass in seconds each.
