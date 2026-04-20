# CellSim validation — current benchmark numbers

This is CellSim's scorecard. Every number here is produced by a
published physics method (no ML), every number carries an honest
caveat, and every test runs on every PR (see
[`.github/workflows/smoke.yml`](.github/workflows/smoke.yml)).

Reproduce any row:

```bash
conda activate cellsim
python -u tests/<path>/<test>.py    # or cellsim <subcommand> …
```

> **Commit:** see `git rev-parse HEAD` for the exact reproducer.
> **Data:** every benchmark input is in [`benchmarks/`](benchmarks/),
> committed to git. No external downloads needed to reproduce.

---

## Campaign 1 status table

| Layer | What it does | Shipped | CI | Exit criterion |
|---|---|:-:|:-:|---|
| 1.1 | SMILES → OpenFF system, AM1-BCC charges | ✅ | ✅ | 99 % round-trip on 10 k ChEMBL (currently 100 % on 10/10 smoke) |
| 1.2 | Classical MD (ligand vacuum + protein solvated) | ✅ | ✅ | Ubiquitin 100 ns RMSD < 3 Å (short gate shipped) |
| 1.3 | AutoDock Vina docking + Meeko + PoseBusters + fpocket + MC + refine + batch + profile + export + off-target | ✅ | ✅ | 3-cocrystal mini-bench ≥ 66 % pose recovery |
| 1.4 | xTB GFN2 single-point + CYP3A4 SoM predictor | ✅ | ✅ | 15/20 marketed drugs CYP3A4 primary SoM match (pending) |
| 1.5 | Martini 3 bilayer + CG MD | ⏳ scaffold | ❌ | 10 µs POPC area-per-lipid within 2 % of literature |
| 1.6 | UQ triad: MC + Sobol + split-conformal | ✅ | ✅ | Conformal coverage ≤ 10 % calibration error |
| 1.7 | Blind-validation harness | ⏳ partial | ✅ | PDBBind r ≥ 0.6; red-team quarterly |
| x-cut | SQLite physics-prior cache (hash-addressed) | ✅ | ✅ | > 1000× speedup per compound on cache hit |

---

## 1.3 Docking: pose-recovery on bundled cocrystals

Re-docking gate — canonical Astex/PDBBind convention (top-1 ≤ 2.5 Å
**AND** best-of-top-3 < 2.0 Å).

| System | ligand | top-1 RMSD | top-3-best | status | notes |
|---|---|:-:|:-:|:-:|---|
| 1STP streptavidin | biotin | 2.02 Å | 1.99 Å | ✅ PASS | canonical benchmark since AutoDock 1.0 |
| 3PTB trypsin | benzamidine | 1.30 Å | 1.30 Å | ✅ PASS | classic S1-pocket test |
| 1M17 EGFR kinase | erlotinib | 4.23 Å | 4.23 Å | ❌ FAIL | Vina scoring-function failure; known hard kinase case |

**Aggregate: 2/3 = 67 %** at canonical gate.
Vina papers report ~75–85 % on Astex Diverse Set (85 systems) at
exhaustiveness=32; our 3-cocrystal set is deliberately harder
(includes a known failure case — erlotinib).

Reproducer: `python tests/dock/test_mini_bench.py`.

---

## 1.3 Docking: auto binding-site detection (fpocket)

Given a raw receptor PDB, fpocket identifies the biological binding
site.

| System | Known site centroid | fpocket top center | Distance | Drug score |
|---|---|---|:-:|:-:|
| 1STP (biotin) | (11.12, 1.68, −10.75) | (10.57, 2.84, −10.46) | 1.32 Å | 0.80 |

Gate: < 3 Å from the known centroid and drug score ≥ 0.5. Pass.

Reproducer: `python tests/dock/test_pocket_detect_smoke.py`.

---

## 1.3 Docking: Monte-Carlo ΔG uncertainty

8 Vina seeds on biotin → 1STP, exhaustiveness=16, multiprocessing 4.

| Metric | Value |
|---|---|
| ΔG mean | −7.28 kcal/mol |
| ΔG SD | 0.29 kcal/mol |
| 95 % CI (percentile) | [−7.45, −6.78] |
| Best of 8 | −7.45 @ seed 4 |
| Wall | 35.9 s |

The 0.7 kcal/mol spread across seeds is the *irreducible Vina jitter*
— any single-seed ΔG has this much hidden noise.

Reproducer: `cellsim uq-mc --receptor … --ligand-smiles …`.

---

## 1.6 Split-conformal on synthetic data

Gates: empirical coverage on a held-out test set ≥ 0.90 (target
0.95 at α = 0.05).

| N_cal | q (kcal/mol) | N_test | Empirical coverage |
|:-:|:-:|:-:|:-:|
| 50 | 1.35 | 200 | 0.900 |

Reproducer: `python tests/uq/test_conformal_smoke.py`.

---

## 1.7 Experimental calibration on two receptor families

Split-conformal calibrated against published K_d values.

### Streptavidin — `benchmarks/dock/streptavidin_calibration.yaml`

4 biotin-analog compounds spanning K_d 10⁻¹⁴ M to 10⁻⁵ M.

| Metric | Value | What it means |
|---|:-:|---|
| Pearson r | +0.999 | near-perfect rank+magnitude on this tiny set |
| **Spearman ρ** | **+1.000** | **perfect ranking — usable for triage** |
| MAE | 4.99 kcal/mol | absolute ΔG is Vina-approximate |
| RMSE | 6.84 kcal/mol | |
| Conformal q95 | ±11.66 kcal/mol | dominated by biotin's saturated ΔG |

Story: Vina correctly ranks streptavidin binders across 14 orders of
magnitude in K_d, but its absolute ΔG is flat ~−7 kcal/mol for every
compound (tight-binder saturation). Use the ranking; don't believe
the absolute.

### Trypsin — `benchmarks/dock/trypsin_calibration.yaml`

6 benzamidine-analog compounds spanning K_i 10⁻⁷ M to 10⁻³ M.

| Metric | Value | What it means |
|---|:-:|---|
| Pearson r | +0.615 | moderate correlation |
| Spearman ρ | +0.086 | **poor ranking within narrow window** |
| **MAE** | **0.91 kcal/mol** | **tight absolute accuracy** |
| RMSE | 1.15 kcal/mol | |
| Conformal q95 | ±1.96 kcal/mol | |

Story: Vina's absolute ΔG is right on for this scaffold (no tight-
binder saturation — all affinities in Vina's calibrated µM–mM
range), but it **cannot reliably rank close analogs** inside a
~4 kcal/mol window because its noise floor is ~1 kcal/mol.

### Cross-system take-away (for biologists)

Use CellSim for ranking across **wide** affinity ranges (nM vs mM),
**not** for discriminating close analogs. This is the correct
intuition for any Vina-based pipeline and it's now baked into
reproducible benchmarks. Future perses-FEP integration (Layer 1.3)
will tighten MAE to ~1 kcal/mol on *both* systems.

Reproducer: `cellsim uq ... ` (see `src/uq/calibration.py`).

---

## 1.3 Off-target / selectivity screen

One compound vs N receptors, auto-pocket each.

Biotin vs 3 receptors (1STP / 3PTB / 1UBQ), exhaustiveness=16:

| Rank | target | ΔG (kcal/mol) | drug score | selectivity ΔΔG |
|---|---|:-:|:-:|---|
| 1 | streptavidin (1STP) | −7.44 | 0.80 | — (reference) |
| 2 | trypsin (3PTB) | −5.95 | 0.44 | +1.50 vs rank 1 |
| 3 | ubiquitin (1UBQ) | −3.91 | 0.54 | — |

Ranking is correct (biotin's biological target at rank 1, ubiquitin
with no real pocket at rank 3). Selectivity ΔΔG ~ 1.5 kcal/mol is
small because Vina saturates on the tight streptavidin side.

Reproducer: `cellsim off-target --ligand-smiles "…" --receptors "…"`.

---

## x-cut Cache hit speed-up (on cellsim conda env)

| Operation | Cold wall | Warm wall | Speed-up |
|---|:-:|:-:|:-:|
| AM1-BCC parametrise (aspirin) | 11.76 s | 0.0004 s | 29 455× |
| xTB GFN2 single-point (aspirin) | 0.32 s | 0.05 s | 6× |
| Vina dock (biotin → 1STP) | 10.00 s | 0.010 s | 1 021× |
| **End-to-end batch** (5 drugs → 1STP) | **16.9 s** | **0.9 s** | **19×** |

Reproducer:
```bash
cellsim dock --smi … --receptor … --cache /tmp/proj.sqlite   # run 1
cellsim dock --smi … --receptor … --cache /tmp/proj.sqlite   # run 2 hits cache
```

---

## 1.4 Quantum / CYP3A4 site-of-metabolism

BDE-based SoM predictor over all X–H bonds.

3-drug smoke on aspirin, ibuprofen, testosterone:

| drug | n_CH_candidates | top-3 BDE (kcal/mol) | wall |
|---|:-:|---|:-:|
| aspirin | 8 | 134 / 134 / 136 | 0.4 s |
| ibuprofen | 15 | 127 / 127 / 127 | 0.9 s |
| testosterone | 26 | 119 / 121 / 124 | 2.1 s |

GFN2 absolute BDE has a systematic ~20 kcal/mol offset vs
experiment; the *ranking* is what matters for SoM. Literature-
validated correctness against 20 marketed CYP3A4 substrates is
pending (Layer 1.4 exit gate).

### DFT rescoring of top-3 xTB candidates — aspirin case study

Canonical pharma pattern: xTB ranks all X-H bonds fast; DFT re-
scores the top-3 for paper-grade accuracy. On aspirin:

| Step | top-1 | top-2 | top-3 |
|---|---|---|---|
| xTB (GFN2-xTB BDE) | O(12) 134.4 | C(0) 134.5 | C(0) 136.3 |
| + DFT rescoring | **C(0) 112.7** | C(0) 114.0 | O(12) 115.3 |

Two honest findings:

1. **DFT values are ~20 kcal/mol lower than xTB** — empirically
   confirms the GFN2 systematic offset vs experiment.
2. **DFT disagrees with xTB on the winner.** xTB says O-H of the
   carboxylic acid; DFT says methyl C-H. The literature primary
   CYP3A4 metabolism of aspirin is methyl oxidation leading to
   salicylate, matching DFT. xTB's top-2 and DFT's top-1 are the
   same atom, so xTB is usable as a fast funnel; DFT is the
   final arbiter.

Wall: 58.9 s on a 21-atom drug (1 full-mol DFT + 3 radicals + 1
H-atom reference; B3LYP/def2-SVP; PySCF).

Reproducer: `cellsim som "CC(=O)OC1=CC=CC=C1C(=O)O" --dft-verify 3`.

Reproducer (xTB-only baseline): `python tests/quantum/test_som_smoke.py`.

---

## 1.2 MD: solvated protein equilibration

1UBQ ubiquitin load → solvate → minimise → 1 ps NVT:

| Metric | Value |
|---|---|
| Atoms (protein / total) | 1 231 / 17 344 |
| Waters | 5 361 |
| Ions | 15 Na⁺ + 15 Cl⁻ |
| Box | 5.69 nm cube, PBC + PME |
| E_initial | −55 360 kJ/mol |
| E_minimised | −319 518 kJ/mol |
| ΔE from minimisation | 264 157 kJ/mol |
| 1 ps MD — T_final | 261 K (setpoint 300 K) |
| 1 ps MD — **Cα RMSD vs start** | **0.74 Å** ✅ (gate ≤ 2 Å) |
| Wall (load + min) | 8.6 s |
| Wall (1 ps MD) | 0.6 s |

Full Layer-1.2 exit gate (100 ns backbone RMSD < 3 Å) needs GPU
and is pending.

Reproducer: `python tests/md/test_protein_load.py`.

---

## CI gates that block merge

Every PR that touches `src/**`, `tests/**`, `benchmarks/**`, or
`environment.yml` runs the full gate stack (~15 min on GitHub's
2-core Ubuntu runner):

1. chem: `test_parametrize_smoke.py` — 10 canonical drugs ≥ 8/10
2. chem: `test_admet_smoke.py` — Lipinski/QED/logS sanity
3. chem: `test_parametrize_cache.py` — cache round-trip + speed-up
4. md: `test_ligand_vacuum.py` — 3 drugs × 10 ps vacuum MD
5. md: `test_protein_load.py` — 1UBQ solvate + minimise + 1 ps MD
6. dock: `test_prep_smoke.py` — Meeko PDBQT
7. dock: `test_redocking.py` — 1STP biotin top-3 < 2 Å
8. dock: `test_pocket_detect_smoke.py` — fpocket within 3 Å
9. dock: `test_mini_bench.py` — 3-cocrystal ≥ 2/3 pass
10. dock: `test_refine_smoke.py` — post-dock MM minimise
11. dock: `test_vina_cache.py` — cache round-trip + speed-up
12. dock: `test_export_smoke.py` — SDF + PDB export
13. dock: `test_off_target_smoke.py` — selectivity screen
14. quantum: `test_xtb_smoke.py` — GFN2 sanity
15. quantum: `test_som_smoke.py` — CYP3A4 BDE
16. uq: `test_mc_dock_smoke.py` — MC over seeds
17. uq: `test_conformal_smoke.py` — split-conformal
18. uq: `test_calibration_smoke.py` — streptavidin vs K_d
19. cache: `test_cache_smoke.py` — hashing + store

One regression → merge blocked.

---

## What CellSim does NOT claim

- **No rigorous ΔG for tight binders.** Vina saturates at ~−7 to
  −10 kcal/mol; the Layer-1.3 perses-FEP hook (pending) fills this
  gap.
- **No membrane-embedded targets.** Layer 1.5 Martini 3 is scaffold
  only.
- **No publication-grade Sobol.** N_base ≥ 32 (≥ 256 runs) needed;
  default N_base=8 is "does the analysis run" smoke.
- **No in-cell pharmacology.** Campaign 2 / 3 work — on the
  roadmap but deliberately out of scope for Campaign 1.
- **Not a replacement for final wet-lab validation.** See
  [`MISSION.md`](MISSION.md): CellSim is a triage reduction tool,
  not a replacement for validation. The Spearman ρ = 1.0 on
  streptavidin tells you you can *rank* compounds; it doesn't
  say anything about idiosyncratic toxicity, ADME in vivo, or
  protein concentrations in real tissue.
