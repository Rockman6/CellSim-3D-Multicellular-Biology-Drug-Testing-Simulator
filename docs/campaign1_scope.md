# Campaign 1 — Atomic → Molecular Foundation (non-AI)

**Amendment:** 2026-04-20. Original Layer 1.5 (MACE-OFF23 ML
potential) removed per the non-AI / physics-only commitment in
`MISSION.md`. Campaign 1 now has 7 layers.

**Scope horizon:** 1–2 years, 1 → 3–5 FTE.

**Current progress (April 2026):** Layers 1.1 / 1.2 / 1.4 / 1.6
shipped MVP, Layer 1.3 docking-side shipped with triage / strain /
shortlist / off-target / DDI-risk surfaces (FEP integration open),
Layer 1.5 scaffolded, Layer 1.7 harness partial (PoseBusters
integrated). See `ROADMAP.md` for the current status table and
`BENCHMARKS.md` for the numbers being enforced by the 34-gate CI.

**Single-sentence deliverable:** `cellsim-chem`, a transparent,
physics-first chemistry engine that screens 10⁴–10⁵ compounds for
protein–ligand binding and basic reactivity, emits uncertainty-
quantified predictions via parameter-sweep + Sobol + Monte-Carlo
methods, and passes PDBBind / ChEMBL / PoseBusters blind validation
at Pearson r ≥ 0.7.

## Layers

Each layer pairs a numeric harness with a minimal real-time viewer
(see `src/render/README.md` and `src/viewer/README.md`). Every
method is auditable end-to-end. No learned surrogates.

### 1.1 Chem foundation — `src/chem/`
RDKit + OpenFF Sage 2.1.0 + AM1-BCC partial charges (via
AmberTools `antechamber` / `sqm`). Round-trip 10 k ChEMBL compounds
into OpenMM systems. Every result carries method + force-field +
tool versions. Viewer: ligand ball-and-stick with per-atom charge
colouring.

### 1.2 Classical MD — `src/md/`
OpenMM driver (CUDA + OpenCL + CPU), AMBER ff14SB / ff19SB for
proteins, TIP3P or TIP3P-FB water, PME electrostatics, Langevin
NVT at seeded temperature, HBonds constraints, rigid water,
equilibration protocols. GROMACS as optional independent-engine
cross-check. Viewer: live protein Cα trace + T / E / Cα-RMSD panel.

### 1.3 Docking + FEP — `src/dock/`, `src/cache/`
**AutoDock Vina as primary docking engine** (empirical, auditable
scoring function). `perses` + `openmmtools` alchemical FEP for
relative binding free energies. SQLite + HDF5 pose / ΔG cache keyed
by `(ligand_hash, receptor_hash, method, ff_version)`. Viewer:
receptor Cα ribbon + top-pose overlay + ΔG bar with Monte-Carlo CI.

**No CNN-scored mode.** Per the 2026-04-20 non-AI amendment in
`MISSION.md`, GNINA CNN scoring is excluded from Campaign 1.
AutoDock Vina's empirical function is the only scoring path;
pose trust is enforced by the UFF-ensemble strain gate
(`src/dock/strain.py`) and PoseBusters geometry checks, both of
which are physics/rule-based.

### 1.4 Quantum — `src/quantum/`
`xtb` GFN2-xTB semi-empirical for reactive fragments, geometry
optimisation, HOMO/LUMO, ESP. PySCF DFT hook for cases where xTB
is insufficient. Reactive-metabolite site prediction against 20
marketed CYP3A4 substrates. Viewer: ESP isosurface + HOMO/LUMO
density clouds on the ligand.

### 1.5 Coarse-grained — `src/cg/`
Martini 3 lipid bilayer + protein elastic network via
`martinize2 / vermouth`. 10 µs of a 200 nm² POPC bilayer; area-per-
lipid within 2 % of literature. Viewer: lipid heads colour-coded +
live area-per-lipid time-series.

### 1.6 UQ scaffold — `src/uq/`
Non-AI uncertainty quantification:

- **Primary: mechanistic / statistical sampling.** Parameter sweeps
  (grid + Latin hypercube) over rate constants and force-field
  parameters; Monte-Carlo ensembles over charges; Sobol global
  sensitivity indices via `SALib` to identify which parameters
  actually move which predictions.
- **Secondary: post-hoc distribution-free statistical bounds.**
  MAPIE conformal prediction is acceptable only as a non-parametric
  wrapper that maps a held-out calibration set to distribution-free
  CI bounds on any physics prediction. It provides *statistical*
  bounds, not *mechanistic* insight, and is documented as such on
  every predictor it wraps.

No deep ensembles of learned potentials. No neural UQ. Viewer:
sensitivity tornado plot + Sobol index bars + calibration curve.

### 1.7 Blind-validation harness — `benchmarks/` + `src/viewer/`
`cellsim-chem --benchmark` drives PDBBind, CASF-2016, PoseBusters,
ChEMBL held-out panels. GitHub Actions CI runs smoke gates on every
PR; longer FEP gates run on rented GPU via `workflow_dispatch`.
Quarterly red-team slot where external chemists submit compounds
designed to break the model; every failure becomes a new regression
benchmark. Viewer: per-benchmark pass/fail tiles + r-scatter +
red-team leaderboard.

## Exit criteria (hard pass/fail)

1. PDBBind refined-set blind pose recovery ≥ 75 % within 2 Å RMSD,
   using AutoDock Vina (no CNN scoring; see MISSION.md).
2. ChEMBL held-out IC50 ranking Pearson r ≥ 0.7 on 5 kinase panels.
3. PoseBusters physical-validity pass rate ≥ 95 %.
4. UQ calibration error ≤ 10 % on held-out set via MAPIE conformal
   wrapper; Sobol sensitivity indices documented for every headline
   prediction.
5. Reactive-metabolite prediction matches literature on ≥ 15/20
   marketed CYP3A4 substrates.
6. Every prediction carries a calibrated uncertainty bar + method
   provenance (enforced by `src/uq/Prediction`). Every rate constant
   cites a PMID or a cached physics-calculation ID.
7. Everything reproducible from a fresh clone +
   `docker compose up`.
8. Every layer's viewer renders correctly on the layer's reference
   scene (gate alongside numeric gates).

Failing any → extend Campaign 1. Do not advance to Campaign 2.

## Build-vs-buy (non-AI)

Open-source integrations (BUY):
RDKit, OpenFF Toolkit / Sage, OpenMM, AMBER ff14SB/ff19SB, CHARMM36m,
TIP3P-FB, PDBFixer, PROPKA, DimorphiteDL, AmberTools (antechamber,
sqm), AutoDock Vina, perses + openmmtools + openff-evaluator,
Martini 3 + martinize2 / vermouth, xtb, PySCF, SALib, MAPIE,
Snakemake, Docker, HDF5, SQLite, Parquet, PDBBind, CASF-2016,
ChEMBL, BindingDB, PoseBusters.

Glue we write (BUILD): pipeline orchestrator, physics-prior cache,
deterministic blind-benchmark harness, UQ envelope + sensitivity
infra, Campaign-1 → Campaign-2 rate-law bridge.

Explicitly excluded: MACE / MACE-OFF23 / NequIP / Allegro / OrbNet /
any ML potential as a drop-in force path; AlphaFold retraining;
end-to-end learned docking; neural scoring as sole evidence; deep
ensembles for UQ.

## Compute

Tier 1–3 (RDKit, OpenFF, docking, xTB, short MD) run locally on a
64 GB Apple-silicon laptop. Tier 4 (µs FEP, long ubiquitin MD,
10 k ChEMBL full gate) run on rented NVIDIA H100 at ~$2/h; planned
budget ~$10 k over Campaign 1 for blind-validation and scale runs.
