# Campaign 1 — Atomic → Molecular Foundation

**Scope horizon:** 1–2 years, 1 → 3–5 FTE.

**Single-sentence deliverable:** `cellsim-chem`, a Tier-4 chemistry
engine that screens 10⁴–10⁵ compounds for protein–ligand binding and
basic reactivity, emits calibrated-uncertainty predictions, and
passes PDBBind / ChEMBL / PoseBusters blind validation at r ≥ 0.7.

## Layers

Each layer pairs a numeric harness with a minimal real-time viewer
(see `src/render/README.md` and `src/viewer/README.md`).

### 1.1 Chem foundation — `src/chem/`
RDKit + OpenFF + AM1-BCC partial charges. Round-trip 10 k ChEMBL
compounds into OpenMM systems. Viewer: ligand ball-and-stick with
per-atom charge colouring.

### 1.2 Classical MD — `src/md/`
OpenMM driver (CUDA + Metal + CPU), ff14SB / ff19SB, TIP3P, PME,
equilibration protocols. Viewer: live protein ribbon + solvent
skin + RMSD gauge.

### 1.3 Docking + FEP — `src/chem/`, `src/md/`, `src/cache/`
GNINA pose prediction, `perses` relative FEP, pose/ΔG cache.
Viewer: receptor + pose overlay, ΔG bar with 95 % CI.

### 1.4 Quantum — `src/quantum/`
`xtb` GFN2 + PySCF DFT hook for reactive fragments. Reactive-
metabolite site prediction. Viewer: ESP isosurface + HOMO/LUMO.

### 1.5 ML potential — `src/md/` (extended)
MACE-OFF23 via `openmm-ml`. Deep-ensemble of 3–5 replicas for
force-level UQ. Viewer: side-by-side MACE vs classical + residual
heatmap.

### 1.6 Coarse-grained — `src/cg/`
Martini 3 bilayer + protein elastic network. Viewer: lipid heads
colour-coded + live area-per-lipid plot.

### 1.7 UQ scaffold — `src/uq/`
MAPIE conformal prediction around docking ΔG; deep-ensemble around
ML-potential forces. Viewer: calibration curve + reliability
diagram.

### 1.8 Blind-validation harness — `benchmarks/` + `src/viewer/`
`cellsim-chem --benchmark` drives PDBBind, CASF-2016, PoseBusters,
ChEMBL held-out panels. GitHub Actions CI with rented H100 for FEP;
refuses merge on pass-rate regression. Quarterly red-team slot.
Viewer: dashboard of per-benchmark tiles + r-scatter + red-team
leaderboard.

## Exit criteria (hard pass/fail)

1. PDBBind refined-set blind pose recovery ≥ 75 % within 2 Å RMSD.
2. ChEMBL held-out IC50 ranking Pearson r ≥ 0.7 on 5 kinase panels.
3. PoseBusters physical-validity pass rate ≥ 95 %.
4. Conformal UQ calibration error ≤ 10 % on held-out set.
5. Reactive-metabolite prediction matches literature on ≥ 15/20
   marketed CYP3A4 substrates.
6. Every prediction carries a calibrated uncertainty bar and
   method provenance (enforced by `src/uq/Prediction`).
7. Everything reproducible from a fresh clone +
   `docker compose up`.
8. Every layer's viewer renders correctly on the layer's
   reference scene (gate alongside numeric gates).

Failing any → extend Campaign 1. Do not advance to Campaign 2.

## Build-vs-buy

See the plan file
(`.claude-plans/make-this-plan-into-virtual-key.md`) for the full
tool table. Core rule: we integrate mature open-source stacks
(RDKit, OpenMM, OpenFF, GNINA, perses, xtb, PySCF, MACE, Martini 3,
MAPIE, Snakemake). We only write the glue, the physics-prior cache,
the blind harness, the UQ envelope, and the Campaign-1 → Campaign-2
bridge.

## Compute

Tier 1–3 (RDKit, OpenFF, docking, xTB, short MD) run locally on the
user's 64 GB Apple-silicon machine. Tier 4 (microsecond FEP, MACE
ensembles) rent NVIDIA H100 at ~$2/h; planned budget ~$10 k over
Campaign 1 for blind-validation runs.
