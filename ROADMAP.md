# CellSim Roadmap

Authoritative strategic plan: **`GOAL`** (5-campaign hierarchical
multi-scale hybrid simulator, 4–7 years).

Approved execution plan for the current restart:
**`.claude-plans/make-this-plan-into-virtual-key.md`** (migration
manifest + Campaign-1 scope + build-vs-buy table).

## Where we are

**Campaign 1, mid-execution, April 2026.**

The pre-restart HeLa / p53 / cisplatin Campaign-2 prototype is
frozen under `OLD/` as a regression snapshot. The new top-level
`src/` tree is populated across Layers 1.1–1.6 with the exceptions
noted below. 34 smoke gates run on every PR (see `README.md`
§"Validation that runs on every PR").

## Campaign 1 — Atomic → Molecular Foundation (Years 1–2, non-AI)

Per the 2026-04-20 non-AI amendment (see `MISSION.md`), Campaign 1
is 7 layers, no ML potentials. Original Layer 1.5 (MACE-OFF23) was
removed; the coarse-grained, UQ, and blind-validation layers
renumbered accordingly.

- **1.1 Chem foundation (RDKit + OpenFF + AM1-BCC)** — DONE.
  Parametrise, ADMET (Lipinski + QED + logS + BBB + hERG + Ames),
  AM1-BCC cache, profile-PNG dashboard.
- **1.2 Classical MD (OpenMM + ff14SB + TIP3P)** — DONE (MVP).
  Vacuum Langevin on ligands, 1UBQ solvate + 1 ps MD. Full
  100-ns ubiquitin gate pending GPU runner.
- **1.3 Docking + FEP (AutoDock Vina primary; non-AI alchemical
  FEP via openmmtools + pymbar)** — Docking DONE; alchemical FEP
  pipeline COMPLETE end-to-end. Docking side covered: Meeko prep,
  re-dock gate, mini-bench, fpocket auto-detect, Vina cache,
  MC-dock, refine, SDF/PDB export, off-target panel, CYP3A4
  inhibition screen, UFF-ensemble strain as a cross-cutting
  pose-trust gate, triage verdict column, strain-gate top-pose
  promotion, shortlist CSV, triage-PNG dashboard, kinase-receptor
  heads-up, "next steps" paste-ready guidance. FEP side covered:
    * Milestone A (hydration): `compute_hydration_dg` ready;
      FreeSolv-12 gate running on M5 Max, awaiting tarball.
    * Milestone B (binding): DDM + hybrid amber14/SMIRNOFF
      builder; `compute_absolute_binding_dg` +
      `compute_relative_binding_ddg`; end-to-end sampled ΔG_bind
      validated on 1ubq+methane (82.5 s CPU).
    * Biologist entry points: `cellsim fep-binding {dg,ddg,bench,
      validate}` + `cellsim fep-report` + `cellsim bench-all`.
      Validate emits Lipinski table + wall-time estimate with
      local-GPU detection; bench has --resume + Ctrl-C handler +
      per-compound ETA + crash-proof incremental CSV; report
      autodiscovers latest run, infers YAML kind, kind-aware
      sign-critical rule + parity PNG labels.
    * 50+ smoke tests CI-wired; `cellsim doctor` step 13 drives
      every benchmarks/fep/*.yaml through `validate` on every PR.
    * Open: sampled Phase-2 runs (streptavidin + EGFR) pending
      Milestone A clearance + GPU time.
- **1.4 Quantum (xtb + PySCF)** — DONE (MVP). GFN2 single-point,
  homolytic C-H BDE ranking for CYP3A4 SoM, optional DFT top-N
  rescore, CYP3A4 heme-accessibility pose-SoM with ensemble pose
  search. 2/3 on the literature validation set (aspirin ✓,
  midazolam ✓, diazepam ✗ — documented as xTB-BDE limitation
  on α-N-methyl sites, not a pose-sampling issue).
- **1.5 Coarse-grained (Martini 3)** — SCAFFOLD only. Stubs raise
  NotImplementedError; martinize2 present in env via `vermouth`.
- **1.6 UQ (Sobol + Monte-Carlo + parameter sweeps; MAPIE
  conformal)** — DONE. MC-dock over seeds, Sobol sensitivity,
  split-conformal quantiles, streptavidin + trypsin + EGFR
  calibration bundles (EGFR is an honest rank-order-failure
  finding on kinases, pinned as a regression witness).
- **1.7 Blind-validation harness (PDBBind / CASF / PoseBusters /
  ChEMBL)** — PARTIAL. PoseBusters integrated; PDBBind / CASF
  bundles not yet populated.

Every layer ships **numeric harness + minimal viewer** in lockstep.
See `docs/campaign1_scope.md` for exit criteria, `BENCHMARKS.md`
for the current numbers, `TUTORIAL.md` §8 for the target-class
reliability table that tells biologists when to trust CellSim
and when to stop.

## Campaigns 2–5

Campaign 2 (Mesoscale Cellular Modules) resumes after Campaign 1
exit. Campaign-2 pathway rates consume Campaign-1 cache output via
`src/bridge/`, so hand-tuning phenomenology becomes physics-prior
citation — the fix for the "closed ontology / circular validation"
critique that triggered the restart.

Campaigns 3 (Whole-Cell + 3D Tissue), 4 (Adversarial Validation +
Active Learning), 5 (Production / Open Deployment / Regulatory):
see `GOAL`.
