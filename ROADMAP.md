# CellSim Roadmap

Authoritative strategic plan: **`GOAL`** (5-campaign hierarchical
multi-scale hybrid simulator, 4–7 years).

Approved execution plan for the current restart:
**`.claude-plans/make-this-plan-into-virtual-key.md`** (migration
manifest + Campaign-1 scope + build-vs-buy table).

## Where we are

**Pre-Campaign-1 restart, April 2026.**

The pre-restart HeLa / p53 / cisplatin Campaign-2 prototype is
frozen under `OLD/`. It still builds and passes its 6 headless
validators as a regression snapshot. The new top-level `src/` tree
is scaffolded but not yet populated; each layer's source is added in
a dedicated PR.

## Campaign 1 — Atomic → Molecular Foundation (Years 1–2, non-AI)

Per the 2026-04-20 non-AI amendment (see `MISSION.md`), Campaign 1
is 7 layers, no ML potentials. Original Layer 1.5 (MACE-OFF23) was
removed; the coarse-grained, UQ, and blind-validation layers
renumbered accordingly.

- 1.1 Chem foundation (RDKit + OpenFF + AM1-BCC)          — DONE
- 1.2 Classical MD (OpenMM + ff14SB + TIP3P)              — DONE (MVP)
- 1.3 Docking + FEP (AutoDock Vina primary + perses FEP)
- 1.4 Quantum (xtb + PySCF)
- 1.5 Coarse-grained (Martini 3)
- 1.6 UQ (Sobol + Monte-Carlo + parameter sweeps; MAPIE conformal
  as post-hoc statistical wrapper only)
- 1.7 Blind-validation harness (PDBBind / CASF / PoseBusters /
  ChEMBL)

Every layer ships **numeric harness + minimal viewer** in lockstep.
See `docs/campaign1_scope.md` for exit criteria.

## Campaigns 2–5

Campaign 2 (Mesoscale Cellular Modules) resumes after Campaign 1
exit. Campaign-2 pathway rates consume Campaign-1 cache output via
`src/bridge/`, so hand-tuning phenomenology becomes physics-prior
citation — the fix for the "closed ontology / circular validation"
critique that triggered the restart.

Campaigns 3 (Whole-Cell + 3D Tissue), 4 (Adversarial Validation +
Active Learning), 5 (Production / Open Deployment / Regulatory):
see `GOAL`.
