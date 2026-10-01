# src/bridge — Campaign-1 → Campaign-2 rate-law emitter

**Cross-campaign bridge, Campaign-1-scope only for now.**

## Scope
Given a cached FEP ΔG (or binding-affinity prediction with CI),
emit the matching Hill / Michaelis-Menten / mass-action rate law
with an uncertainty bar that Campaign-2's pathway models consume as
a physics-derived prior.

This is the module that retires the professor's "closed ontology /
circular validation" critique: Campaign-2 pathway rates will no
longer be hand-tuned; they will cite `src/bridge/` output with
method provenance and CI.

## What we write
Essentially all of this (no upstream equivalent):
- `binding_to_hill(deltaG, receptor_copies, ligand_conc) → (K_d, n, CI)`.
- `affinity_to_michaelis(kcat_pred, Km_pred, …) → (V_max, K_M, CI)`.
- Propagation of CI through standard rate-law compositions.
- Emission format that Campaign-2 config loaders can read directly
  (YAML with provenance block).

## Accuracy, not just precision (2026-07)
Every `RateLawPrior` carries `accuracy_kcalmol` + `accuracy_basis` +
`trust`, and the K_d CI is driven by the LARGER of the input
reproducibility σ and the measured target-class error. Rationale: the
input σ is *precision* ("does the calc agree with itself?"), which on
biotin/streptavidin was ±0.29 kcal/mol while the *true* error was
~11.8 — a Campaign-2 cell model handed the tight bar would trust a
wrong rate. Rules (see `_resolve_accuracy`):

- **Docking-sourced** ΔG uses `benchmarks/dock/reliability_table.yaml`
  via `src/uq/reliability.reliability_for(receptor)`. Untrustworthy
  classes widen the CI and set `trust` accordingly.
- **FEP-sourced** ΔG does *not* borrow the docking table (different
  method) and CellSim's absolute FEP is not yet GPU-validated, so it is
  `trust="uncalibrated"`, `accuracy_kcalmol=None` — flagged, never faked.
- **Unknown receptor** → `uncalibrated` (NEVER GUESS). Consumers must
  check `prior.accuracy_known` / `prior.trust` before treating a K_d as
  a calibrated absolute value.

## Exit test
For any compound in the cache, calling `bridge.rate_law(cid, target)`
emits a valid rate-law record. When Campaign 2 resumes, its first
pathway model ingests ≥ 10 bridged priors without manual tuning.

## Non-goal for Campaign 1
No actual pathway simulation lives here. That is Campaign 2.
