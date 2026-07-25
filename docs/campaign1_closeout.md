# Campaign 1 — close-out scorecard + molecule→cell bridge

**Date:** 2026-07. **Purpose:** record the verified state of Campaign 1's
eight exit criteria, the decisions taken on the ones the physics-only
method cannot meet, and scope the bridge from Campaign-1 molecular
outputs to a Campaign-2 cell model — the actual end goal.

This file is a summary. The authoritative per-criterion text (with the
recorded amendments) lives in `campaign1_scope.md`; the decision
constants in `provenance_audit.md`; the numbers in `BENCHMARKS.md`.

## Scorecard (verified)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| 1 | Pose recovery ≥ 75 % (best-of-top-3) | ✅ **PASS — 86.7 %** (13/15) | `benchmarks/pdbbind/blind_set_results.summary.json` |
| 3 | PoseBusters physical validity ≥ 95 % | ✅ **PASS — 100 %** (15/15) | same summary |
| 4 | UQ calibration + Sobol | 🟡 **substantively met** | coverage measured (`run_uq_coverage.py`), per-class reliability table, Sobol insensitivity (`sobol_sensitivity.json`). Gap: n < 19/class ⇒ not yet a *guaranteed* conformal bound |
| 6 | Provenance / no magic numbers | ✅ **PASS** | `provenance_audit.md` |
| 2 | Kinase IC50 ranking r ≥ 0.7 | ⛔ **bounded limitation** | Vina cannot rank kinases (EGFR `rank_order_only`, negative Spearman). Re-scoped: docking generates, FEP ranks. Deferred to FEP milestone (needs GPU) |
| 5 | CYP3A4 SoM ≥ 15/20 | ⛔ **bounded limitation** | ~2/3 sampled; systematic xTB-BDE blind spot on N-dealkylation. Tool re-scoped to **advisory** |
| 7 | Docker fresh-clone reproducibility | ⚪ **deferred engineering** | no Dockerfile yet; not science |
| 8 | Per-layer viewers render | ⚪ **deferred engineering** | `src/render/`, `src/viewer/` are README-only |

### Reading of the scorecard

The four criteria that test the **core molecular science** — does docking
*find* the right pose (1), are the poses *physically real* (3), are the
*uncertainties honest* (4), is every number *traceable* (6) — are met.
That is the foundation Campaign 2 stands on, and it is sound.

The two ⛔ criteria are **genuine, characterised limitations of the
physics-only fast path**, not unfinished work:

- **#2 kinase ranking.** Fast empirical docking rank-orders potent
  ATP-site binders poorly. This is a known property, the same one behind
  criterion 1's top-1 vs top-3 gap. The pipeline's answer is division of
  labour: **Vina generates poses, alchemical FEP ranks them.** Kinase
  ranking is an FEP-milestone deliverable, gated on GPU time.
- **#5 CYP3A4 metabolism.** GFN2-xTB C–H bond-dissociation energies do
  not resolve N-dealkylation sites, so tertiary-amine substrates are
  systematically mis-ranked. The site-of-metabolism tool is therefore
  **advisory** — it flags plausible C–H oxidation sites and is documented
  as not a validated 15/20 predictor.

Both are recorded honestly with a known cause and a known (deferred) fix,
rather than hidden behind a tuned number.

The two ⚪ criteria (#7 Docker, #8 viewers) are **engineering, not
science.** They are deferred so effort goes to the real goal — the cell
simulation — and are tracked here as remaining work, not as met.

**Decision (2026-07):** Campaign 1's molecular-science foundation is
validated to the bar its methods can honestly reach. Rather than
rabbit-hole on #2/#5 (research problems) or #7/#8 (engineering), begin
bridging Campaign-1 molecular outputs toward the Campaign-2 cell model,
and carry #2/#5/#7/#8 forward as documented, non-blocking remaining work.

## The molecule → cell bridge

The whole point of Campaign 1 is to feed Campaign 2 (Mesoscale Cellular
Modules) with **physics-derived rate constants that carry provenance and
calibrated uncertainty**, so cell models cite a binding calculation
instead of hand-tuning a rate — the fix for the "circular validation"
critique that triggered the project restart.

### What exists today — `src/bridge/`

`RateLawPrior` (a provenance-tracked record) plus two emitters:

- `binding_to_hill(dG, σ)` — absolute binding ΔG → K_d = exp(ΔG/RT) →
  Hill-equation prior; the ΔG 1σ propagates multiplicatively to a 95 %
  CI on K_d. Units are SI molar so a Campaign-2 ODE plugs them in
  directly.
- `affinity_to_michaelis(kcat, K_M)` — enzyme-kinetics → Michaelis-Menten
  prior (pass-through; FEP does not emit kcat directly).

Each record carries `source` / `method` / `temperature_K` / timestamp so
the audit trail says exactly which calculation produced each rate.

### Gaps to reach a cell (the real next work)

**Gap A — no consumer.** The bridge *emits* rate-law priors; nothing yet
*consumes* them. Campaign 2 needs the first cell-level model: a small
reaction/ODE network whose rate constants are `RateLawPrior` records, so
we can watch a molecular ΔG propagate to a cell-level readout with its
uncertainty intact. This is the next concrete build.

**Gap B — the bridge propagates precision, not accuracy.**
`binding_to_hill` turns the MBAR/seed σ into the K_d CI. But per
criterion 4, that σ is *reproducibility*, not *correctness*: the real
error is target-class dependent (0.9 – 5 kcal/mol; biotin's seed bar
understated true error ~41×). A rate-law prior that carries only the
MBAR σ would hand Campaign 2 an **over-confident K_d** — the exact trap
fixed one layer down in docking. **The bridge must attach the
target-class reliability verdict** (`src/uq/reliability.py`) to every
`RateLawPrior`, and widen or flag the CI when the source calculation is
in an untrustworthy class. Until then, a cell model could inherit a tight
bar on a wrong rate. *(High-value, self-contained fix; do first.)*

**Gap C — only two edge types exist, and the missing ones are exactly
the questions already raised about the cell:**

- **Two drugs interacting / competing at one site** ("will two drugs
  stick and change the ruler?") — needs a *competitive-binding* edge:
  two `RateLawPrior` K_d's feeding a shared-site occupancy term, and
  (for covalent/complex formation) a reaction that makes a new species.
  Not yet modelled.
- **Concentration varying across the membrane** ("measurements differ by
  how much drug is around that part") — needs a *spatial concentration
  field* + a *membrane-permeability* edge. Layer 1.5 (Martini bilayer)
  is the physics source for the permeability coefficient; the bridge
  would emit a transport rate, and the Campaign-2 model would carry
  per-compartment concentration. Not yet wired.
- **Drug inside vs outside the cell** — needs at least two
  *compartments* (extracellular / cytoplasm) with the permeability edge
  connecting them, so occupancy is computed at the *local* concentration,
  not a single bulk number. Not yet modelled.

These three are the substance of Campaign 2's first module and map
one-to-one onto real questions about how a drug actually acts on a cell.
They are named here so the roadmap is explicit, not so they are claimed
as built.

### Suggested order

1. **Gap B** — attach reliability to `RateLawPrior` (small, high-value,
   applies the criterion-4 lesson; keeps the cell model honest from
   day one).
2. **Gap A** — a minimal single-compartment occupancy model consuming a
   `RateLawPrior`, with the uncertainty carried end-to-end and rendered
   — this also gives criterion 8 its first real viewer.
3. **Gap C** — add compartments + permeability (uses Layer 1.5) and
   competitive/multi-drug edges, in that order.
