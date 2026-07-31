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

**Gap C — the missing edge types are exactly the questions already
raised about the cell. All three now have a first, honest model:**

- **Concentration varying across the membrane / drug inside vs outside**
  ("measurements differ by how much drug is around that part",
  "drugs inside the cell") — ✅ `src/cell/compartments.py`. Two
  compartments (extracellular reservoir ↔ cytoplasm), passive permeation
  `C_in(t) = C_out(1 − e^{−t/τ})`, `τ = V/(P·A) = r/3P`. Occupancy is
  computed at the *local* `C_in`, and the naive bulk-`C_out` reading is
  kept alongside so the barrier's effect is explicit. The permeability P
  carries provenance + `trust` (would come from a Martini PMF/diffusion
  profile, Layer 1.5 — not wired yet, so `uncalibrated` until then).
- **Two drugs competing at one site** ("will two drugs change the
  ruler?") — ✅ `src/cell/competition.py`. Single-site competitive
  equilibrium `θᵢ = (Lᵢ/Kdᵢ)/(1 + Σ Lⱼ/Kdⱼ)`; rigorous interval CI over
  the K_d CIs (θ is monotonic in each K_d).
- **Two drugs sticking to form a new species** ("stick and form a new
  drug, how strongly they stick") — ✅ `src/cell/complexation.py`.
  Bimolecular `A + B ⇌ AB` solved as the standard quadratic; the "how
  strongly" ruler is the pair's K_d prior; complex CI is exact over K_d.

Every one of these carries the `trust` verdict through to the readout, so
a cell-level number built on an untrustworthy K_d (or an uncalibrated
permeability) is never mistaken for decision-grade.

**Time-course dynamics — ✅ `src/cell/dynamics.py`.** The coupled
permeation + reversible-binding ODE, integrated with a stiff-aware
adaptive solver (LSODA). Kinetics come from thermodynamics: K_d fixes
k_off = k_on·K_d, with k_on a documented physics input (diffusion-limited,
Berg–Purcell-capped). Rigor anchors: (1) the ODE steady state PROVABLY
equals the analytic Hill occupancy at L_out (consistency test, exact);
(2) mass is conserved to ~1e-15; (3) the EQUILIBRIUM trust depends on K_d
alone and is kept separate from the uncalibrated k_on *timescale*, so an
unmeasured on-rate never taints the equilibrium verdict.

**Joint multi-parameter uncertainty — ✅ `src/cell/montecarlo.py`.**
`monte_carlo_propagate` samples each prior's binding ΔG (Gaussian; the
EFFECTIVE σ = max of precision and measured accuracy, so it carries
accuracy) → K_d = exp(ΔG/RT), runs any cell readout, and returns the
output distribution + percentile CI + worst-of-inputs trust. Deterministic
given a seed. Rigor: for a single-prior monotone readout its CI AGREES
with the exact interval-arithmetic CI the analytic modules compute (mutual
cross-check); for multiple correlated inputs / non-monotone readouts —
where the per-module interval is not valid — the MC is the correct tool.

**pH partitioning / ion trapping — ✅ `src/cell/partitioning.py`.** Exact
Henderson-Hasselbalch: only the neutral form crosses, so at steady state
the TOTAL concentration differs between compartments of different pH —
weak bases accumulate in acidic compartments (lysosomes), weak acids in
basic ones. `accumulation R = f_neutral(pH_out)/f_neutral(pH_in)`, exact
given (pKa, pHs, ion type). Cited physiological pH constants
(`COMPARTMENT_PH`). This is the first module where C_in ≠ C_out at steady
state — demo: a pKa-9 base targeting a lysosomal enzyme accumulates 489×,
lifting predicted occupancy from a naive 6.8 % (blood conc) to 97 %.

**Active efflux transport — ✅ `src/cell/transport.py`.** A saturable
Michaelis-Menten efflux pump (P-gp / ABCB1) that holds the intracellular
concentration BELOW C_out — the central mechanism of multidrug
resistance. Steady state solves `k_perm(C_out − C_in) = V_max·C_in/(K_m +
C_in)` in closed form (positive root of the quadratic). Because the pump
is saturable, its effect is dose-dependent: demo shows C_in held 17.7×
low at 1 nM (0 % occupancy — drug defeated) but the pump 99 % saturated at
10 µM so C_in recovers (97 % occupancy). K_m can come from the
`bridge.affinity_to_michaelis` prior (now it has a consumer).

**Still not modelled (next):** an intracellular *binding sink* buffering
the free concentration, and diprotic partitioning (chloroquine-class,
traps even harder). Composition of the transport / trapping / binding
effects into one steady state is also open (each is exact alone).

### Order taken

1. **Gap B — DONE.** Reliability attached to `RateLawPrior`.
2. **Gap A — DONE.** Single-compartment occupancy readout.
3. **Gap C — DONE (first models).** Compartments + permeability,
   competitive binding, and drug–drug complexation, all with trust +
   CI carried end-to-end.

A real viewer (criterion 8) rendering one of these scenes is the natural
next step; the numeric layer it would visualise now exists.
