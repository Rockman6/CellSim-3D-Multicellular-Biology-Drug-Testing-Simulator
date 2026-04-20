# src/uq — Uncertainty quantification (non-AI)

**Campaign 1, Layer 1.6.**

## Scope
Every Campaign-1 prediction emits `(value, 95 % CI, method,
provenance)` — nothing ships naked. All UQ is mechanistic or
statistical; no neural ensembles.

## Methods

**Primary — mechanistic / statistical sampling:**
- **Parameter sweeps:** grid and Latin-hypercube scans over rate
  constants, force-field parameters, charges, geometric priors.
  Outputs a full prediction distribution, not a single value.
- **Monte-Carlo:** per-prediction MC over parameter distributions
  with seeded RNG; reports mean ± SD + quantile CIs.
- **Sobol global sensitivity** via **SALib**: Sobol first-order
  and total-effect indices identify which parameters actually move
  each headline prediction. Enables honest statements like "ΔG
  depends 60 % on force-field choice, 30 % on pose sampling,
  10 % on charge model."

**Secondary — post-hoc distribution-free statistical bounds:**
- **MAPIE conformal prediction** (`MAPIE`, BSD-3). Non-parametric
  wrapper that maps a held-out calibration set to distribution-free
  CI bounds on any predictor. Used *only* to attach statistical
  bounds to physics predictions; does not itself fit a model.
  Documented on every wrapped predictor as "statistical bound, not
  mechanistic insight — still depends on underlying calibration
  set quality."

## Explicitly not used
- Deep ensembles of learned potentials.
- Bayesian neural nets.
- Any neural surrogate for UQ.

## Upstream tools (build-vs-buy: BUY)
- **SALib** (MIT) — Sobol, Morris, FAST sensitivity analysis.
- **MAPIE** (BSD-3) — conformal prediction (regression +
  classification) as non-parametric post-hoc wrapper.
- **NumPy / SciPy** — parameter sweep drivers + quantile math.

## What we write
- `Prediction` dataclass that enforces CI + method on every emission.
- Parameter-sweep and Monte-Carlo drivers consuming any physics
  callable (`chem`, `md`, `dock`, `quantum`).
- SALib wrappers that produce Sobol indices per prediction class
  (docking ΔG, FEP ΔΔG, MD-derived observables).
- MAPIE wrapper keyed to PDBBind / ChEMBL calibration splits.
- Calibration-curve + tornado-sensitivity reporter for the blind
  harness.

## Exit test
- Conformal calibration curves on held-out ChEMBL show nominal-vs-
  empirical coverage within 5 %.
- Every headline prediction class (docking ΔG, FEP ΔΔG, CYP site
  score) ships with Sobol first-order indices published alongside.
- No layer's output bypasses the `Prediction` envelope.
