# src/uq — Uncertainty quantification

**Campaign 1, Layer 1.7.**

## Scope
Every Campaign-1 prediction emits `(value, 95 % CI, method,
provenance)` — nothing is allowed to ship naked. Conformal
prediction for docking ΔG; deep ensembles for ML-potential forces.

## Upstream tools (build-vs-buy: BUY)
- **MAPIE** (BSD-3) — conformal prediction (regression + classification).
- **NumPy / PyTorch** — hand-rolled deep-ensemble wrapper (trivial).

## What we write
- `Prediction` dataclass that enforces CI + method on every emission.
- MAPIE wrapper keyed to PDBBind/ChEMBL calibration splits.
- Deep-ensemble trainer for 3–5 MACE replicas (Layer 1.5 integration).
- Calibration-curve reporter for the blind harness.

## Exit test
Conformal calibration curves on held-out ChEMBL show nominal-vs-
empirical coverage within 5 %. No layer's output bypasses the
`Prediction` envelope.
