# src/cell — Campaign-2 single-cell drug disposition

**The first Campaign-2 module.** Consumes `src/bridge` rate-law priors
(physics-derived K_d + accuracy + trust) and answers cell-level questions:
how much drug gets inside, what it does there, and how sure we are.

## Modules

| File | Physics |
|---|---|
| `occupancy.py` | Hill occupancy θ = L^n/(K_d^n + L^n) — fraction of target bound |
| `compartments.py` | Passive permeation, C_in(t) = C_out(1−e^{−t/τ}), τ = V/(P·A) = r/3P |
| `partitioning.py` | pH ion trapping (Henderson–Hasselbalch); weak bases accumulate in lysosomes |
| `transport.py` | Saturable efflux pump (P-gp) — multidrug resistance |
| `binding_sink.py` | Nonspecific intracellular binding — free vs total drug |
| `steady_state.py` | **Composes** permeation + pH + efflux + sink into one steady state |
| `dynamics.py` | Permeation + reversible binding ODE (kinetics from thermodynamics) |
| `disposition.py` | **Full-system transient** — converges to `steady_state` |
| `competition.py` | Two drugs contending for one site |
| `complexation.py` | Two drugs binding each other (A + B ⇌ AB) |
| `montecarlo.py` | Joint uncertainty propagation through any readout |
| `viewer.py` | Criterion-8 reference scene |

## Two invariants everything obeys

1. **Trust rides through.** Every readout inherits the weakest input's
   verdict. A cell-level number built on a `do_not_trust_absolute` K_d or
   an uncalibrated permeability is never reported as decision-grade. The
   uncertainty is *accuracy*-based (measured target-class error), not just
   reproducibility — see `src/bridge`.
2. **Composed models reduce to their parts.** `steady_state` provably
   reduces to each standalone module in its limit, and `disposition`
   converges to `steady_state`. These are tests, not claims.

## Reference scene (criterion 8)

```bash
python src/cell/viewer.py benchmarks/dock/cell_disposition_scene.png
```

Four panels, all computed from the validated engine above:

- **A** dose–response with a Monte-Carlo CI band, passive vs P-gp efflux;
- **B** transient on log time — the binding sink as a *capacitor* (same
  plateau, t½ shifted ~3 decades by sink capacity);
- **C** pH ion trapping across compartments (weak base vs weak acid);
- **D** composed steady state: permeation 28 % → +lysosomal trapping 99 %
  → +P-gp efflux 0.1 %.

Palette is Okabe–Ito (colourblind-safe), assigned in fixed order.

## Not modelled yet

Diprotic partitioning (chloroquine-class); multi-cell / spatial gradients;
metabolism-driven clearance.
