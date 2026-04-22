# CellSim Mission

CellSim's mission is to build the world's best open-source,
multi-scale, uncertainty-aware cell simulator — reducing
triage-stage wet-lab drug-discovery burden by 10–100× on modelled
pathways, with calibrated uncertainty on every prediction,
validated continuously against external blind wet-lab data.

We do not claim replacement of wet-lab experiments. We claim
**triage reduction with honest uncertainty**. Every prediction
ships with method provenance and a calibrated confidence interval.
Every benchmark is publicly reproducible from a fresh clone. Every
claim is backed by external blind tests, not internal retrospective
fits.

## Primary metric

**Wet-lab compounds avoided per year** — tracked via partner-lab
collaborations. This is the only metric that matters for the
humanitarian goal of accelerating therapies.

## No black-box / no AI surrogates

CellSim is strictly physics-first and non-AI by design. This is a
load-bearing scientific commitment, not a preference, and it
governs every layer:

1. **All core processes are mechanistic.** Every biological or
   chemical step uses a deterministic or well-defined stochastic
   representation: ODEs, Gillespie SSA, tau-leaping, flux-balance
   LP, agent-based rules with physics-informed forces, mass-action
   kinetics, PDE solvers. No learned surrogate stands in for a
   mechanism.
2. **All forces are physics-derived or empirically parametrised.**
   Force fields (AMBER ff14SB, CHARMM36m, OpenFF Sage, Martini 3),
   semi-empirical / DFT (xtb, PySCF), alchemical FEP via
   `openmmtools.alchemy` + `pymbar` (non-AI MBAR estimator; the
   perses path was evaluated and not required — openmmtools
   primitives + a custom DDM driver cover the Milestone A/B
   scope) — yes. Neural potentials (MACE / OrbNet / NequIP /
   Allegro), CNN scoring (GNINA default mode) as the primary
   path — no.
3. **All uncertainty is deterministic or statistical, not
   neural.** Parameter sweeps, Monte-Carlo over rate constants and
   charges, Sobol global sensitivity indices (SALib) are the
   primary UQ. MAPIE conformal prediction is acceptable only as a
   *post-hoc non-parametric statistical wrapper* that maps a
   calibration set to distribution-free CI bounds; it provides
   statistical bounds, not mechanistic insight, and must be
   documented as such.
4. **Every prediction carries provenance and uncertainty** from
   the deterministic or statistical path above. No naked values.
5. **ML is allowed only as an accelerator of a well-validated
   subroutine**, and only when (a) a physics reference path
   actually runs alongside and is retained, (b) the ML output is
   cross-checked against the physics path on every production
   prediction, and (c) the prediction's provenance records which
   path produced it and the residual vs the reference path. ML is
   never the primary evidence. Never without an audit trail.
   Example: GNINA's CNN-scored pose prediction may ship as an
   explicitly labeled "fast-guess" mode alongside AutoDock Vina
   (primary); the CNN score is displayed next to the Vina score
   but never replaces it for final predictions.

Every rate constant cites a PMID or a specific physics calculation
hash from `src/cache/`. A grep for the rate value in the repo
should land on a literature citation or a cached calculation ID,
never on a hand-tuned magic number.

## Roadmap

Five overlapping campaigns over 4–7 years; see `GOAL` and
`ROADMAP.md`. We ship a usable deliverable at the end of every
campaign so the project has value even if later campaigns stall.
