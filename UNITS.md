# CellSim Unit Conventions

This file documents the physical units used in the simulation. Every
constant, variable, and rate should match these conventions.

## Time
- **Real time (`dt`)**: seconds of wall-clock time
- **Biological time (`bio_dt`)**: seconds of simulated biology
  - Conversion: `bio_dt = dt * BIO_TIME_SCALE`
  - Default `BIO_TIME_SCALE = 180` (1 real-second = 3 bio-minutes)
  - A 20-hour HeLa cell cycle completes in `20h * 3600 / 180 = 400 real-seconds`
- **Engine substeps**: each engine (Metabolism, Signaling, HH) internally
  subdivides `bio_dt` to meet its own stability limit.

## Space
- **Sim-units**: effectively micrometres (µm).
  - `CELL_RADIUS_BASE = 2.15` sim-units is a deliberately compressed size
    for viewport visibility. In "true" µm a HeLa cell is ~10 µm radius.
  - `SCENE_BOUND = 55` sim-units = a 110-µm-wide culture dish (small).
- Grid diffusion uses dimensionless lattice coefficients (not cm²/s).
  Diffusion speeds are tuned for visual plausibility, not quantitative match.

## Concentrations
- **Metabolism** (`Metabolism.h::MetabolitePool`): millimolar (mM).
  - ATP ≈ 3.0 mM, glucose ≈ 5 mM (physiological).
- **Signaling** (`Signaling.h::SignalingState`): micromolar (µM).
  - Cytosolic Ca²⁺ rest ≈ 0.1 µM, activated ≈ 1 µM.
- **Ion transport** (`MembraneTransport.h::IonState`): millimolar (mM).
- **Drug** (`Simulation.h::Drug`): micromolar (µM), matching published IC50s.

## Energy
- **Free energy (ΔG)**: kJ/mol.
- **Boltzmann temperature (kT)**: kJ/mol at 310 K ≈ 2.577 kJ/mol.
- **Membrane potential (Vm)**: millivolts (mV). Rest ≈ -70 mV (modern convention).

## Forces
- **Intracellular physics** (`IntracellularPhysics.h`): sim-units per second
  for velocities, sim-unit/s² for acceleration. Not calibrated to pN.
- **Cell-cell contact** (`Simulation.h`): Hertz contact, `F = k·overlap^1.5`
  where `k = HERTZ_STIFFNESS`. Again sim-unit scaled, not mN/m.

## Counts
- **Telomere length**: base pairs (bp). Loss per division 50-100 bp (Harley 1990).
- **DNA damage**: lesions per cell per bio-second. ~5000 depurinations/day
  (Lindahl 1993) → 0.058 lesions/s spontaneous.
- **Cell count**: integer cells in the scene.

## What is **dimensionless**
- Pathway "activity" fractions (GPCR_active, MEK_active, etc.) in Signaling.
  Range 0–1. Do NOT multiply these as if they were concentrations.
- CDK cyclin levels in `Simulation.h::CDKState` — a phenomenological
  oscillator, not biochemical concentrations. Use for phase classification only.
- `biomass`, `size`, `mitoHealth` — dimensionless health indicators, 0–2.3.

## Known inconsistencies (to fix)

1. `DIFF_O2_COEFF` in `Constants.h` is a grid-update coefficient, NOT a real
   diffusion constant. Labeled as such.
2. `HERTZ_STIFFNESS` is in sim-unit force per (sim-unit overlap)^1.5, not N/m.
3. `Metabolism.h` and `MembraneTransport.h` both use their own RT/F prefactor
   (log10 vs ln form). Both are correct but confusing.
4. Two ATP "sources of truth": `Simulation::updateMetabolism` (sets c.ATP 0-100)
   and `MetabolismEngine.pool.ATP` (mM). The cellBio engine overwrites cell[0]
   ATP every frame via `getATPpercent()`.
