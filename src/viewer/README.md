# src/viewer — Per-layer visual drivers

**Cross-cutting, one file per Campaign-1 layer.**

## Scope
Thin binding between a layer's streaming output (trajectory frames,
ΔG updates, calibration numbers) and `src/render/` primitives.
Invoked via `cellsim-chem --view <layer>`.

## Contents (to be populated layer-by-layer)
- `viewer_chem.cpp` — 1.1 ligand ball-and-stick.
- `viewer_md.cpp` — 1.2 protein trajectory.
- `viewer_dock.cpp` — 1.3 docking + ΔG.
- `viewer_quantum.cpp` — 1.4 ESP + HOMO/LUMO.
- `viewer_ml.cpp` — 1.5 MACE vs classical.
- `viewer_cg.cpp` — 1.6 Martini bilayer.
- `viewer_uq.cpp` — 1.7 calibration.
- `viewer_bench.cpp` — 1.8 blind-bench dashboard.

Each file is < 500 LOC and depends only on `src/render/` + the
shared-memory trajectory stream from the layer's numeric engine.
