# src/viewer — Per-layer visual drivers

**Cross-cutting, one file per Campaign-1 layer.**

## Scope
Thin binding between a layer's streaming output (trajectory frames,
ΔG updates, calibration numbers) and `src/render/` primitives.
Invoked via `cellsim-chem --view <layer>`.

## Contents (to be populated layer-by-layer)
Per the 2026-04-20 non-AI amendment: 7 layers, no ML potential.
- `viewer_chem.cpp` — 1.1 ligand ball-and-stick.
- `viewer_md.cpp` — 1.2 protein trajectory.
- `viewer_dock.cpp` — 1.3 docking + ΔG (Vina primary; GNINA CNN
  as explicitly labeled fast-guess).
- `viewer_quantum.cpp` — 1.4 ESP + HOMO/LUMO.
- `viewer_cg.cpp` — 1.5 Martini bilayer.
- `viewer_uq.cpp` — 1.6 Sobol tornado + conformal calibration.
- `viewer_bench.cpp` — 1.7 blind-bench dashboard.

Each file is < 500 LOC and depends only on `src/render/` + the
shared-memory trajectory stream from the layer's numeric engine.
