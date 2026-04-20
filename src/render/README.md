# src/render — Generic 3D rendering pipeline

**Cross-cutting, used by every Campaign-1 viewer.**

## Scope
Window + input + camera + Metal shader pipeline, stripped of any
cell-specific geometry. Provides molecular-primitive draws
(spheres, cylinders, ribbons, isosurfaces) plus a Dear ImGui panel
stack for live readouts and plots.

The OLD/ tree contains the original cell-specific render code
(`OLD/src/render/`, `OLD/src/gpu/`, `OLD/src/main.mm`). We salvage
its shader pipeline, camera, and MetalContext pattern; everything
that touched `SimCell` stays in OLD/.

## Upstream tools
- **Metal** (Apple system framework) — primary backend.
- **Dear ImGui** (MIT) — UI widgets, plots, buttons.
- **GLFW** (zlib) — cross-platform window / input (later).
- **GLM / simd** — 3D math.

## What we write
- `render/core/` — window, camera, MetalContext, shader loader.
- `render/mol/` — molecular primitives (sphere instancer, stick
  instancer, ribbon mesh generator, isosurface marcher).
- `render/plot/` — live time-series / scatter / heatmap widgets.

## Per-layer viewer targets (bound to Campaign-1 layers)
Per the 2026-04-20 non-AI amendment: 7 layers, no ML potential.
- 1.1 ligand ball-and-stick with per-atom charge colouring.
- 1.2 protein ribbon + solvent skin, RMSD gauge.
- 1.3 receptor + docked pose overlay, ΔG bar with CI (AutoDock Vina
  primary; GNINA CNN shown as labeled fast-guess only).
- 1.4 ESP isosurface + HOMO/LUMO density.
- 1.5 Martini bilayer + area-per-lipid live plot.
- 1.6 sensitivity tornado plot + Sobol index bars + calibration
  curve (SALib + MAPIE outputs).
- 1.7 blind-bench dashboard (pass/fail tiles + r-scatter + red-team
  leaderboard).

## Exit test
Each layer's numeric harness and its viewer are gated together: a
layer does not pass until the viewer renders correctly on the
layer's reference scene.
