# OLD — the 2026 C++/Metal prototype (frozen)

This is CellSim v1.0–v1.2 (March–April 2026): a real-time 3D cell
simulator in C++20 with an Apple Metal renderer, an ImGui panel and drug
PK/PD. It was frozen here on 2026-04-20 when the project restarted.

**Status.** Kept as a reference implementation of the pathway topology
and its literature citations, and as a regression snapshot. Its
interface is retired; a new web interface is planned on top of the
Python engine (see `../docs/PLAN.md`). The cell cycle, the p53 axis and
the apoptosis cascade have been ported to `../cellsim/cell/engine.py`,
which is where new biology goes.

## Build the headless validators

No vcpkg or GUI dependencies are needed for the six validators. From the
repository root:

```bash
cmake -S OLD -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/cellsim_blind_benchmarks            # BM01–BM08
./build/cellsim_p53_pulsing
./build/cellsim_cisplatin_arrest
./build/cellsim_cisplatin_apoptosis
```

The GUI target (`CellSimOLD`) builds only when vcpkg provides glfw3,
glm, imgui and implot; Metal makes it macOS-only.

Re-run on 2026-10-01: 8/8 blind benchmarks pass at seeds 0x1 and
0xDEADBEEF, and the three single validators pass. What those passes do
not show is recorded in `../docs/VALIDATION.md`: no cell ever dies under
cisplatin (caspase-3 stays at 0.04 under saturated MOMP), and the p53
pulse period is 3.0 h against 5.5 h in the literature. Both are fixed in
the Python port.

## Layout

| Path | What it is |
|---|---|
| `src/` | the simulator (`simulation/Simulation.h` is the core; `biochem/` the pathway modules; `gpu/`, `render/` the Metal renderer) |
| `tests/` | the six headless validators |
| `scripts/` | calibration, validation, capture and molecule-generation scripts for this prototype. Run them from the repository root after configuring `build/` from `OLD/` as above. Unmaintained. |
| `docs/` | the LaTeX technical reference and equations book for this code |
| `UNITS.md` | the prototype's unit conventions and known inconsistencies |

Assets (`../assets/`: GLB organelle meshes, Metal shaders, molecule and
protein files) and the bioagent catalogue (`../data/bioagents/`) are
still read from the repository root by this build.
