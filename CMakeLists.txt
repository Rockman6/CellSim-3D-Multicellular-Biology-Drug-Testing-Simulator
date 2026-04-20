cmake_minimum_required(VERSION 3.24)
project(CellSim VERSION 2.0 LANGUAGES CXX)

# ══════════════════════════════════════════════════════════════════════════
#  CellSim — Hierarchical Multi-Scale Hybrid Cell Simulator
#
#  Restart root. Campaign-1 (Atomic → Molecular Foundation) target set
#  lives below. The pre-restart Campaign-2 prototype (HeLa agent model +
#  p53/MDM2/apoptosis + cisplatin/Nutlin benchmarks) is preserved in
#  OLD/ as a frozen regression snapshot.
#
#  See:
#    GOAL                    — 5-campaign roadmap
#    MISSION.md              — mission sentence
#    ROADMAP.md              — where we are, where we're going
#    docs/campaign1_scope.md — Campaign 1 scope + exit criteria
# ══════════════════════════════════════════════════════════════════════════

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(CELLSIM_BUILD_OLD "Build the Campaign-2 prototype snapshot in OLD/" ON)
option(CELLSIM_BUILD_CAMPAIGN1 "Build Campaign-1 targets (currently stubs)" ON)

# ── Campaign-1 scaffolding ──────────────────────────────────────────────
#
#   src/chem     — RDKit + OpenFF + charges         (Layer 1.1)
#   src/md       — OpenMM + MACE driver             (Layer 1.2, 1.5)
#   src/quantum  — xtb + PySCF                      (Layer 1.4)
#   src/cg       — Martini 3                        (Layer 1.6)
#   src/cache    — physics-prior cache (SQLite+HDF5)
#   src/uq       — conformal + deep ensembles       (Layer 1.7)
#   src/bridge   — Campaign-1 → Campaign-2 rate-law emitter
#   src/render   — generic Metal/ImGui pipeline     (viewers)
#   src/viewer   — per-layer visual drivers
#   src/core     — physics-neutral utilities (RNG, telemetry, constants)
#
# Each layer ships binary + viewer as siblings; see src/*/README.md.
# No source files exist yet — PRs populating each layer will add them
# and register their target here. This file is intentionally minimal.

if(CELLSIM_BUILD_CAMPAIGN1)
    message(STATUS "CellSim Campaign-1: scaffolding present; no targets built yet.")
    message(STATUS "  Next PR implements Layer 1.1 (RDKit/OpenFF) in src/chem/.")
endif()

# ── OLD/ regression snapshot ────────────────────────────────────────────
if(CELLSIM_BUILD_OLD AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/OLD/CMakeLists.txt")
    message(STATUS "CellSim OLD: building Campaign-2 prototype snapshot.")
    add_subdirectory(OLD)
endif()

message(STATUS "CellSim root: project v${PROJECT_VERSION} (restart).")
