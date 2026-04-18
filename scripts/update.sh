#!/usr/bin/env bash
# ════════════════════════════════════════════════════════════════════════
#  CellSim incremental update
#  ----------------------------------------------------------------------
#  Pulls only the latest CHANGED files (source + shaders + assets meta)
#  from the upstream GitHub repo and rebuilds. No full re-clone needed.
#
#  Usage:
#      scripts/update.sh                # update + rebuild
#      scripts/update.sh --constants    # only refresh Constants.h + tuning,
#                                       # then rebuild
#      scripts/update.sh --check        # dry-run: show what would change
#
#  Designed so frequent constant tweaks (rate calibrations, biology tunes,
#  thresholds) ship as a tiny patch — no need to re-download GLB assets,
#  PDFs, or the molecule cache that already lives on the user's machine.
# ════════════════════════════════════════════════════════════════════════

set -e
REPO_URL="https://github.com/Rockman6/CellSim-3D-Multicellular-Biology-Drug-Testing-Simulator.git"
BRANCH="main"
MODE="${1:---all}"

cd "$(dirname "$0")/.."   # repo root

# Make sure we're in a git checkout. If not, the user cloned manually
# without git — bail with helpful message.
if [ ! -d .git ]; then
    echo "[update] not a git repo — please re-clone with:"
    echo "         git clone $REPO_URL"
    exit 1
fi

case "$MODE" in
    --constants)
        echo "[update] fetching only constants + tuning files..."
        git fetch origin "$BRANCH"
        # Sparse-pull only the constants / biology rate files. Everything
        # else (assets, GLBs, shaders) stays local — these are the small
        # files that change most often during calibration.
        git checkout "origin/$BRANCH" -- \
            src/core/Constants.h \
            src/simulation/Simulation.h \
            src/simulation/CentralDogma.h \
            src/simulation/CellCycleProgram.h \
            src/simulation/MediumField.h \
            src/simulation/TelemetryLog.h
        ;;
    --check)
        echo "[update] dry-run — files that would change:"
        git fetch origin "$BRANCH"
        git diff --stat HEAD "origin/$BRANCH"
        exit 0
        ;;
    --all|*)
        echo "[update] pulling latest source from origin/$BRANCH..."
        git pull --ff-only origin "$BRANCH"
        ;;
esac

# Rebuild — incremental, only changed files recompile.
echo "[update] rebuilding (incremental)..."
cmake --build build -j

echo ""
echo "[update] done. Launch with: ./build/CellSim"
echo "         or run validation: ./build/cellsim_headless 48 60"
