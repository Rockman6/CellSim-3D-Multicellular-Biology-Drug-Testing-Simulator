#!/usr/bin/env bash
# CellSim — Milestone B Phase-2 EGFR kinase inhibitor GPU runner.
#
# Runs the 6-compound binding_egfr.yaml benchmark — same chemical
# series the docking layer got Spearman -0.49 on (anti-correlated
# with experiment). This is the test that retires the prof's
# "physics-FEP vs Vina on kinases" critique once it samples on
# GPU.
#
# Compounds (in increasing IC50 order):
#   erlotinib, AG-1478, lapatinib, gefitinib,
#   4-anilinoquinazoline, tyrphostin AG-494
# IC50 spans 2 nM -> 1.1 µM (~4 kcal/mol dynamic range).
#
# Target hardware: rented NVIDIA H100 / A100 (~$2-3/h cloud spot)
# or M-series Mac with Metal backend.
#
# Expected wall time at production parameters (11 windows × 25 000
# prod + 2 500 equil steps × 2 legs × 6 compounds, 110k-atom complex):
#   - H100:   ≈  2-3 h total
#   - A100:   ≈  4-5 h total
#   - M5 Max: ≈ 30-50 min total (per cellsim fep-binding validate
#                                 wall-time estimate)
#
# Lapatinib is the slowest (40 heavy atoms; AM1-BCC quartic scaling)
# at ~5 min scaffold + dominant share of sample time. Other 5
# compounds are 17-31 heavy atoms.
#
# Crash-recovery: --resume on bench preserves completed compounds
# across Ctrl-C / power loss / OOM kills. Restart with same script.

set -euo pipefail

if [ ! -f "benchmarks/fep/binding_egfr.yaml" ]; then
    echo "error: run from CellSim repo root (not $(pwd))" >&2
    exit 2
fi

# Activate env
if ! command -v conda >/dev/null 2>&1; then
    echo "error: conda/mamba not found. Install miniforge first:" >&2
    echo "  brew install miniforge" >&2
    exit 2
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cellsim

# Override via env: `OUT_DIR=run/fep/egfr_20260523_1234 bash <this>`
# resumes a crashed run via --resume on the bench step below
# (compounds with non-empty dG_pred_kcalmol are kept).
if [ -z "${OUT_DIR:-}" ]; then
    STAMP=$(date +%Y%m%d_%H%M%S)
    OUT_DIR="run/fep/egfr_${STAMP}"
fi
mkdir -p "${OUT_DIR}"
if [ -f "${OUT_DIR}/egfr_results.csv" ]; then
    echo "[run_binding] resuming from existing CSV in ${OUT_DIR}"
fi
CSV="${OUT_DIR}/egfr_results.csv"

# Header block — mirror to env.log so `cellsim fep-report` extracts
# `git commit:` for provenance. Previously stdout-only.
{
    echo "============================================================"
    echo "CellSim — EGFR kinase binding FEP (Milestone B flagship)"
    echo "============================================================"
    echo "  started:     $(date)"
    echo "  machine:     $(uname -a)"
    echo "  ram:         $(sysctl -n hw.memsize 2>/dev/null | awk '{print $1/1024/1024/1024 " GB"}' || echo '?')"
    echo "  git commit:  $(git rev-parse HEAD)"
    echo "  git ref:     $(git describe --tags --always 2>/dev/null || echo '?')"
    echo "  output:      ${OUT_DIR}/"
    echo ""
} | tee "${OUT_DIR}/env.log"

# Env + platform report.
python -c "
import openmm, openmmtools, openmmforcefields, openff.toolkit, pymbar, sys
from openmm import Platform
print('python            ', sys.version.split()[0])
print('openmm            ', openmm.__version__)
print('openmmtools       ', openmmtools.__version__)
print('openmmforcefields ', openmmforcefields.__version__)
print('openff-toolkit    ', openff.toolkit.__version__)
print('pymbar            ', pymbar.__version__)
print()
print('Available OpenMM platforms (fastest first):')
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print(f'  {p.getName()} (speed {p.getSpeed()})')
print()
print('Pipeline will prefer Metal -> CUDA -> OpenCL -> CPU.')
" 2>&1 | tee -a "${OUT_DIR}/env.log"
echo ""

echo "=== cellsim doctor ==="
./scripts/cellsim doctor 2>&1 | tee "${OUT_DIR}/doctor.log" | tail -5
DOCTOR_EXIT=${PIPESTATUS[0]}
if [ $DOCTOR_EXIT -ne 0 ]; then
    echo "error: cellsim doctor failed. Check ${OUT_DIR}/doctor.log" >&2
    exit 3
fi
echo ""

echo "=== fep-binding validate ==="
./scripts/cellsim fep-binding validate \
    benchmarks/fep/binding_egfr.yaml 2>&1 | tee -a "${OUT_DIR}/run.log"
VALIDATE_EXIT=${PIPESTATUS[0]}
if [ $VALIDATE_EXIT -ne 0 ]; then
    echo "error: validate FAILED (hard error). Fix YAML before sampling." >&2
    exit 4
fi
echo ""

echo "=== EGFR binding bench at $(date) ==="
echo "Parameters: 11 windows × 25 000 prod + 2 500 equil steps × 2 legs"
echo "            × 6 EGFR inhibitors on 1m17"
echo "            FF: amber14-all + tip3pfb + SMIRNOFF ligand"
echo "            Note: lapatinib (40 hvy) dominates wall-time"
echo ""

time ./scripts/cellsim fep-binding bench \
    benchmarks/fep/binding_egfr.yaml \
    --sample \
    --n-windows 11 \
    --production-steps 25000 \
    --equilibration-steps 2500 \
    --sample-stride 250 \
    --padding 1.2 \
    --force-field-path amber14 \
    --out-csv "${CSV}" \
    --resume \
    2>&1 | tee -a "${OUT_DIR}/run.log"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "============================================================"
echo "Bench finished at $(date)  (exit ${EXIT_CODE})"
echo "============================================================"
ls -la "${OUT_DIR}/"
echo ""

# Local verdict via fep-report. Binding YAML → 2.0 kcal/mol gate +
# every-ΔG-negative sign rule. Kendall τ ≥ 0.6 is the load-bearing
# rank-correlation claim against the docking-failure series.
echo "============================================================"
echo "Local verdict (cellsim fep-report)"
echo "============================================================"
if ./scripts/cellsim fep-report "${CSV}" \
        --yaml benchmarks/fep/binding_egfr.yaml \
        --out-dir "${OUT_DIR}" 2>&1 | tee -a "${OUT_DIR}/run.log"; then
    VERDICT_EXIT=0
else
    VERDICT_EXIT=$?
fi
echo ""

TARBALL="egfr_gpu_${STAMP}.tar.gz"
tar czf "${TARBALL}" -C run/fep "egfr_${STAMP}/"
echo "Created: ${TARBALL}  ($(du -h ${TARBALL} | cut -f1))"
echo ""
echo "Send this tarball back to Henry. The headline question for"
echo "the prof: did Kendall τ on this 6-compound EGFR series clear"
echo "the ≥ 0.6 gate (i.e. did physics-FEP rescue the docking-layer"
echo "Spearman -0.49 failure)?"
echo ""

FINAL_EXIT=${VERDICT_EXIT:-${EXIT_CODE}}
if command -v terminal-notifier >/dev/null 2>&1; then
    STATUS=$([ $FINAL_EXIT -eq 0 ] && echo "PASS" || echo "FAIL")
    terminal-notifier \
        -title "CellSim EGFR binding FEP" \
        -subtitle "${STATUS} (exit ${FINAL_EXIT})" \
        -message "Tarball ready: ${TARBALL}. Send it to Henry." \
        -sound Glass \
        -activate com.apple.Terminal 2>/dev/null || true
fi

exit ${FINAL_EXIT}
