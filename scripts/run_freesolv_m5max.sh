#!/usr/bin/env bash
# CellSim — Milestone A FreeSolv FEP gate runner for Apple-silicon Mac.
# One command, walk away, come back to a tarball to send back to Henry.
#
# Target hardware: M1 Pro / Max / Ultra, M2, M3, M4, M5 Max. ≥ 32 GB
# unified memory recommended; 16 GB works but may swap.
#
# Expected wall time (from CPU grid extrapolation + OpenCL 5× speedup):
#   - 4-6 hours on M5 Max
#   - 6-10 hours on M1/M2 base
# Plugged in + display-sleep-only (not full sleep) required.
#
# The script runs the Milestone A gate per the professor's acceptance
# criteria: all 12 FreeSolv compounds at production parameters
# (50 ps equilibration + 50 ps production per window × 11 windows ×
# 2 legs = ~2.2 ns simulated MD per compound). Emits:
#   - per-compound predicted ΔG_hyd with MBAR uncertainty
#   - per-compound per-window GHMC acceptance rates
#   - aggregate MAE vs experimental (FreeSolv gate: MAE ≤ 1.5 kcal/mol)

set -euo pipefail

if [ ! -f "benchmarks/fep/freesolv_12.yaml" ]; then
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

# Output dir, timestamped so repeat runs don't collide. Override
# via env: `OUT_DIR=run/fep/20260523_100401 bash <this>` resumes
# a crashed run from where it left off (any compound with a
# non-empty dG_pred_kcalmol in freesolv_results.csv is skipped).
if [ -z "${OUT_DIR:-}" ]; then
    STAMP=$(date +%Y%m%d_%H%M%S)
    OUT_DIR="run/fep/${STAMP}"
fi
mkdir -p "${OUT_DIR}"
if [ -f "${OUT_DIR}/freesolv_results.csv" ]; then
    echo "[run_freesolv] resuming from existing CSV in ${OUT_DIR}"
fi

# Header block — mirror to env.log so `cellsim fep-report` can
# extract the `git commit:` line for provenance. Earlier versions
# of this script sent the header to stdout only, which left
# report.env_log metadata empty on the tarball.
{
    echo "============================================================"
    echo "CellSim — FreeSolv FEP gate (Milestone A)"
    echo "============================================================"
    echo "  started:     $(date)"
    echo "  machine:     $(uname -a)"
    echo "  ram:         $(sysctl -n hw.memsize 2>/dev/null | awk '{print $1/1024/1024/1024 " GB"}' || echo '?')"
    echo "  git commit:  $(git rev-parse HEAD)"
    echo "  git ref:     $(git describe --tags --always 2>/dev/null || echo '?')"
    echo "  output:      ${OUT_DIR}/"
    echo ""
} | tee "${OUT_DIR}/env.log"

# Env + platform report — critical for reproducibility. Appended
# to env.log (not overwriting the header we just wrote).
python -c "
import openmm, openmmtools, openff.toolkit, pymbar, sys
from openmm import Platform
print('python       ', sys.version.split()[0])
print('openmm       ', openmm.__version__)
print('openmmtools  ', openmmtools.__version__)
print('openff-toolkit', openff.toolkit.__version__)
print('pymbar       ', pymbar.__version__)
print()
print('Available OpenMM platforms (fastest first):')
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print(f'  {p.getName()} (speed {p.getSpeed()})')
print()
print('Pipeline will prefer Metal → OpenCL → CUDA → CPU.')
" 2>&1 | tee -a "${OUT_DIR}/env.log"
echo ""

# Cellsim doctor — fail fast if env is broken
echo "=== cellsim doctor ==="
./scripts/cellsim doctor 2>&1 | tee "${OUT_DIR}/doctor.log" | tail -5
DOCTOR_EXIT=${PIPESTATUS[0]}
if [ $DOCTOR_EXIT -ne 0 ]; then
    echo "" >&2
    echo "error: cellsim doctor failed. Check ${OUT_DIR}/doctor.log" >&2
    echo "Do NOT proceed with the gate run until doctor passes." >&2
    exit 3
fi
echo ""

# Run the full FreeSolv gate
echo "=== FreeSolv FEP gate starting at $(date) ==="
echo "Parameters: 11 windows × 25 000 production steps (50 ps at 2 fs)"
echo "            ×  5 000 equilibration steps (10 ps)"
echo "            × 2 legs (vacuum + TIP3P-solvated)"
echo "            × 12 FreeSolv compounds (methane → acetamide)"
echo ""

time ./scripts/cellsim fep-freesolv \
    benchmarks/fep/freesolv_12.yaml \
    --n-windows 11 \
    --production-steps 25000 \
    --equilibration-steps 5000 \
    --sample-stride 250 \
    --out-csv "${OUT_DIR}/freesolv_results.csv" \
    --resume \
    2>&1 | tee -a "${OUT_DIR}/run.log"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "============================================================"
echo "Finished at $(date)  (exit ${EXIT_CODE})"
echo "  0      = MAE ≤ 1.5 kcal/mol (gate passed)"
echo "  non-0  = gate failed or some compounds NaN'd"
echo "============================================================"
echo ""
ls -la "${OUT_DIR}/"
echo ""

# Local verdict preview — run cellsim fep-report on the result so
# you see PASS/FAIL/partial before WeChat'ing the tarball. Also
# writes report.md + parity.png INTO the output dir so the tarball
# Henry receives has the verdict baked in. Uses --yaml to auto-set
# expected-rows so a killed-early run is flagged 'partial'.
echo "============================================================"
echo "Local verdict preview (cellsim fep-report)"
echo "============================================================"
if ./scripts/cellsim fep-report "${OUT_DIR}" \
        --yaml benchmarks/fep/freesolv_12.yaml \
        --out-dir "${OUT_DIR}" 2>&1 | tee -a "${OUT_DIR}/run.log"; then
    VERDICT_EXIT=0
else
    VERDICT_EXIT=$?
fi
echo ""

# Bundle for transfer — after the analyser so report.md + parity.png
# are inside the tarball.
TARBALL="freesolv_m5max_${STAMP}.tar.gz"
tar czf "${TARBALL}" -C run/fep "${STAMP}/"
echo "Created: ${TARBALL}  ($(du -h ${TARBALL} | cut -f1))"
echo ""
echo "Send this tarball back to Henry (WeChat / email / Google Drive /"
echo "AirDrop). It now contains report.md + parity.png with the verdict"
echo "ALREADY computed — Henry just opens them to read the numbers."
echo ""

# Notify on completion (silent if terminal-notifier isn't installed).
# Report-exit takes precedence over the gate's raw exit: a partial
# run (exit 1 from --expected check) is a more useful status than
# 'gate passed per compound present'.
FINAL_EXIT=${VERDICT_EXIT:-${EXIT_CODE}}
if command -v terminal-notifier >/dev/null 2>&1; then
    STATUS=$([ $FINAL_EXIT -eq 0 ] && echo "PASS" || echo "FAIL")
    terminal-notifier \
        -title "CellSim FreeSolv FEP" \
        -subtitle "${STATUS} (exit ${FINAL_EXIT})" \
        -message "Tarball ready: ${TARBALL}. Send it to Henry." \
        -sound Glass \
        -activate com.apple.Terminal 2>/dev/null || true
fi

exit ${FINAL_EXIT}
