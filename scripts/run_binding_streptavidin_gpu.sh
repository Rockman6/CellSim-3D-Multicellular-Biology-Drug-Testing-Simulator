#!/usr/bin/env bash
# CellSim — Milestone B Phase-2 streptavidin-biotin GPU runner.
#
# Runs the 4-compound binding_streptavidin.yaml benchmark on a CUDA
# / OpenCL / Metal GPU and emits a tarball containing the per-
# compound ΔG_bind CSV, the auto-generated PASS/FAIL report,
# the parity PNG, and the env / doctor / run logs.
#
# Target hardware:
#   - rented NVIDIA H100 / A100 (~$2/h cloud spot)
#   - or M-series Mac with Metal backend
#
# Expected wall time at production parameters (11 windows × 25 000
# prod + 2 500 equil steps × 2 legs per compound × 4 compounds):
#   - H100:        ≈  1 h total
#   - A100:        ≈  2 h total
#   - M5 Max:      ≈ 15-30 min total (GPU detect; mostly waters)
#
# Crash-recovery: --resume on the bench command means a Ctrl-C or
# power loss preserves completed compounds; restart with the same
# script and they're skipped automatically.
#
# Pre-flight: cellsim doctor must pass; bench-all must show this
# YAML clean. The script aborts on either failure.

set -euo pipefail

if [ ! -f "benchmarks/fep/binding_streptavidin.yaml" ]; then
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

# Output dir, timestamped (resume across runs by re-using the same
# CSV path). New STAMP per launch but a stable CSV name lets a
# Ctrl-C'd run resume by re-launching with the same date — see
# --resume note in the bench section below.
STAMP=$(date +%Y%m%d_%H%M%S)
OUT_DIR="run/fep/streptavidin_${STAMP}"
mkdir -p "${OUT_DIR}"
CSV="${OUT_DIR}/streptavidin_results.csv"

# Header block — mirror to env.log so `cellsim fep-report` extracts
# `git commit:` for provenance. Previously stdout-only, which made
# report.git_commit=None on the tarball.
{
    echo "============================================================"
    echo "CellSim — Streptavidin binding FEP (Milestone B)"
    echo "============================================================"
    echo "  started:     $(date)"
    echo "  machine:     $(uname -a)"
    echo "  ram:         $(sysctl -n hw.memsize 2>/dev/null | awk '{print $1/1024/1024/1024 " GB"}' || echo '?')"
    echo "  git commit:  $(git rev-parse HEAD)"
    echo "  git ref:     $(git describe --tags --always 2>/dev/null || echo '?')"
    echo "  output:      ${OUT_DIR}/"
    echo ""
} | tee "${OUT_DIR}/env.log"

# Env + platform report (critical for reproducibility).
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

# Cellsim doctor — fail fast if env is broken.
echo "=== cellsim doctor ==="
./scripts/cellsim doctor 2>&1 | tee "${OUT_DIR}/doctor.log" | tail -5
DOCTOR_EXIT=${PIPESTATUS[0]}
if [ $DOCTOR_EXIT -ne 0 ]; then
    echo "error: cellsim doctor failed. Check ${OUT_DIR}/doctor.log" >&2
    echo "Do NOT proceed with the binding run until doctor passes." >&2
    exit 3
fi
echo ""

# Pre-flight validate — refuses to run if the YAML has hard issues
# (typo'd SMILES, missing receptor, charged ligand without correction).
echo "=== fep-binding validate ==="
./scripts/cellsim fep-binding validate \
    benchmarks/fep/binding_streptavidin.yaml 2>&1 | tee -a "${OUT_DIR}/run.log"
VALIDATE_EXIT=${PIPESTATUS[0]}
if [ $VALIDATE_EXIT -ne 0 ]; then
    echo "error: validate FAILED. Fix the YAML before sampling." >&2
    exit 4
fi
echo ""

# Run the binding bench. --resume means a previously-killed run
# (same CSV path) picks up where it left off.
echo "=== binding bench at $(date) ==="
echo "Parameters: 11 windows × 25 000 prod + 2 500 equil steps × 2 legs"
echo "            × 4 streptavidin compounds (biotin / desthiobiotin /"
echo "              2-iminobiotin / biotin methyl ester)"
echo "            FF: amber14-all + tip3pfb + SMIRNOFF ligand"
echo ""

time ./scripts/cellsim fep-binding bench \
    benchmarks/fep/binding_streptavidin.yaml \
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

# Local verdict via fep-report (auto-infers binding gate 2.0 + sign rule).
echo "============================================================"
echo "Local verdict (cellsim fep-report)"
echo "============================================================"
if ./scripts/cellsim fep-report "${CSV}" \
        --yaml benchmarks/fep/binding_streptavidin.yaml \
        --out-dir "${OUT_DIR}" 2>&1 | tee -a "${OUT_DIR}/run.log"; then
    VERDICT_EXIT=0
else
    VERDICT_EXIT=$?
fi
echo ""

# Bundle for transfer (after report so the verdict markdown +
# parity PNG land inside the tarball).
TARBALL="streptavidin_gpu_${STAMP}.tar.gz"
tar czf "${TARBALL}" -C run/fep "streptavidin_${STAMP}/"
echo "Created: ${TARBALL}  ($(du -h ${TARBALL} | cut -f1))"
echo ""
echo "Send this tarball back to Henry. Inside:"
echo "  - report.md           (PASS/FAIL verdict, ranked residuals)"
echo "  - parity.png          (predicted vs expt ΔG_bind, ±2 kcal/mol gate)"
echo "  - table.csv           (normalised, with within-σ + sign-correct cols)"
echo "  - streptavidin_results.csv (raw bench output)"
echo "  - env.log, doctor.log, run.log (provenance)"
echo ""

# Notify on completion.
FINAL_EXIT=${VERDICT_EXIT:-${EXIT_CODE}}
if command -v terminal-notifier >/dev/null 2>&1; then
    STATUS=$([ $FINAL_EXIT -eq 0 ] && echo "PASS" || echo "FAIL")
    terminal-notifier \
        -title "CellSim Streptavidin FEP" \
        -subtitle "${STATUS} (exit ${FINAL_EXIT})" \
        -message "Tarball ready: ${TARBALL}. Send it to Henry." \
        -sound Glass \
        -activate com.apple.Terminal 2>/dev/null || true
fi

exit ${FINAL_EXIT}
