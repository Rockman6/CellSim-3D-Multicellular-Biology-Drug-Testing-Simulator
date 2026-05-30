#!/bin/bash
# finalize_run.sh — post-run convenience for a FreeSolv / binding bench.
#
# Usage:
#   bash scripts/finalize_run.sh <run-dir>
#   bash scripts/finalize_run.sh <run-dir> --prof-name "Dr. Lee"
#
# What it does:
#   1. cellsim fep-report <run-dir> --out-dir <run-dir>  → report.md + table.csv + parity.png
#   2. python scripts/csv_tldr.py <run-dir>/*.csv         → one-line status to stdout
#   3. python scripts/fill_prof_email.py <run-dir> ...    → email draft → <run-dir>/prof_email_draft.txt
#   4. tar -czf <run-dir>.tar.gz <run-dir>                → send-ready tarball alongside the dir
#
# All four are idempotent. Re-running after a partial run that closed
# with --resume regenerates everything from the current CSV.

set -euo pipefail

if [ -z "${1:-}" ]; then
    echo "usage: $0 <run-dir> [--prof-name NAME] [--hardware H] [--platform P]" >&2
    exit 1
fi

RUN_DIR="${1}"
shift

if [ ! -d "${RUN_DIR}" ]; then
    echo "error: ${RUN_DIR} is not a directory" >&2
    exit 1
fi

CSV_PATH=$(ls "${RUN_DIR}"/*results.csv 2>/dev/null | head -1)
if [ -z "${CSV_PATH}" ]; then
    echo "error: no *results.csv under ${RUN_DIR}" >&2
    exit 1
fi

# Auto-detect hardware from env.log so the prof email isn't lying
# about M5 Max when the run was on CPU laptop.
HW_DEFAULT="rented GPU"
PLAT_DEFAULT="CUDA"
if [ -f "${RUN_DIR}/env.log" ]; then
    if grep -q "MacBook Pro\|M-series\|MacBook" "${RUN_DIR}/env.log"; then
        HW_DEFAULT="$(grep -m1 '^  machine:' "${RUN_DIR}/env.log" | sed 's|^  machine: *||' | head -1)"
        HW_DEFAULT="${HW_DEFAULT:0:80}"   # truncate uname spew
    fi
    if grep -q "Metal" "${RUN_DIR}/env.log"; then
        PLAT_DEFAULT="Metal"
    elif grep -q "OpenCL" "${RUN_DIR}/env.log"; then
        PLAT_DEFAULT="OpenCL"
    elif grep -q "CPU (speed" "${RUN_DIR}/env.log"; then
        PLAT_DEFAULT="CPU"
    fi
fi

PROF_NAME="[Prof]"
HARDWARE="${HW_DEFAULT}"
PLATFORM="${PLAT_DEFAULT}"

# Pass-through args
while [ $# -gt 0 ]; do
    case "$1" in
        --prof-name) PROF_NAME="$2"; shift 2 ;;
        --hardware)  HARDWARE="$2"; shift 2 ;;
        --platform)  PLATFORM="$2"; shift 2 ;;
        *) echo "warning: unknown arg $1, ignoring" >&2; shift ;;
    esac
done

echo "[finalize_run] run dir : ${RUN_DIR}"
echo "[finalize_run] csv     : ${CSV_PATH}"
echo "[finalize_run] prof    : ${PROF_NAME}"
echo "[finalize_run] hw      : ${HARDWARE}"
echo "[finalize_run] platform: ${PLATFORM}"
echo ""

# 1. fep-report
echo "=== cellsim fep-report ==="
./scripts/cellsim fep-report "${RUN_DIR}" --out-dir "${RUN_DIR}" 2>&1 | tail -20
echo ""

# 2. csv_tldr — one-line status
echo "=== csv_tldr ==="
python scripts/csv_tldr.py "${CSV_PATH}" || true
echo ""

# 3. fill_prof_email — tolerate exit 2 (warn = missing fields, draft
# still written; biologist inspects manually before sending).
echo "=== fill_prof_email ==="
set +e
python scripts/fill_prof_email.py "${RUN_DIR}" \
    --prof-name "${PROF_NAME}" \
    --hardware "${HARDWARE}" \
    --platform "${PLATFORM}" \
    --out "${RUN_DIR}/prof_email_draft.txt" 2>&1 | tail -3
FILL_RC=$?
set -e
if [ $FILL_RC -eq 0 ]; then
    echo "wrote ${RUN_DIR}/prof_email_draft.txt"
elif [ $FILL_RC -eq 2 ]; then
    echo "wrote ${RUN_DIR}/prof_email_draft.txt  (with <missing> warnings; inspect before sending)"
else
    echo "fill_prof_email FAILED (exit $FILL_RC) — see message above" >&2
fi
echo ""

# 4. tarball alongside the dir
TARBALL="${RUN_DIR}.tar.gz"
echo "=== tarball ==="
tar -czf "${TARBALL}" -C "$(dirname "${RUN_DIR}")" "$(basename "${RUN_DIR}")"
SIZE_MB=$(du -m "${TARBALL}" | awk '{print $1}')
echo "wrote ${TARBALL}  (${SIZE_MB} MB)"
echo ""

echo "=== DONE ==="
echo "Verdict:  $(cat "${RUN_DIR}/report.md" | head -1)"
echo "Tarball:  ${TARBALL}"
echo "Email:    ${RUN_DIR}/prof_email_draft.txt"
