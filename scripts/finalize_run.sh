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
# about M5 Max when the run was on CPU laptop. Detect rules:
#
# - HW string: derive a short friendly label from the `machine:` line
#   ("Apple MacBook Pro" / "Apple Mac" / "Linux node"), NOT the full
#   uname spew which truncates ugly in the email.
# - Platform: scan the run.log for "FEP sampling platform: <NAME>"
#   which is what the pipeline ACTUALLY selected. Fall back to the
#   env.log "Available platforms" block (selecting the fastest-listed
#   real backend), never the "Pipeline will prefer" preference line
#   (which lists Metal / OpenCL even on builds that don't have them).
HW_DEFAULT="unspecified host"
PLAT_DEFAULT="unspecified"
if [ -f "${RUN_DIR}/env.log" ]; then
    MACHINE_LINE="$(grep -m1 '^  machine:' "${RUN_DIR}/env.log" || true)"
    if echo "${MACHINE_LINE}" | grep -q "MacBook Pro"; then
        HW_DEFAULT="Apple MacBook Pro"
    elif echo "${MACHINE_LINE}" | grep -qi "Darwin.*arm64"; then
        HW_DEFAULT="Apple Silicon Mac"
    elif echo "${MACHINE_LINE}" | grep -qi "Linux"; then
        HW_DEFAULT="Linux node"
    fi
fi
if [ -f "${RUN_DIR}/run.log" ] && grep -q "FEP sampling platform:" "${RUN_DIR}/run.log"; then
    PLAT_DEFAULT="$(grep -m1 'FEP sampling platform:' "${RUN_DIR}/run.log" | sed -E 's|.*FEP sampling platform:[[:space:]]+([A-Za-z]+).*|\1|')"
elif [ -f "${RUN_DIR}/env.log" ]; then
    # Parse the "Available platforms" block, take the highest-speed
    # non-Reference entry (Reference is OpenMM's "always present
    # baseline" and is never selected by the pipeline unless it's
    # the only option). Use sed instead of gawk-only match() so this
    # works under macOS / BSD awk.
    AVAIL=$(awk '/Available OpenMM platforms/,/Pipeline will prefer/' "${RUN_DIR}/env.log" \
            | sed -nE 's/.*[[:space:]]+([A-Za-z]+) \(speed ([0-9.]+)\).*/\2 \1/p' \
            | grep -vE " Reference$" \
            | sort -nr | head -1 | awk '{print $2}')
    if [ -z "${AVAIL}" ]; then
        AVAIL="Reference"   # truly no other backend available
    fi
    PLAT_DEFAULT="${AVAIL}"
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
