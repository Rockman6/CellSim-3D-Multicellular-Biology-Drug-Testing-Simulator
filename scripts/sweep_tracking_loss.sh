#!/bin/bash
# Sweep TRACKING_LOSS_PROB_PER_BIOSEC to find the value that minimises
# max(seq01_err, seq02_err). Both CTC Fluo-N2DL-HeLa reference curves
# should fall below 7.5 % mean |rel_err|.
set -e
cd "$(dirname "$0")/.."
CONSTS="src/core/Constants.h"

apply_rate() {
    local rate=$1
    sed -i.tmp "s/^constexpr float TRACKING_LOSS_PROB_PER_BIOSEC.*$/constexpr float TRACKING_LOSS_PROB_PER_BIOSEC = ${rate}f;  \/\/ swept/" "$CONSTS"
    rm -f "${CONSTS}.tmp"
}

err_for() {
    local init=$1 ref=$2
    # Average over 5 runs — TRACKING_LOSS is stochastic so single runs
    # vary. Use a short-loop average for stability.
    local sum=0
    for i in 1 2 3 4 5; do
        e=$(./build/cellsim_headless 46 60 "$init" "$ref" 2>&1 \
            | awk -F= '/mean \|relative error\|/ { gsub(/[^0-9.]/,"",$2); print $2; exit }')
        sum=$(python3 -c "print($sum + float('$e'))")
    done
    python3 -c "print(round($sum/5.0, 2))"
}

printf "%-14s  %-8s  %-8s  %s\n" "rate" "seq01" "seq02" "max"
printf "%-14s  %-8s  %-8s  %s\n" "----" "-----" "-----" "---"

best_max=1e9
best_rate=0
for rate in 3.0e-6 6.0e-6 9.0e-6 1.2e-5 1.5e-5 1.8e-5 2.1e-5 2.4e-5; do
    apply_rate "$rate"
    cmake --build build -j > /dev/null 2>&1 || { echo "build fail at $rate"; continue; }
    e1=$(err_for 43 data/reference/growth_curves/ctc_hela_cellcount_seq01.csv)
    e2=$(err_for 125 data/reference/growth_curves/ctc_hela_cellcount_seq02.csv)
    mx=$(python3 -c "print(max(float('$e1'), float('$e2')))")
    printf "%-14s  %-8s  %-8s  %s\n" "$rate" "$e1" "$e2" "$mx"
    if python3 -c "import sys; sys.exit(0 if $mx < $best_max else 1)"; then
        best_max=$mx
        best_rate=$rate
    fi
done

printf "\n[BEST] rate=%s  max(rel_err)=%s%%\n" "$best_rate" "$best_max"
apply_rate "$best_rate"
cmake --build build -j > /dev/null 2>&1
