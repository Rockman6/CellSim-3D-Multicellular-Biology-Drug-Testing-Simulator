#!/bin/bash
# Fine-tune sim rates against the calibration-target dataset
# (CTC Fluo-N2DL-HeLa seq02, 125 cells → 363 over 46 h). Sweeps
# SLOW_DT_SCALE and MECH_P21_COUPLING on a fine grid and picks the
# combination that minimises mean |rel_err| below 7.5 % (target <5 %).
set -e
cd "$(dirname "$0")/.."

CONSTS="src/core/Constants.h"
BACKUP="$CONSTS.bak"
[ -f "$BACKUP" ] || cp "$CONSTS" "$BACKUP"

orig_slow=$(grep "^constexpr float SLOW_DT_SCALE" "$BACKUP" | awk '{print $5}' | tr -d ';f')
orig_medium=$(grep "^constexpr float MEDIUM_DT_SCALE" "$BACKUP" | awk '{print $5}' | tr -d ';f')
ratio_medium=$(python3 -c "print($orig_medium / $orig_slow)")

apply_knobs() {
    local slow=$1 mech=$2
    local medium
    medium=$(python3 -c "print(round($slow * $ratio_medium, 5))")
    sed -i.tmp "s/^constexpr float SLOW_DT_SCALE.*$/constexpr float SLOW_DT_SCALE   = ${slow}f;/" "$CONSTS"
    sed -i.tmp "s/^constexpr float MEDIUM_DT_SCALE.*$/constexpr float MEDIUM_DT_SCALE = ${medium}f;/" "$CONSTS"
    sed -i.tmp "s/^constexpr float MECH_P21_COUPLING.*$/constexpr float MECH_P21_COUPLING    = ${mech}f;/" "$CONSTS"
    rm -f "${CONSTS}.tmp"
}

err_seq02() {
    ./build/cellsim_headless 46 60 125 \
        data/reference/growth_curves/ctc_hela_cellcount_seq02.csv 2>&1 \
        | awk -F= '/mean \|relative error\|/ { gsub(/[^0-9.]/,"",$2); print $2; exit }'
}

best_err=1e9
best_slow=$orig_slow
best_mech=0.018

printf "\n%-8s %-8s  %s\n" "SLOW" "MECH" "seq02_err%"
for slow in 0.055 0.053 0.052 0.051 0.050 0.049 0.048; do
    for mech in 0.018 0.010 0.005 0.002 0.001 0.000; do
        apply_knobs "$slow" "$mech"
        cmake --build build -j 2>&1 > /dev/null 2>&1 || { echo "  build FAIL"; continue; }
        err=$(err_seq02)
        printf "%-8s %-8s  %s\n" "$slow" "$mech" "$err"
        is_better=$(python3 -c "print('yes' if float('$err') < $best_err else 'no')")
        if [ "$is_better" = "yes" ]; then
            best_err=$err
            best_slow=$slow
            best_mech=$mech
        fi
    done
done

printf "\n[BEST] SLOW=%s MECH=%s seq02_err=%s%%\n" "$best_slow" "$best_mech" "$best_err"
apply_knobs "$best_slow" "$best_mech"
cmake --build build -j > /dev/null 2>&1
echo
echo "Verifying across all datasets..."
./scripts/validate_unseen.sh 2>&1 | grep -E "==|rel_err\s*=|doubling_time sim|doubling_time ref|\+[0-9]"
