#!/bin/bash
# Headless validation against datasets we DIDN'T calibrate against.
# The only fair cross-dataset metric is log-phase doubling time — our
# 550-µm simulated dish saturates at ~500 cells, while lab flasks/wells
# support 10³-10⁶ cells. So we compare the growth RATE (doublings/hour)
# during the linear-on-log phase, not absolute counts.
#
# Produces per-dataset sim vs reference doubling-time comparison.
set -e
cd "$(dirname "$0")/.."

run_one() {
    local label="$1"; local init="$2"; local bh="$3"
    local ref="$4";   local schema="$5"
    printf "\n== %s ==\n" "$label"
    printf "   init=%d, bio_h=%d, schema=%s\n" "$init" "$bh" "$schema"
    ./build/cellsim_headless "$bh" 60 "$init" /dev/null 2>&1 | tail -1
    python3 scripts/compare_against.py logs/headless_doubling.csv "$ref" "$schema" \
        2>&1 | grep -E "doubling_time|rel_err|Δ" || true
}

run_one "CTC Fluo-N2DL seq01 (sparse, slow)" 43 46  \
    data/reference/growth_curves/ctc_hela_cellcount_seq01.csv ctc_cellcount

run_one "CTC Fluo-N2DL seq02 (dense, fast) [CALIBRATION TARGET]" 125 46  \
    data/reference/growth_curves/ctc_hela_cellcount_seq02.csv ctc_cellcount

run_one "CTC DIC-C2DH seq01 (very sparse, DIC)" 10 14  \
    data/reference/growth_curves/ctc_hela_dic_cellcount_seq01.csv ctc_cellcount

run_one "CTC DIC-C2DH seq02 (very sparse, DIC)" 10 14  \
    data/reference/growth_curves/ctc_hela_dic_cellcount_seq02.csv ctc_cellcount

run_one "HeLa viability MTT (Cong 2013, 2000 seed/96-well)" 200 96 \
    data/reference/growth_curves/hela_viability_timecourse.csv viability

run_one "HeLa growth curve (ATCC / Freshney / Eagle 1955, T75)" 125 48 \
    data/reference/growth_curves/hela_growth_curve.csv growth_curve

run_one "HeLa clonal expansion (Cellosaurus / Chao 2019)" 125 120 \
    data/reference/growth_curves/hela_single_cell_origin.csv single_cell
