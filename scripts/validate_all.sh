#!/bin/bash
# Runs one headless simulation per reference (seeded with that ref's init
# cell count) and compares the curve. Skips refs without a cell_count col.
set -e
bio_h="${1:-46}"
printf "%-46s  %6s  %8s  %s\n" "dataset" "init" "final_h" "mean_rel_err"
printf "%-46s  %6s  %8s  %s\n" "-------" "----" "-------" "------------"
for ref in data/reference/growth_curves/ctc_hela_*cellcount*.csv; do
  [ -f "$ref" ] || continue
  name=$(basename "$ref" .csv)
  init=$(awk -F, '!/^#/ && !/^frame/ && NR>=2 { print $3; exit }' "$ref")
  last_h=$(awk -F, '!/^#/ && !/^frame/ { h=$2 } END { print h }' "$ref")
  bh=$(awk "BEGIN{printf \"%d\", ($last_h<$bio_h?$last_h:$bio_h)+0.5}")
  err=$(./build/cellsim_headless "$bh" 60 "$init" "$ref" 2>&1 \
        | grep "mean |relative error|" | awk '{ print $6 }')
  printf "%-46s  %6s  %8s  %s\n" "$name" "$init" "$bh" "$err"
done
