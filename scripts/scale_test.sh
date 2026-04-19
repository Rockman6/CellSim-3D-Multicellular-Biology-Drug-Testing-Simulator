#!/bin/bash
# Scale test harness for CellSim headless. Runs a matrix of (init_cells,
# bio_hours) combinations and reports wall-time / final-cells / deaths
# so you can pick a run size for the 64GB machine.
#
# Usage: ./scripts/scale_test.sh
set -e
cd "$(dirname "$0")/.."

printf "%-6s %-6s %-6s %-6s %-6s %-8s\n" "init" "bio_h" "wall_s" "final" "divs" "deaths"
printf "%-6s %-6s %-6s %-6s %-6s %-8s\n" "----" "-----" "------" "-----" "----" "------"
for cfg in "125 46" "250 46" "500 48" "500 96" "800 48" "800 120" "1200 48"; do
  set -- $cfg
  init=$1; bh=$2
  t0=$(python3 -c "import time; print(time.time())")
  out=$(./build/cellsim_headless "$bh" 60 "$init" /dev/null 2>&1)
  t1=$(python3 -c "import time; print(time.time())")
  wall=$(awk "BEGIN{printf \"%.1f\", $t1-$t0}")
  line=$(echo "$out" | grep "final cells")
  final=$(echo "$line"  | sed -n 's/.*final cells=\([0-9]*\).*/\1/p')
  divs=$(echo "$line"   | sed -n 's/.*divisions=\([0-9]*\).*/\1/p')
  deaths=$(echo "$line" | sed -n 's/.*deaths=\([0-9]*\).*/\1/p')
  printf "%-6s %-6s %-6s %-6s %-6s %-8s\n" "$init" "$bh" "$wall" "$final" "$divs" "$deaths"
done
