#!/usr/bin/env bash
# Regenerate the viewer's demo runs (not committed: they are outputs).
# Usage: bash web/demo/make_demos.sh [out_dir]      (~6 min on a laptop)
set -euo pipefail
OUT="${1:-$(dirname "$0")}"
mkdir -p "$OUT"
dish() { python -m cellsim.cell.dish "$@"; }
dish --line A549 --geometry monolayer --hours 120 --cells 120 --grid 48 --every 4 --seed 1 \
     --out "$OUT/a549-monolayer-growth.jsonl"
dish --line A549 --geometry monolayer --drug cisplatin --schedule "24:20,48:0" --hours 120 \
     --cells 120 --grid 48 --every 4 --seed 1 --out "$OUT/a549-monolayer-cisplatin.jsonl"
dish --line DLD-1 --geometry spheroid --hours 144 --cells 6000 --grid 52 --every 8 --seed 2 \
     --out "$OUT/dld1-spheroid-growth.jsonl"
dish --line DLD-1 --geometry spheroid --drug doxorubicin --schedule "72:2,96:0" --hours 144 \
     --cells 6000 --grid 52 --every 8 --seed 2 --out "$OUT/dld1-spheroid-doxorubicin.jsonl"
python -m cellsim.cell.stream --line A549 --drug paclitaxel --schedule "0:0.06,24:0" --hours 96 \
     --cells 128 --every 2 --seed 1 --out "$OUT/a549-population-paclitaxel.jsonl"
for f in "$OUT"/*.jsonl; do gzip -f "$f"; done
ls -la "$OUT"
