#!/bin/bash
# Derive frame-by-frame cell count from Cell Tracking Challenge man_track.txt.
# Usage: ./derive_ctc_counts.sh <man_track.txt> <minutes_per_frame> <out_csv> <comment>
set -e
in="$1"; mpf="$2"; out="$3"; comment="$4"
nframes=$(awk '{ if ($3 > m) m = $3 } END { print m }' "$in")
total=$(wc -l < "$in")
divs=$(awk '$4 != 0' "$in" | wc -l | tr -d ' ')
{
  echo "# ${comment}"
  echo "# Imaging: ${mpf} min/frame, $((nframes+1)) frames = $(awk "BEGIN{print ($nframes+1)*$mpf/60}") hours"
  echo "# Total cells tracked: ${total}"
  echo "# Division events: ${divs}"
  echo "frame,hours,cell_count,divisions_this_frame"
  for (( f=0; f<=$nframes; f++ )); do
    count=$(awk -v f="$f" '$2 <= f && $3 >= f' "$in" | wc -l | tr -d ' ')
    # Births at frame f = tracks whose first_frame == f AND parent != 0
    births=$(awk -v f="$f" '$2 == f && $4 != 0' "$in" | wc -l | tr -d ' ')
    # Divisions this frame = new pair of daughters / 2
    divs_frame=$(( births / 2 ))
    h=$(awk "BEGIN{print $f*$mpf/60}")
    echo "$f,$h,$count,$divs_frame"
  done
} > "$out"
echo "wrote $out"
