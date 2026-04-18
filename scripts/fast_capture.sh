#!/bin/bash
# Fast capture using cached window ID
WID=${WID:-167880}
DIR="/Users/henry/CellSim/logs/screenshots"
FNAME="${1:-f_$(date +%s%N | cut -c1-13).png}"
screencapture -l "$WID" "$DIR/$FNAME" 2>/dev/null
echo "$DIR/$FNAME"
