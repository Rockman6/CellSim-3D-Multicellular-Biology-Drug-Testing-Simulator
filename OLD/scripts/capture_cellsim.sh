#!/bin/bash
# Capture the CellSim main window to a PNG (works even with other apps on top)
DIR="$(cd "$(dirname "$0")/.." && pwd)/logs/screenshots"
mkdir -p "$DIR"

# Get CellSim main window ID (the large onscreen one)
WID=$(swift -e '
import CoreGraphics
let ws = CGWindowListCopyWindowInfo([.optionAll], kCGNullWindowID) as! [[String: Any]]
for w in ws {
    if let o = w["kCGWindowOwnerName"] as? String, o.contains("CellSim"),
       let b = w["kCGWindowBounds"] as? [String: Any],
       let h = b["Height"] as? Int, h > 100,
       let onScreen = w["kCGWindowIsOnscreen"] as? Bool, onScreen {
        print(w["kCGWindowNumber"]!)
        break
    }
}
' 2>/dev/null)

if [ -z "$WID" ]; then
    echo "CellSim window not found"
    exit 1
fi

FNAME="${1:-capture_$(date +%Y%m%d_%H%M%S).png}"
OUT="$DIR/$FNAME"
screencapture -l "$WID" "$OUT" 2>/dev/null
if [ -f "$OUT" ]; then
    echo "$OUT"
    # Cleanup: keep only last 50 screenshots
    cd "$DIR" && ls -t *.png 2>/dev/null | tail -n +51 | xargs rm -f 2>/dev/null
else
    echo "Capture failed"
    exit 1
fi
