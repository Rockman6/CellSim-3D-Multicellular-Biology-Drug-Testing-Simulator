#!/usr/bin/env python3
"""Cell Tracking Challenge HeLa: counts and cell-cycle times, from ground truth.

Fluo-N2DL-HeLa is HeLa H2B-GFP imaged every 30 min for 46 h (two
sequences, 92 frames each), with expert-curated lineage tracks
(Ulman et al. 2017 Nat Methods 14:1141; Maska et al. 2014 Bioinformatics
30:1609; imaging from Neumann et al. 2010 Nature 464:721). The tracks
(`man_track.txt`: label, first frame, last frame, parent) give two
things a dish can be checked against:

* the number of cells in the field at each frame, and
* the time from a cell's birth to its own division, for every cell whose
  whole cycle falls inside the movie.

Complete cycles are biased short, because a 46 h window cannot contain
a long one, so a Kaplan-Meier estimate that treats cells still undivided
at the end of their track as censored is reported beside them. Cells
also enter and leave the field, which the count trajectory includes and
a closed simulated field does not.

Only these derived numbers are committed
(benchmarks/cell/ctc_hela_reference.json); the dataset itself is
downloaded from the challenge on first use (199 MB, kept in data/ctc/,
which is git-ignored).

Usage:
    python scripts/ctc_reference.py
"""
from __future__ import annotations

import json
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA = REPO_ROOT / "data" / "ctc"
ZIP = DATA / "Fluo-N2DL-HeLa.zip"
URL = "https://data.celltrackingchallenge.net/training-datasets/Fluo-N2DL-HeLa.zip"
OUT_JSON = REPO_ROOT / "benchmarks" / "cell" / "ctc_hela_reference.json"
FRAME_H = 0.5          # 30 min between frames
PIXEL_UM = 0.645       # Fluo-N2DL-HeLa pixel size
FIELD_PX = (1100, 700)


def load_tracks(seq: str) -> np.ndarray:
    member = f"Fluo-N2DL-HeLa/{seq}_GT/TRA/man_track.txt"
    path = DATA / member
    if not path.exists():
        if not ZIP.exists():
            DATA.mkdir(parents=True, exist_ok=True)
            print(f"downloading {URL} (199 MB) ...", flush=True)
            urllib.request.urlretrieve(URL, ZIP)
        with zipfile.ZipFile(ZIP) as z:
            z.extract(member, DATA)
    return np.loadtxt(path, dtype=int, ndmin=2)


def kaplan_meier(durations: np.ndarray, observed: np.ndarray) -> list[tuple[float, float]]:
    """Survival curve of 'not yet divided' against age (hours)."""
    order = np.argsort(durations)
    d, e = durations[order], observed[order]
    S, out = 1.0, []
    for t in np.unique(d):
        at_risk = np.sum(d >= t)
        events = np.sum(e[d == t])
        if events:
            S *= 1.0 - events / at_risk
        out.append((float(t), float(S)))
    return out


def analyse(seq: str) -> dict:
    T = load_tracks(seq)
    label, first, last, parent = T.T
    n_frames = int(last.max()) + 1
    counts = [int(np.sum((first <= f) & (last >= f))) for f in range(n_frames)]
    divided = np.isin(label, parent[parent > 0])
    born = parent > 0
    length_h = (last - first + 1) * FRAME_H
    complete = length_h[born & divided]
    censored = length_h[born & ~divided]
    km = kaplan_meier(np.concatenate([complete, censored]),
                      np.concatenate([np.ones(len(complete)), np.zeros(len(censored))]))
    km_at = {f"{a:g}": next((S for t, S in reversed(km) if t <= a), 1.0)
             for a in (15, 20, 25, 30, 35, 40)}
    # sisters: both daughters of one division completed their cycles
    dur = dict(zip(label.tolist(), length_h.tolist()))
    complete_ids = set(label[born & divided].tolist())
    pairs = []
    for mum in np.unique(parent[parent > 0]):
        kids = label[parent == mum]
        if len(kids) == 2 and all(k in complete_ids for k in kids):
            pairs.append((dur[int(kids[0])], dur[int(kids[1])]))
    sister_r = (float(np.corrcoef(np.array(pairs).T)[0, 1]) if len(pairs) >= 5 else None)
    t = np.arange(n_frames) * FRAME_H
    rate = float(np.polyfit(t, np.log(counts), 1)[0])
    return {
        "frames": n_frames, "frame_interval_h": FRAME_H,
        "field_um": [FIELD_PX[0] * PIXEL_UM, FIELD_PX[1] * PIXEL_UM],
        "counts": counts, "divisions": int(len(np.unique(parent[parent > 0]))),
        "count_doubling_time_h": float(np.log(2) / rate),
        "complete_cycles_h": complete.tolist(),
        "complete_cycle_mean_h": float(complete.mean()),
        "complete_cycle_cv": float(complete.std(ddof=1) / complete.mean()),
        "n_complete": int(len(complete)), "n_censored": int(len(censored)),
        "km_undivided_at_age_h": km_at,
        "sister_pairs": len(pairs), "sister_correlation": sister_r,
    }


def main() -> int:
    report = {
        "dataset": "Cell Tracking Challenge Fluo-N2DL-HeLa (training ground truth, TRA)",
        "cite": ["Ulman V et al. 2017 Nat Methods 14:1141",
                 "Maska M et al. 2014 Bioinformatics 30:1609",
                 "Neumann B et al. 2010 Nature 464:721 (imaging)"],
        "source_url": URL, "sequences": {},
    }
    for seq in ("01", "02"):
        r = analyse(seq)
        report["sequences"][seq] = r
        print(f"sequence {seq}: {r['counts'][0]} -> {r['counts'][-1]} cells in "
              f"{(r['frames'] - 1) * FRAME_H:g} h (count doubling {r['count_doubling_time_h']:.1f} h); "
              f"{r['n_complete']} complete cycles, mean {r['complete_cycle_mean_h']:.1f} h, "
              f"CV {r['complete_cycle_cv']:.2f}; Kaplan-Meier undivided at 20/30 h "
              f"{r['km_undivided_at_age_h']['20']:.2f}/{r['km_undivided_at_age_h']['30']:.2f}; "
              f"sister correlation {r['sister_correlation']} (n={r['sister_pairs']})")
    OUT_JSON.write_text(json.dumps(report, indent=1) + "\n")
    print(f"wrote {OUT_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
