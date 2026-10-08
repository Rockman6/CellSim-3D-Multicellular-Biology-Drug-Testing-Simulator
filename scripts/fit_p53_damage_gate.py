#!/usr/bin/env python3
"""Fit `Params.p53_apoptosis_damage_free` — how much killing p53 does
without DNA-damage signalling.

What it is fitted to. Tovar et al. 2006 (PNAS 103:1888) held p53
wild-type lines in nutlin-3 for 72 h and counted annexin-V-positive
cells: 7.6 ± 1.1 % of A549 and 8.8 ± 2.1 % of HCT116, where the
MDM2-amplified SJSA-1 reached 88 %. In A549 and HCT116, blocking MDM2
arrests the cycle and rarely kills. Before this gate the engine killed
75 % of both at 10 µM by 72 h. The concentration Tovar used was not
re-read from the paper for this fit (the full text was not reachable);
10 µM, the top of the GDSC screen, is assumed and 5 µM is reported
beside it.

How. Nutlin's own fitted constant (its partition coefficient, the one
thing `scripts/validate_gdsc.py` tunes for it) depends on the gate,
because a drug that arrests needs a different exposure to halve a 72-h
cell count than one that kills. So for each candidate gate value the
partition is refitted to GDSC on the usual fit lines (A549, MCF7), and
only then is the 72-h dead fraction read off. The chosen value is the
one whose dead fraction, averaged over A549 and HCT116, is nearest
Tovar's 8.2 %.

What it is NOT fitted to: anything about cisplatin, doxorubicin or the
other DNA-damaging drugs. Those are refitted afterwards by
`validate_gdsc.py` with the gate in place, and the HCT116 cisplatin
time course of Paek et al. 2016 stays a held-out check.

Usage:
    python scripts/fit_p53_damage_gate.py            # ~15 min, 10 cores
    python scripts/fit_p53_damage_gate.py --quick
"""
from __future__ import annotations

import argparse
import dataclasses
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from cellsim.cell.engine import Params, calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import CELL_LINES, get_drug  # noqa: E402
import validate_gdsc as vg  # noqa: E402

TOVAR = {"A549": 0.076, "HCT116": 0.088}      # annexin-V positive after 72 h of nutlin-3
TARGET = sum(TOVAR.values()) / len(TOVAR)
DOSES_UM = (10.0, 5.0)                        # 10 assumed, 5 reported beside it


def dead_fraction(line: str, drug, conc: float, params: Params, *, n_cells: int,
                  seeds=(1, 2, 3)) -> float:
    """Dead weight over all weight at 72 h — what an annexin count of the
    whole culture sees — averaged over seeds."""
    ln = CELL_LINES[line]
    k = calibrate_cycle_scale(ln, params)
    vals = []
    for s in seeds:
        r = simulate(ln, drug, conc, t_end_h=72.0, n_cells=n_cells, k_cyc=k, seed=s, p=params)
        dead, alive = float(r.deaths[-1, 0]), float(r.alive_weight[-1, 0])
        vals.append(dead / (dead + alive))
    return float(np.mean(vals))


def evaluate(gate: float, quick: bool) -> dict:
    params = dataclasses.replace(Params(), p53_apoptosis_damage_free=gate)
    drug = get_drug("nutlin-3a")
    targets = vg.load_targets(vg.REFERENCE_CSV)
    n_cells = 16 if quick else 32
    iters = 8 if quick else 10
    lines = sorted(set(CELL_LINES) & set(vg.FIT_LINES) | {"HCT116"})
    k_cyc = {name: calibrate_cycle_scale(CELL_LINES[name], params) for name in lines}
    gains = []
    for line in vg.FIT_LINES:
        t = targets.get(("nutlin-3a", line))
        if t is None or t.kind != "value":
            continue
        g, _ = vg.fit_gain_for_line(line, drug, t.center, n_cells=n_cells, k_cyc=k_cyc,
                                    iters=iters, params=params)
        if math.isfinite(g):
            gains.append(g)
    if not gains:
        return {"gate": gate, "partition": math.nan}
    partition = float(math.exp(np.mean(np.log(gains))))
    fitted = dataclasses.replace(drug, partition=partition)
    out = {"gate": gate, "partition": partition}
    for line in TOVAR:
        for c in DOSES_UM:
            out[f"{line}@{c:g}"] = dead_fraction(line, fitted, c, params,
                                                 n_cells=64 if quick else 192)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--gates", default="1.0,0.5,0.3,0.2,0.1,0.05,0.02,0.0")
    a = ap.parse_args(argv)
    gates = [float(x) for x in a.gates.split(",")]
    with ProcessPoolExecutor(max_workers=min(len(gates), 8)) as ex:
        rows = list(ex.map(evaluate, gates, [a.quick] * len(gates)))
    print(f"Tovar 2006, 72 h nutlin-3: A549 {TOVAR['A549']:.1%}, HCT116 {TOVAR['HCT116']:.1%} "
          f"(mean {TARGET:.1%})\n")
    hdr = "  ".join(f"{k:>12s}" for k in rows[0] if k not in ("gate", "partition"))
    print(f"{'gate':>6s} {'partition':>10s}  {hdr}")
    best = None
    for r in rows:
        cols = [k for k in r if k not in ("gate", "partition")]
        print(f"{r['gate']:6.2f} {r['partition']:10.4g}  " + "  ".join(f"{r[k]:12.1%}" for k in cols))
        if cols:
            err = abs(np.mean([r[f"{ln}@{DOSES_UM[0]:g}"] for ln in TOVAR]) - TARGET)
            if best is None or err < best[0]:
                best = (err, r)
    if best:
        print(f"\nnearest Tovar at {DOSES_UM[0]:g} uM: gate = {best[1]['gate']}, "
              f"partition = {best[1]['partition']:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
