"""cellsim.cell.stream — the engine's output contract for any interface.

A UI never runs biology; it renders this stream. One run produces a
header record followed by one snapshot per record tick, each a plain
JSON-serialisable dict, so the same output feeds a web viewer, a
notebook, or a file on disk (JSON Lines via `write_jsonl`).

Schema ``cellsim.cell.stream/v1``:

    header   {"schema", "line", "drug", "schedule", "t_end_h", "dt_h",
              "record_every_h", "n_cells", "seed", "params", "doubling_time_h",
              "p53_functional", "fields"}
    snapshot {"t_h", "drug_uM",
              "population": {"alive_weight", "dead_weight", "fold_change",
                             "phase_frac": {"G1","S","G2","M"},
                             "mean": {"p53","Puma","C3","D","Cin_uM"}},
              "cells": {"alive","weight","phase","generation",
                        "p53","Puma","C3","D","Cin_uM","in_mitosis"}}

`cells` holds the representative cells (one entry each, in a fixed
order for the whole run); `weight` is how many real cells each one
stands for, so population totals are weight-sums over alive cells.
Phases: 0 = G1, 1 = S, 2 = G2, 3 = M.
"""
from __future__ import annotations

import dataclasses
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from cellsim.cell.engine import IX, Params, Population, _phase, calibrate_cycle_scale, simulate
from cellsim.cell.library import CellLine, Drug

SCHEMA = "cellsim.cell.stream/v1"
CELL_FIELDS = ("alive", "weight", "phase", "generation", "p53", "Puma", "C3", "D",
               "Cin_uM", "in_mitosis")


@dataclass(frozen=True)
class Schedule:
    """Piecewise-constant extracellular concentration.

    ``segments`` is a sequence of (start_h, concentration_uM), sorted by
    start time; the concentration holds until the next segment starts.
    A wash-out is a segment at 0 µM. Before the first segment the
    concentration is 0.
    """
    segments: tuple[tuple[float, float], ...]

    def __post_init__(self):
        starts = [s for s, _ in self.segments]
        if starts != sorted(starts):
            raise ValueError("schedule segments must be sorted by start time")
        if any(c < 0 for _, c in self.segments):
            raise ValueError("concentrations must be non-negative")

    @classmethod
    def constant(cls, conc_uM: float) -> "Schedule":
        return cls(((0.0, float(conc_uM)),))

    @classmethod
    def parse(cls, text: str) -> "Schedule":
        """'0:10,24:0' -> 10 µM from 0 h, washed out at 24 h."""
        segs = []
        for part in text.split(","):
            part = part.strip()
            if not part:
                continue
            t, c = part.split(":")
            segs.append((float(t), float(c)))
        return cls(tuple(segs))

    def __call__(self, t_h: float) -> float:
        conc = 0.0
        for start, c in self.segments:
            if t_h >= start - 1e-12:
                conc = c
            else:
                break
        return conc


def _snapshot(pop: Population, t_h: float, drug_uM: float, n0: float) -> dict:
    Y = pop.Y
    alive = pop.alive
    aw = pop.weight * alive
    ph = np.where(pop.in_M, 3, _phase(Y))
    tot = max(aw.sum(), 1e-12)

    def mean(col: str) -> float:
        return float(np.average(Y[:, IX[col]], weights=aw)) if aw.sum() > 0 else float("nan")

    return {
        "t_h": round(float(t_h), 6),
        "drug_uM": float(drug_uM),
        "population": {
            "alive_weight": float(aw.sum()),
            "dead_weight": float((pop.weight * ~alive).sum()),
            "fold_change": float(aw.sum() / n0),
            "phase_frac": {name: float(aw[ph == q].sum() / tot)
                           for q, name in enumerate(("G1", "S", "G2", "M"))},
            "mean": {"p53": mean("p53"), "Puma": mean("Puma"), "C3": mean("C3"),
                     "D": mean("D"), "Cin_uM": mean("Cin")},
        },
        "cells": {
            "alive": alive.astype(int).tolist(),
            "weight": pop.weight.tolist(),
            "phase": ph.astype(int).tolist(),
            "generation": pop.generation.astype(int).tolist(),
            "p53": np.round(Y[:, IX["p53"]], 5).tolist(),
            "Puma": np.round(Y[:, IX["Puma"]], 5).tolist(),
            "C3": np.round(Y[:, IX["C3"]], 5).tolist(),
            "D": np.round(Y[:, IX["D"]], 5).tolist(),
            "Cin_uM": np.round(Y[:, IX["Cin"]], 8).tolist(),
            "in_mitosis": pop.in_M.astype(int).tolist(),
        },
    }


def run_stream(line: CellLine, drug: Optional[Drug], schedule: Schedule, *,
               t_end_h: float = 72.0, n_cells: int = 64, dt_h: float = 0.02,
               record_every_h: float = 1.0, seed: int = 1,
               p: Params = Params()) -> list[dict]:
    """Run one population under a dosing schedule; return [header, *snapshots]."""
    header = {
        "schema": SCHEMA,
        "line": line.name, "drug": None if drug is None else drug.name,
        "schedule": [list(s) for s in schedule.segments],
        "t_end_h": t_end_h, "dt_h": dt_h, "record_every_h": record_every_h,
        "n_cells": n_cells, "seed": seed, "params": dataclasses.asdict(p),
        "doubling_time_h": line.doubling_time_h, "p53_functional": line.p53_functional,
        "fields": list(CELL_FIELDS),
    }
    snaps: list[dict] = []
    n0 = float(n_cells)

    def on_record(pop: Population, t_h: float) -> None:
        snaps.append(_snapshot(pop, t_h, schedule(t_h), n0))

    simulate(line, drug, 0.0, t_end_h=t_end_h, n_cells=n_cells, dt_h=dt_h, seed=seed, p=p,
             k_cyc=calibrate_cycle_scale(line, p), record_every_h=record_every_h,
             dose_fn=schedule, on_record=on_record)
    return [header, *snaps]


def write_jsonl(records: Sequence[dict], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for r in records:
            f.write(json.dumps(r, separators=(",", ":")) + "\n")


def main(argv: Optional[list[str]] = None) -> int:
    import argparse
    from cellsim.cell.library import get_drug, get_line

    ap = argparse.ArgumentParser(prog="cellsim cell-stream", description=(
        "Run one cell population under a dosing schedule and write the per-tick "
        "state stream (JSON Lines, schema " + SCHEMA + ")."))
    ap.add_argument("--line", default="A549")
    ap.add_argument("--drug", default="cisplatin", help="drug name, or 'none' for untreated")
    ap.add_argument("--schedule", default="0:10",
                    help="'start_h:conc_uM,...' piecewise constant, e.g. '0:10,24:0' (wash-out at 24 h)")
    ap.add_argument("--hours", type=float, default=72.0)
    ap.add_argument("--cells", type=int, default=64)
    ap.add_argument("--every", type=float, default=1.0, help="record interval, hours")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True, help="output .jsonl")
    a = ap.parse_args(argv)
    drug = None if a.drug.lower() == "none" else get_drug(a.drug)
    recs = run_stream(get_line(a.line), drug, Schedule.parse(a.schedule), t_end_h=a.hours,
                      n_cells=a.cells, record_every_h=a.every, seed=a.seed)
    write_jsonl(recs, a.out)
    last = recs[-1]["population"]
    print(f"wrote {len(recs) - 1} snapshots to {a.out}  "
          f"(final fold change x{last['fold_change']:.2f}, dead weight {last['dead_weight']:.1f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
