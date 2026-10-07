"""cellsim.cell.trace — record everything a cell does, for watching.

`simulate` keeps only colony-level summaries: how many cells are alive,
mean p53, the fraction in each phase. Underneath, every simulated cell
carries about twenty molecular species that change every few minutes,
and divisions and deaths are discrete events. This records all of it
for a handful of cells so it can be shown — the drug arriving, DNA
damage building, p53 rising, the cyclins cycling, and the death cascade
firing or not.

Nothing here is computed differently from `simulate`; it uses the
engine's own `on_record` hook to copy the state out, so what is shown
is exactly what the engine does.

What the engine does NOT model is worth saying beside any picture of
it: individual molecules (each species is a lumped, normalised level),
the great majority of the proteome, and intracellular metabolism. The
species below are the ones the engine tracks, and only those.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Union

import numpy as np

from cellsim.cell.engine import IX, NAMES, N_BASE, Params, calibrate_cycle_scale, simulate
from cellsim.cell.library import CellLine, Drug, get_drug, get_line

__all__ = ["trace_cells", "Trace", "SPECIES", "PHASES"]

PHASES = ("G1", "S", "G2", "M")

# What each tracked species is, in words a non-specialist can follow.
# The keys are the engine's own names, so nothing is renamed on the way.
SPECIES = {
    "Cin":       "drug inside the cell",
    "D":         "DNA damage",
    "ATM":       "ATM — the damage sensor",
    "p53":       "p53 — the alarm",
    "MDM2_mRNA": "MDM2 message",
    "MDM2":      "MDM2 — switches p53 off",
    "p21":       "p21 — the brake on division",
    "CycD":      "cyclin D — start the cycle",
    "Rb":        "Rb — the restriction gate",
    "E2F":       "E2F — opens the gate",
    "CycE":      "cyclin E — enter S phase",
    "CycA":      "cyclin A — copy DNA",
    "CycB":      "cyclin B — divide",
    "Puma":      "PUMA — the death signal",
    "Bax":       "BAX — punches the mitochondria",
    "MOMP":      "mitochondria opened",
    "CytC":      "cytochrome c released",
    "Smac":      "Smac released",
    "C9":        "caspase-9",
    "C3":        "caspase-3 — the executioner",
    "Exec":      "dismantling the cell",
}


@dataclass
class Trace:
    """Per-cell histories, ready to send to a page as JSON."""
    t_h: list
    line: str
    drug: Optional[str]
    conc_uM: float
    species: dict                       # name -> description
    cells: list = field(default_factory=list)   # one dict per traced cell

    def to_dict(self) -> dict:
        return {"t_h": self.t_h, "line": self.line, "drug": self.drug,
                "conc_uM": self.conc_uM, "species": self.species, "cells": self.cells,
                "schema": "cellsim.cell.trace/v1"}


def _phase(row: np.ndarray, in_m: bool) -> str:
    """The engine's own phase rule, read off a cell's cyclins."""
    if in_m:
        return "M"
    if row[IX["CycA"]] > 0.3 and row[IX["CycB"]] < 0.25:
        return "S"
    if row[IX["CycB"]] >= 0.25:
        return "G2"
    return "G1"


def trace_cells(line: Union[str, CellLine], drug: Optional[Union[str, Drug]] = None,
                conc_uM: float = 0.0, *, hours: float = 72.0, n_cells: int = 6,
                every_h: float = 0.25, exposure_h: Optional[float] = None,
                seed: int = 1, params: Params = Params()) -> Trace:
    """Simulate `n_cells` and record every species of each, every `every_h`.

    Divisions and deaths are events between records, so each is detected
    from what changed: a division as the cell's generation counter going
    up (its cyclins reset to G1 at the same moment), a death as it ceasing
    to be alive. A cell held in mitosis by a spindle poison is flagged
    from the engine's own arrest clock.
    """
    ln = line if isinstance(line, CellLine) else get_line(line)
    dg = None if drug is None else (drug if isinstance(drug, Drug) else get_drug(drug))
    k = calibrate_cycle_scale(ln, params)
    names = list(NAMES)
    times: list = []
    hist = {i: {"state": [], "phase": [], "alive": [], "arrested": []} for i in range(n_cells)}
    events = {i: [] for i in range(n_cells)}
    prev = {"gen": None, "alive": None, "arrest": None}

    def on_record(pop, t_h):
        times.append(round(float(t_h), 4))
        gen = pop.generation[:n_cells].copy()
        alive = pop.alive[:n_cells].copy()
        arrest = (pop.arrest_h[:n_cells] > 0)
        for i in range(n_cells):
            row = pop.Y[i]
            st = {n: float(row[IX[n]]) for n in names if n in IX}
            if pop.Y.shape[1] > N_BASE:
                st["Cin"] = float(row[N_BASE])
            hist[i]["state"].append(st)
            hist[i]["phase"].append(_phase(row, bool(pop.in_M[i])) if alive[i] else "dead")
            hist[i]["alive"].append(bool(alive[i]))
            hist[i]["arrested"].append(bool(arrest[i]) and bool(alive[i]))
            if prev["gen"] is not None:
                if gen[i] > prev["gen"][i]:
                    events[i].append({"t_h": times[-1], "kind": "divided",
                                      "text": f"divided (generation {int(gen[i])})"})
                if prev["alive"][i] and not alive[i]:
                    events[i].append({"t_h": times[-1], "kind": "died",
                                      "text": "died by apoptosis"})
                if arrest[i] and not prev["arrest"][i] and alive[i]:
                    events[i].append({"t_h": times[-1], "kind": "arrested",
                                      "text": "stuck in mitosis"})
        prev["gen"], prev["alive"], prev["arrest"] = gen, alive, arrest

    def dose(t):
        if exposure_h is not None and t >= exposure_h:
            return 0.0
        return conc_uM

    simulate(ln, dg, conc_uM if dg is not None else 0.0, t_end_h=hours, n_cells=n_cells,
             k_cyc=k, p=params, seed=seed, record_every_h=every_h,
             dose_fn=dose if dg is not None else None, on_record=on_record)

    cells = []
    for i in range(n_cells):
        series = {n: [s.get(n, 0.0) for s in hist[i]["state"]] for n in SPECIES}
        cells.append({"index": i, "series": series, "phase": hist[i]["phase"],
                      "alive": hist[i]["alive"], "arrested": hist[i]["arrested"],
                      "events": events[i]})
    return Trace(t_h=times, line=ln.name, drug=None if dg is None else dg.name,
                 conc_uM=float(conc_uM), species=dict(SPECIES), cells=cells)
