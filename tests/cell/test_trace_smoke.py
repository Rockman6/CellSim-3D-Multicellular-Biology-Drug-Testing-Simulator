#!/usr/bin/env python3
"""Gates for cellsim.cell.trace, which records what each cell does so the
"Inside a cell" page can show it.

The recorder adds nothing to the engine — it copies state out through
simulate's own on_record hook — so these check that what it reports is
the engine's behaviour and not an artefact of reading it: untreated cells
divide and live, a lethal dose kills, events come in time order, and the
species the page draws are all present.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _tiers import run  # noqa: E402
from cellsim.cell.trace import SPECIES, trace_cells  # noqa: E402


def test_untreated_cells_divide_and_none_die():
    tr = trace_cells("A549", None, 0.0, hours=72, n_cells=4)
    for c in tr.cells:
        kinds = [e["kind"] for e in c["events"]]
        assert "died" not in kinds, f"cell {c['index']} died with no drug: {kinds}"
        assert kinds.count("divided") >= 2, (
            f"A549 divides about every 22 h, so 72 h should hold two or more: {kinds}")
        assert all(c["alive"]), "an untreated cell went dead"


def test_a_lethal_dose_kills_and_the_death_cascade_fires():
    tr = trace_cells("A549", "cisplatin", 40.0, hours=72, n_cells=4)
    for c in tr.cells:
        died = [e for e in c["events"] if e["kind"] == "died"]
        assert died, f"cell {c['index']} survived 40 uM cisplatin for 72 h"
        assert max(c["series"]["C3"]) > 0.5, "death should come with caspase-3 active"
        assert not c["alive"][-1]


def test_a_spindle_poison_traps_cells_in_mitosis():
    tr = trace_cells("A549", "paclitaxel", 0.1, hours=48, n_cells=4)
    assert any(any(c["arrested"]) for c in tr.cells), "paclitaxel should arrest some cells"
    assert any(e["kind"] == "arrested" for c in tr.cells for e in c["events"])


def test_events_come_in_time_order_and_inside_the_run():
    tr = trace_cells("A549", "cisplatin", 20.0, hours=72, n_cells=4)
    end = tr.t_h[-1]
    for c in tr.cells:
        ts = [e["t_h"] for e in c["events"]]
        assert ts == sorted(ts), ts
        assert all(0 <= t <= end + 1e-9 for t in ts), ts
        deaths = [e for e in c["events"] if e["kind"] == "died"]
        assert len(deaths) <= 1, "a cell can die only once"
        if deaths:
            assert all(e["t_h"] <= deaths[0]["t_h"] for e in c["events"]), (
                "nothing can happen to a cell after it has died")


def test_every_species_the_page_draws_is_recorded():
    tr = trace_cells("A549", "cisplatin", 10.0, hours=12, n_cells=2)
    n = len(tr.t_h)
    for c in tr.cells:
        for name in SPECIES:
            assert name in c["series"], f"{name} missing"
            assert len(c["series"][name]) == n, f"{name} has the wrong length"
        assert len(c["phase"]) == n and len(c["alive"]) == n


if __name__ == "__main__":
    raise SystemExit(run(globals(), "trace"))
