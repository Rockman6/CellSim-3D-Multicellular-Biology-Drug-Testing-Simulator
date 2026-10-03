#!/usr/bin/env python3
"""The engine -> UI contract (cellsim.cell.stream, schema v1)."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.stream import CELL_FIELDS, SCHEMA, Schedule, run_stream, write_jsonl  # noqa: E402

A549 = get_line("A549")


def test_schedule_is_piecewise_constant_and_parses():
    s = Schedule.parse("0:10, 24:0")
    assert s(0.0) == 10 and s(23.99) == 10 and s(24.0) == 0 and s(70) == 0
    assert Schedule.parse("6:1")(3.0) == 0.0, "zero before the first segment"
    try:
        Schedule(((5.0, 1.0), (1.0, 2.0)))
    except ValueError:
        pass
    else:
        raise AssertionError("unsorted segments must be rejected")


def test_stream_has_header_and_one_snapshot_per_tick():
    recs = run_stream(A549, None, Schedule.constant(0.0), t_end_h=12, n_cells=16,
                      record_every_h=2.0)
    header, snaps = recs[0], recs[1:]
    assert header["schema"] == SCHEMA and header["line"] == "A549"
    assert header["fields"] == list(CELL_FIELDS)
    assert len(snaps) == 12 // 2 + 1, len(snaps)
    assert [s["t_h"] for s in snaps] == [0, 2, 4, 6, 8, 10, 12]
    for s in snaps:
        for f in CELL_FIELDS:
            assert len(s["cells"][f]) == 16, f
        pf = s["population"]["phase_frac"]
        assert abs(sum(pf.values()) - 1.0) < 1e-6


def test_stream_round_trips_through_json_lines():
    recs = run_stream(A549, get_drug("cisplatin"), Schedule.constant(5.0), t_end_h=6,
                      n_cells=8, record_every_h=3.0)
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "run.jsonl"
        write_jsonl(recs, path)
        back = [json.loads(line) for line in path.read_text().splitlines()]
    assert back == json.loads(json.dumps(recs)), "stream must be plain JSON"


def test_recording_does_not_perturb_the_simulation():
    k = calibrate_cycle_scale(A549)
    a = simulate(A549, get_drug("cisplatin"), 20.0, t_end_h=24, n_cells=32, k_cyc=k)
    seen = []
    b = simulate(A549, get_drug("cisplatin"), 20.0, t_end_h=24, n_cells=32, k_cyc=k,
                 on_record=lambda pop, t: seen.append(t))
    assert np.allclose(a.alive_weight, b.alive_weight) and len(seen) == 25


def test_wash_out_rescues_compared_with_continuous_exposure():
    """Drug removed at 6 h: the colony must end larger than under the same
    dose held for 72 h (damage is repaired once uptake stops). Potency is
    pinned to the GDSC-fitted order of magnitude (IC50 ~10 µM on A549) so
    the test does not depend on library defaults."""
    import dataclasses
    drug = dataclasses.replace(get_drug("cisplatin"), k_damage_per_uM_h=0.001)
    cont = run_stream(A549, drug, Schedule.constant(30.0), t_end_h=72, n_cells=32)
    pulse = run_stream(A549, drug, Schedule.parse("0:30,6:0"), t_end_h=72, n_cells=32)
    fc_cont = cont[-1]["population"]["fold_change"]
    fc_pulse = pulse[-1]["population"]["fold_change"]
    assert pulse[-1]["drug_uM"] == 0.0 and cont[-1]["drug_uM"] == 30.0
    assert fc_pulse > fc_cont, f"wash-out x{fc_pulse:.2f} not above continuous x{fc_cont:.2f}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} stream gates pass")
    raise SystemExit(1 if failed else 0)
