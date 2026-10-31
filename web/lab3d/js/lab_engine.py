"""The lab's bridge to the Cell Sim engine.

Runs `cellsim.cell.dish.run_dish_stream` for one culture dish: the cell line,
the drugs dropped into it and when, and how long it has been incubated. The
dish is re-simulated from time zero with the same seed every time, so adding
a drug later replays the same cells up to that moment. Returns the population
readout for every record and the full cell snapshot at the dish's current time.
"""
import json

from cellsim.cell.dish import run_dish_stream
from cellsim.cell.library import get_drug, get_line


def lab_run_dish(cfg_json: str) -> str:
    cfg = json.loads(cfg_json)
    names = list(cfg.get("drugs") or [])
    events = sorted(cfg.get("events") or [], key=lambda e: e["t"])

    def dose(t_h):
        cur = {}
        for e in events:
            if t_h + 1e-9 >= e["t"]:
                cur = e["conc"]
        v = [float(cur.get(n, 0.0)) for n in names]
        return v if len(v) > 1 else v[0]

    drugs = [get_drug(n) for n in names]
    drug = None if not drugs else (drugs[0] if len(drugs) == 1 else tuple(drugs))
    t_now = float(cfg["t"])
    every = float(cfg.get("every", 2.0))
    stream = run_dish_stream(get_line(cfg["line"]), drug, dose if names else None, geometry="monolayer",
                             t_end_h=max(t_now, every), n_seed=int(cfg["n_seed"]),
                             grid_sites=int(cfg["grid"]), record_every_h=every, seed=int(cfg["seed"]))
    header, snaps = stream[0], stream[1:]
    snap = min(snaps, key=lambda s: abs(s["t_h"] - t_now))
    series = [{"t_h": s["t_h"], "drug_uM": s.get("drug_uM"), **s["population"]} for s in snaps if s["t_h"] <= t_now + 1e-6]
    keep = {k: header[k] for k in ("schema", "line", "drug", "grid_shape", "spacing_um", "doubling_time_h", "p53_functional", "seed")}
    return json.dumps({"header": keep, "snap": snap, "series": series})
