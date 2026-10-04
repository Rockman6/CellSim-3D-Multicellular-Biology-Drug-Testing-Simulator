#!/usr/bin/env python3
"""The cell simulator must install and run from pip alone.

A cell biologist should not have to install RDKit, OpenMM, AutoDock Vina
and AmberTools to simulate cells. Everything under `cellsim.cell` needs
only NumPy and SciPy, and this test holds that true: it builds a bare
virtual environment, installs the project into it with pip, and runs a
dose-response, a dish and the two CLIs there — with no conda environment
and no repository on `sys.path`.

It is the gate behind the README's "pip install cellsim" instruction. If
someone imports RDKit from `cellsim/cell/` or from the CLI's import
path, this fails.

Run:
    python tests/test_pip_install.py
    python tests/test_pip_install.py --keep   # leave the venv for inspection
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import venv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

PROBE = r'''
import json, sys, numpy as np
out = {"python": sys.version.split()[0]}
# nothing from the repository: the installed package only
assert not any("Cell Sim" in p for p in sys.path if p), sys.path
from cellsim.cell.engine import Params, dose_response, ic50, ic50_from_curve, simulate
from cellsim.cell.library import get_drug, get_line
from cellsim.cell.dish import DishParams, run_dish
from cellsim.cell.stream import Schedule, run_stream
from cellsim.cell.selection import pulsed, summarize, survival_map
line, drug = get_line("A549"), get_drug("cisplatin")
viab, _ = dose_response(line, drug, np.array([3.0, 10.0, 30.0]), n_cells_per_conc=16, t_end_h=24)
out["viability"] = [round(float(v), 3) for v in viab[1:]]
out["ic50_uM"] = round(float(ic50(line, drug, guess_uM=10.0, n_cells_per_conc=12, t_end_h=24)), 3)
recs = run_stream(line, drug, Schedule.parse("0:10,12:0"), t_end_h=24, n_cells=8, record_every_h=12)
out["stream_snapshots"] = len(recs) - 1
res = run_dish(line, t_end_h=12, n_seed=20, grid_sites=16, record_every_h=12)
out["dish_live"] = res.record.n_live
res = run_dish(get_line("DLD-1"), t_end_h=6, n_seed=400, grid_sites=20, record_every_h=6,
               dp=DishParams(geometry="spheroid"))
out["spheroid_min_o2"] = round(res.record.min_o2_mmHg[-1], 2)
# the modules that need the chemistry stack must NOT be reachable, or the
# dependency claim is wrong somewhere
import importlib
try:
    importlib.import_module("rdkit")
    out["rdkit_present"] = True
except ImportError:
    out["rdkit_present"] = False
print("PROBE" + json.dumps(out))
'''


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def main(keep: bool = False) -> int:
    tmp = Path(tempfile.mkdtemp(prefix="cellsim-pip-"))
    env_dir = tmp / "venv"
    print(f"[pip] building a bare virtual environment in {env_dir}", flush=True)
    venv.EnvBuilder(with_pip=True, clear=True).create(env_dir)
    py = env_dir / ("Scripts" if os.name == "nt" else "bin") / "python"
    cellsim_exe = env_dir / ("Scripts" if os.name == "nt" else "bin") / "cellsim"

    print("[pip] pip install <repo>", flush=True)
    r = run([str(py), "-m", "pip", "install", "--quiet", str(REPO_ROOT)])
    if r.returncode != 0:
        print(r.stdout[-3000:] + r.stderr[-3000:])
        print("FAIL: pip install failed")
        return 1

    # The probe runs from a directory that is not the repository, so an
    # accidental relative import cannot rescue a missing dependency.
    r = run([str(py), "-c", PROBE], cwd=str(tmp))
    if r.returncode != 0:
        print(r.stdout[-3000:] + r.stderr[-3000:])
        print("FAIL: the installed package could not run the cell simulator")
        return 1
    line = next(ln for ln in r.stdout.splitlines() if ln.startswith("PROBE"))
    out = json.loads(line[len("PROBE"):])
    print(f"[pip] python {out['python']}: viability {out['viability']}, "
          f"IC50 {out['ic50_uM']} uM, {out['stream_snapshots']} stream snapshots, "
          f"dish {out['dish_live']}, spheroid min O2 {out['spheroid_min_o2']} mmHg")

    failures = []
    if out["rdkit_present"]:
        failures.append("RDKit was present, so this run did not test a bare environment")
    if not (0.9 < out["viability"][0] <= 1.05 and out["viability"][-1] < 0.8):
        failures.append(f"dose-response looks wrong: {out['viability']}")
    if out["stream_snapshots"] != 3:
        failures.append(f"expected 3 stream snapshots, got {out['stream_snapshots']}")
    if not (out["dish_live"][0] == 20 and out["dish_live"][-1] >= 20):
        failures.append(f"dish did not grow: {out['dish_live']}")
    if not out["spheroid_min_o2"] < 100.0:
        failures.append("spheroid showed no oxygen gradient")

    for sub in (["cell-sim", "--help"], ["dish", "--help"], ["cell-stream", "--help"]):
        r = run([str(cellsim_exe), *sub], cwd=str(tmp))
        if r.returncode != 0 or "usage:" not in r.stdout:
            failures.append(f"`cellsim {' '.join(sub)}` failed: {r.stderr[-300:]}")

    if not keep:
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        print(f"[pip] kept {tmp}")
    for f in failures:
        print(f"  FAIL {f}")
    print("PASS: the cell simulator installs and runs with pip alone"
          if not failures else f"{len(failures)} pip-install gates FAILED")
    return 1 if failures else 0


def test_pip_install_cell_simulator():
    """`pip install cellsim` must be enough to simulate cells."""
    assert main() == 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--keep", action="store_true", help="leave the venv in place")
    sys.exit(main(ap.parse_args().keep))
