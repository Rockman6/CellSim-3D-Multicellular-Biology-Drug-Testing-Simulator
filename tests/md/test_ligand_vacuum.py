"""Layer 1.2 smoke test: vacuum Langevin MD of canonical drugs.

Gate:
    - Every compound that parametrised in Layer 1.1 must also
      integrate for 5 000 × 2 fs = 10 ps without NaN, with final
      temperature within 50 K of 300 K setpoint and RMSD < 10 Å.

Defaults: full 10-compound set with a ≥ 8 pass gate. CI uses a
smaller subset (`--max 3 --gate 2`) so the run fits into the PR
budget (AM1-BCC is ~20 s/compound).

Requires the `cellsim` conda env (OpenMM + OpenFF + AmberTools
+ AMBERHOME set).

Run:
    conda activate cellsim
    python tests/md/test_ligand_vacuum.py                 # full
    python tests/md/test_ligand_vacuum.py --max 3 --gate 2   # CI
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.md import simulate_ligand  # noqa: E402

SMI_FILE = REPO_ROOT / "benchmarks" / "chembl" / "smoke_10.smi"


def load_smi(path: Path, max_n: int | None = None) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split(None, 1)
        entries.append((parts[0], parts[1] if len(parts) > 1 else parts[0]))
        if max_n is not None and len(entries) >= max_n:
            break
    return entries


def run(max_n: int | None, gate: int) -> int:
    entries = load_smi(SMI_FILE, max_n=max_n)
    results = []
    t0 = time.time()
    for smi, name in entries:
        print(f"  [md] {name:25s} {smi}", flush=True)
        r = simulate_ligand(smi, n_steps=5000, dt_fs=2.0,
                            save_every=500, temperature_K=300.0,
                            friction_per_ps=2.0)
        results.append((name, smi, r))
        print(f"      -> {r.summary()}", flush=True)
    dt = time.time() - t0
    ok_count = sum(1 for _, _, r in results if r.ok)
    n = len(entries)

    print(f"\n[md] pass-rate {ok_count}/{n}  gate={gate}  wall={dt:.1f} s")
    if ok_count < gate:
        print(f"FAIL: {ok_count}/{n} < {gate}")
        return 1
    print("PASS")
    return 0


def test_ligand_vacuum_md():
    """≥ 8/10 canonical drugs must run 10 ps MD without blowing up."""
    assert run(None, 8) == 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max", type=int, default=None,
                    help="cap compound count (default: all 10)")
    ap.add_argument("--gate", type=int, default=8,
                    help="minimum passes required (default 8)")
    args = ap.parse_args()
    sys.exit(run(args.max, args.gate))
