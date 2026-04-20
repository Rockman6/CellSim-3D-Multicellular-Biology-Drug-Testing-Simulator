"""Layer 1.2 smoke test: vacuum Langevin MD of 10 canonical drugs.

Gate:
    - Every compound that parametrised in Layer 1.1 must also
      integrate for 5 000 × 2 fs = 10 ps without NaN, with final
      temperature within 50 K of 300 K setpoint and RMSD < 10 Å.

This is the Layer-1.1 -> Layer-1.2 smoke bridge. Requires the
`cellsim` conda env (OpenMM + OpenFF + AmberTools + AMBERHOME set).

Run:
    conda activate cellsim
    python tests/md/test_ligand_vacuum.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.md import simulate_ligand  # noqa: E402

SMI_FILE = REPO_ROOT / "benchmarks" / "chembl" / "smoke_10.smi"


def load_smi(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split(None, 1)
        entries.append((parts[0], parts[1] if len(parts) > 1 else parts[0]))
    return entries


def test_ligand_vacuum_md():
    """≥ 8/10 canonical drugs must run 10 ps MD without blowing up."""
    entries = load_smi(SMI_FILE)
    results = []
    t0 = time.time()
    for smi, name in entries:
        print(f"  [md] {name:25s} {smi}", flush=True)
        r = simulate_ligand(smi, n_steps=5000, dt_fs=2.0,
                            save_every=500, temperature_K=300.0)
        results.append((name, smi, r))
        print(f"      -> {r.summary()}", flush=True)
    dt = time.time() - t0
    ok_count = sum(1 for _, _, r in results if r.ok)

    print(f"\n[md] pass-rate {ok_count}/{len(entries)}  wall={dt:.1f} s")
    assert ok_count >= 8, (
        f"MD smoke {ok_count}/{len(entries)} < 8/10; "
        "see stdout for per-compound diagnosis")


if __name__ == "__main__":
    try:
        test_ligand_vacuum_md()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
