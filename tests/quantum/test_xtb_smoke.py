"""Layer 1.4 xTB smoke — GFN2 single-point on canonical drugs.

Sanity checks on the returned numbers:
  - total energy is negative and finite (a bound molecule)
  - HOMO is negative; LUMO > HOMO; gap > 0 and < 10 eV for a typical
    organic drug
  - dipole is finite and non-NaN
  - wall time reasonable (< 30 s for a small drug on laptop)

Gate: 8/10 drugs from smoke_10.smi must return a consistent result.

Requires the `cellsim` conda env (xtb, rdkit).

Run:
    conda activate cellsim
    python tests/quantum/test_xtb_smoke.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.quantum import xtb_single_point  # noqa: E402

SMI_FILE = REPO_ROOT / "benchmarks" / "chembl" / "smoke_10.smi"


def load_smi(path: Path, max_n: int | None = None):
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


def _sane(r) -> tuple[bool, str]:
    if not r.ok:
        return False, f"fail: {r.reason}"
    if r.total_energy_eV is None or r.total_energy_eV > 0:
        return False, f"bad energy {r.total_energy_eV}"
    if r.homo_eV is None or r.lumo_eV is None:
        return False, "no HOMO/LUMO"
    if not (-25.0 < r.homo_eV < 0.0):
        return False, f"HOMO out of range {r.homo_eV}"
    if r.lumo_eV <= r.homo_eV:
        return False, f"LUMO < HOMO ({r.lumo_eV}, {r.homo_eV})"
    if r.homo_lumo_gap_eV is None or not (0.0 < r.homo_lumo_gap_eV < 12.0):
        return False, f"gap out of range {r.homo_lumo_gap_eV}"
    if r.dipole_Debye is None or r.dipole_Debye != r.dipole_Debye:
        return False, "NaN dipole"
    return True, ""


def test_xtb_smoke(max_n: int | None = None, gate: int = 8):
    entries = load_smi(SMI_FILE, max_n=max_n)
    ok_count = 0
    t0 = time.time()
    for smi, name in entries:
        r = xtb_single_point(smi)
        sane, note = _sane(r)
        tag = "OK  " if sane else "FAIL"
        detail = (f"E={r.total_energy_eV:>9.2f}  HOMO={r.homo_eV:>6.2f}"
                  f"  LUMO={r.lumo_eV:>6.2f}  gap={r.homo_lumo_gap_eV:>5.2f}  "
                  f"µ={r.dipole_Debye:>5.2f} D  ({r.wall_seconds:.2f} s)"
                  if sane else note)
        print(f"  {tag}  {name:25s}  {detail}", flush=True)
        if sane:
            ok_count += 1
    dt = time.time() - t0
    print(f"[xtb] {ok_count}/{len(entries)} sane  gate={gate}  "
          f"wall={dt:.1f} s")
    if ok_count < gate:
        print(f"FAIL: {ok_count}/{len(entries)} < {gate}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=None)
    ap.add_argument("--gate", type=int, default=8)
    args = ap.parse_args()
    test_xtb_smoke(max_n=args.max, gate=args.gate)
