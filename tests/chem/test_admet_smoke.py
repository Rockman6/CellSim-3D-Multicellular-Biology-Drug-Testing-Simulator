"""Layer 1.3-ish ADMET descriptor smoke (non-AI).

Runs compute_admet on smoke_10.smi and sanity-checks every
returned record:
  - ok=True
  - MW / logP / TPSA are finite numbers in plausible ranges
  - Ro5 flag is consistent with the violation count
  - QED in [0, 1]
  - logS finite, classification non-empty

Gate: all 10 drugs must return ok + sane numbers. No Ro5 pass/fail
gate per compound — the test validates the infrastructure, not
the drug-likeness of these specific compounds.

Run:
    conda activate cellsim
    python tests/chem/test_admet_smoke.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.chem import compute_admet  # noqa: E402


SMI_FILE = REPO_ROOT / "benchmarks" / "chembl" / "smoke_10.smi"


def _sane(r) -> tuple[bool, str]:
    if not r.ok:
        return False, r.reason
    if r.MW is None or not (10.0 < r.MW < 2000.0):
        return False, f"MW out of range {r.MW}"
    if r.logP is None or math.isnan(r.logP) or abs(r.logP) > 20.0:
        return False, f"logP out of range {r.logP}"
    if r.tpsa is None or r.tpsa < 0.0 or r.tpsa > 500.0:
        return False, f"TPSA out of range {r.tpsa}"
    if r.hba is None or r.hbd is None or r.rotb is None:
        return False, "missing H-bond / rotb counts"
    if r.ro5_pass is None or r.ro5_violations is None:
        return False, "missing Ro5 flags"
    if r.ro5_pass != (r.ro5_violations == 0):
        return False, (f"Ro5 inconsistency: pass={r.ro5_pass} "
                       f"violations={r.ro5_violations}")
    if r.qed is None or not (0.0 <= r.qed <= 1.0):
        return False, f"QED out of [0,1]: {r.qed}"
    if r.logS_ESOL is None or math.isnan(r.logS_ESOL):
        return False, "logS is NaN / missing"
    if not r.solubility_class:
        return False, "solubility_class empty"
    return True, ""


def load_smi(path: Path):
    entries: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split(None, 1)
        entries.append((parts[0], parts[1] if len(parts) > 1 else parts[0]))
    return entries


def test_admet_smoke():
    entries = load_smi(SMI_FILE)
    ok = 0
    for smi, name in entries:
        r = compute_admet(smi)
        sane, note = _sane(r)
        if sane:
            ok += 1
            print(f"  OK  {name:25s} {r.summary()[7:]}")  # strip "[OK]   "
        else:
            print(f"  FAIL {name:25s} {note}")
    print(f"[admet] {ok}/{len(entries)}  gate=10")
    if ok < len(entries):
        print(f"FAIL: {ok}/{len(entries)}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    test_admet_smoke()
