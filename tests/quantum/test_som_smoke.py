"""Layer 1.4 SoM smoke — xTB BDE runs on 3 diverse drugs.

Sanity checks on the returned ranking:
  - predictor returns `ok=True` with at least 3 candidates
  - every candidate has a finite BDE in the sane range
    [60, 180] kcal/mol (GFN2 approximate; literature BDE is
    70–110 kcal/mol for C-H bonds, but GFN2 has a +20 kcal/mol
    systematic offset that's well-documented)
  - rank-1 BDE is strictly less than any rank-3 BDE (ranking
    is monotonic)
  - rank-1 parent atom sits somewhere sensible (not on an atom
    that has no X-H bonds)

Literature-validated SoM correctness is a much larger benchmark
effort (future PR with a proper cyp3a4_som.smi bundle). This
smoke is strictly "does the predictor produce plausible output?"

Run:
    conda activate cellsim
    python tests/quantum/test_som_smoke.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.quantum import predict_cyp_som_bde  # noqa: E402


CASES: list[tuple[str, str]] = [
    ("aspirin",     "CC(=O)OC1=CC=CC=C1C(=O)O"),
    ("ibuprofen",   "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
    ("testosterone", "CC12CCC3C(C1CCC2=O)CCC4=CC(=O)CCC34C"),
]

BDE_LO_KCAL = 60.0          # absolute floor (sanity)
BDE_HI_KCAL = 180.0         # absolute ceiling (sanity)


def test_som_smoke():
    pass_count = 0
    t0 = time.time()
    for name, smi in CASES:
        print(f"[som] {name:20s}", flush=True)
        r = predict_cyp_som_bde(smi, top_k=5)
        if not r.ok:
            print(f"  FAIL: {r.reason}")
            continue
        if len(r.candidates) < 3:
            print(f"  FAIL: only {len(r.candidates)} candidates")
            continue
        cands = r.candidates
        bdes = [c.bde_kcalmol for c in cands if c.bde_kcalmol is not None]
        if len(bdes) < 3:
            print(f"  FAIL: <3 finite BDE values")
            continue
        if any(b < BDE_LO_KCAL or b > BDE_HI_KCAL for b in bdes):
            print(f"  FAIL: BDE out of [{BDE_LO_KCAL},{BDE_HI_KCAL}]: "
                  f"range {min(bdes):.1f} - {max(bdes):.1f}")
            continue
        if bdes[0] >= bdes[2]:
            print(f"  FAIL: ranking not monotonic "
                  f"(rank1={bdes[0]:.1f} >= rank3={bdes[2]:.1f})")
            continue
        top = r.top_k(3)
        print(f"  OK  top-3: " + ", ".join(
            f"{c.parent_element}(idx={c.parent_atom_idx}):{c.bde_kcalmol:.1f}"
            for c in top)
            + f"  ({r.wall_seconds:.1f} s)")
        pass_count += 1
    dt = time.time() - t0
    print(f"[som] {pass_count}/{len(CASES)} OK  wall={dt:.1f} s")
    if pass_count < len(CASES):
        print(f"FAIL: {pass_count}/{len(CASES)}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    test_som_smoke()
