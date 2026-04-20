"""Layer 1.1 scale test — parametrise N compounds, report pass-rate.

Consumes a `.smi` file (tab- or whitespace-separated, one SMILES per
non-comment line, optional ID in column 2) and drives
`parametrize_smiles` over the whole set, optionally in parallel.
Emits a pass-rate summary and a failure taxonomy.

Gates (auto-selected by env):
  - Full tier (OpenFF present): ≥ 99 % parameterise to openmm.System.
  - RDKit tier (OpenFF absent):  ≥ 99 % embed + 3D coords.

Layer 1.1 exit criterion is the full tier on 10 k ChEMBL compounds.
This file is the infrastructure; the 10 k run uses a dataset
fetched by `scripts/fetch_chembl_sample.py --n 10000` (not in-tree,
because `data/` is gitignored).

Usage:
    # quick
    python tests/chem/test_parametrize_scale.py --smi benchmarks/chembl/smoke_10.smi

    # 1 k
    python scripts/fetch_chembl_sample.py --n 1000
    python tests/chem/test_parametrize_scale.py \\
        --smi data/chembl/chembl_1000.smi --workers 4

    # full 10 k (hours on CPU; AM1-BCC dominates)
    python scripts/fetch_chembl_sample.py --n 10000
    python tests/chem/test_parametrize_scale.py \\
        --smi data/chembl/chembl_10000.smi --workers 8 --charge am1bcc
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import time
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.chem import parametrize_smiles  # noqa: E402


def load_smi(path: Path, max_n: int | None = None) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split(None, 1)
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else smi
        entries.append((smi, name))
        if max_n is not None and len(entries) >= max_n:
            break
    return entries


def _worker(args):
    smi, charge, ff, seed = args
    r = parametrize_smiles(smi, charge_method=charge,
                           ff_name=ff, random_seed=seed)
    return (r.smiles, r.ok, r.reason, r.formula, r.n_atoms,
            (r.positions_nm is not None))


def _summarize_reason(reason: str) -> str:
    """Short bucket key for the failure taxonomy."""
    if not reason:
        return ""
    first = reason.split(":")[0].strip()
    return first[:80]


def run(smi_path: Path, *, charge: str, ff: str, workers: int,
        seed: int, max_n: int | None, gate_pct: float) -> int:
    entries = load_smi(smi_path, max_n=max_n)
    n = len(entries)
    if n == 0:
        print(f"[scale] no compounds in {smi_path}", file=sys.stderr)
        return 2

    print(f"[scale] {n} compounds from {smi_path}  "
          f"charge={charge}  ff={ff}  workers={workers}",
          file=sys.stderr)
    t0 = time.time()

    tasks = [(smi, charge, ff, seed) for smi, _ in entries]
    results: list = []
    last_heartbeat = t0

    def _log(idx: int, r) -> None:
        nonlocal last_heartbeat
        tag = "OK  " if r[1] else "FAIL"
        short_reason = r[2].split("[")[0].split(":")[0].strip()[:40] \
            if r[2] else ""
        print(f"  [{idx:>5d}/{n}] {tag}  {r[3] or '-':<22s}  "
              f"{short_reason}", flush=True)
        last_heartbeat = time.time()

    if workers == 1:
        for i, t in enumerate(tasks, 1):
            r = _worker(t)
            results.append(r)
            _log(i, r)
    else:
        with mp.Pool(workers) as pool:
            for i, r in enumerate(
                    pool.imap_unordered(_worker, tasks, chunksize=1), 1):
                results.append(r)
                _log(i, r)

    dt = time.time() - t0

    full_ok = sum(1 for _, ok, *_ in results if ok)
    rdkit_ok = sum(1 for _, ok, _, _, _, has_geom in results if ok or has_geom)
    reasons = Counter(_summarize_reason(r[2]) for r in results if not r[1])

    # Pick tier
    try:
        import openff.toolkit  # noqa: F401
        tier = "full"
        passed = full_ok
    except ImportError:
        tier = "rdkit"
        passed = rdkit_ok

    pct = passed / n * 100.0
    per_sec = n / dt if dt > 0 else float("inf")

    print(f"\n[scale] tier={tier}  pass={passed}/{n} "
          f"({pct:.2f} %)  gate={gate_pct:.2f} %  "
          f"throughput={per_sec:.2f} compounds/s  "
          f"wall={dt:.1f} s")
    if reasons:
        print("[scale] failure taxonomy (top 10):")
        for r, c in reasons.most_common(10):
            print(f"  {c:5d}  {r}")

    ok = pct >= gate_pct
    print(f"[scale] verdict: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smi", type=Path, required=True,
                    help="input .smi file (tab or whitespace separated)")
    ap.add_argument("--charge", default="am1bcc",
                    help="charge method (am1bcc|gasteiger|mmff94|zeros|"
                         "formal_charge)")
    ap.add_argument("--ff", default="openff-2.1.0.offxml")
    ap.add_argument("--workers", type=int, default=1,
                    help="multiprocessing workers (1 = sequential)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max", type=int, default=None,
                    help="cap number of compounds consumed")
    ap.add_argument("--gate", type=float, default=99.0,
                    help="pass-rate gate %% (default 99.0)")
    args = ap.parse_args()
    sys.exit(run(args.smi, charge=args.charge, ff=args.ff,
                 workers=args.workers, seed=args.seed,
                 max_n=args.max, gate_pct=args.gate))


if __name__ == "__main__":
    main()
