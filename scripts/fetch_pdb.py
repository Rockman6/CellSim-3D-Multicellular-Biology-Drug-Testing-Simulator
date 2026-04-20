#!/usr/bin/env python3
"""Fetch a PDB structure from RCSB and stash it under data/md/.

Usage:
    python scripts/fetch_pdb.py 1UBQ
    python scripts/fetch_pdb.py 1L2Y --out benchmarks/md/1l2y.pdb

`data/` is gitignored; `benchmarks/` is not — pass `--out` when the
PDB should ship as a regression gate input.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.request import Request, urlopen


def fetch_pdb(pdb_id: str, out: Path) -> None:
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    req = Request(url, headers={"User-Agent": "CellSim/1.0"})
    with urlopen(req, timeout=60) as resp:
        data = resp.read()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(data)
    print(f"[pdb] wrote {out}  ({len(data)/1024:.1f} KB)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb_id")
    ap.add_argument("--out", type=Path,
                    help="output path (default data/md/<id>.pdb)")
    args = ap.parse_args()

    out = args.out or Path(f"data/md/{args.pdb_id.lower()}.pdb")
    fetch_pdb(args.pdb_id, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
