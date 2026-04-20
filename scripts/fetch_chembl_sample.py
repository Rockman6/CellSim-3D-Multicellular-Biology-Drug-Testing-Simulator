#!/usr/bin/env python3
"""Fetch a drug-like subset from ChEMBL for Layer 1.1 scale validation.

Uses the ChEMBL REST API (https://www.ebi.ac.uk/chembl/api/data/).
Filters for neutral, drug-like small molecules so AM1-BCC has a
fighting chance: MW 150–700, ≥ 1 ring, ≤ 60 heavy atoms, valid
canonical SMILES.

Writes `data/chembl/chembl_{N}.smi` (one SMILES + ChEMBL-ID per line,
tab-separated). `data/` is gitignored so the bulk download never
lands in a commit.

Typical use:
    python scripts/fetch_chembl_sample.py --n 1000
    python scripts/fetch_chembl_sample.py --n 10000 --out data/chembl/chembl_10k.smi
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PAGE_SIZE = 1000


def fetch_page(offset: int, limit: int) -> dict:
    params = {
        "limit": min(limit, PAGE_SIZE),
        "offset": offset,
        # Drug-like filters (structure type = molecule; small-molecule).
        "molecule_type": "Small molecule",
        "molecule_properties__mw_freebase__gte": 150,
        "molecule_properties__mw_freebase__lte": 700,
        "molecule_properties__heavy_atoms__lte": 60,
        "molecule_properties__num_ro5_violations__lte": 1,
    }
    url = f"{CHEMBL}?{urlencode(params)}"
    req = Request(url, headers={"User-Agent": "CellSim/1.0"})
    with urlopen(req, timeout=60) as resp:
        return json.load(resp)


def extract_rows(page: dict) -> list[tuple[str, str]]:
    rows = []
    for m in page.get("molecules", []):
        cid = m.get("molecule_chembl_id")
        structs = m.get("molecule_structures") or {}
        smi = structs.get("canonical_smiles")
        if cid and smi and "." not in smi:  # skip salt/mixture rows
            rows.append((smi, cid))
    return rows


def fetch_n(n: int, sleep_s: float = 0.2) -> list[tuple[str, str]]:
    collected: list[tuple[str, str]] = []
    seen_cids: set[str] = set()
    offset = 0
    while len(collected) < n:
        print(f"  [chembl] offset={offset} collected={len(collected)}",
              file=sys.stderr)
        try:
            page = fetch_page(offset, PAGE_SIZE)
        except Exception as e:
            print(f"  [chembl] page fetch failed at offset={offset}: {e}",
                  file=sys.stderr)
            time.sleep(2.0)
            continue
        rows = extract_rows(page)
        if not rows:
            break
        for smi, cid in rows:
            if cid not in seen_cids:
                seen_cids.add(cid)
                collected.append((smi, cid))
                if len(collected) >= n:
                    break
        offset += PAGE_SIZE
        time.sleep(sleep_s)
    return collected[:n]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=1000,
                    help="number of compounds to fetch (default 1000)")
    ap.add_argument("--out", type=Path,
                    help="output path (default "
                         "data/chembl/chembl_{n}.smi)")
    args = ap.parse_args()

    out = args.out or Path(
        f"data/chembl/chembl_{args.n}.smi")
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = fetch_n(args.n)
    with out.open("w") as f:
        f.write(f"# ChEMBL drug-like subset, n={len(rows)}\n")
        f.write("# Columns: smiles\tchembl_id\n")
        for smi, cid in rows:
            f.write(f"{smi}\t{cid}\n")
    print(f"[chembl] wrote {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()
