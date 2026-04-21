"""Receptor PDB-ID auto-fetch resolver.

Extracted from src/dock/batch.py so dock-one, off-target, cyp-
inhibit, and any future receptor-taking subcommand can share the
same 'is this a file path or a 4-letter PDB ID?' logic without
duplication.
"""

from __future__ import annotations

import sys
from pathlib import Path


def resolve_receptor(receptor: str) -> str:
    """If `receptor` is a 4-char alphanumeric PDB ID that doesn't
    exist as a file, download it from RCSB into
    `data/receptors/<id>.pdb` and return that path. Otherwise
    return the argument unchanged.

    Never raises on a network failure — returns the original
    string and prints a diagnostic. Caller code downstream will
    emit its usual file-not-found error, preserving the existing
    behaviour for users who meant a local path that happens to
    be 4 chars long.
    """
    p = Path(receptor)
    if p.exists():
        return str(p)
    token = receptor.strip()
    if not (len(token) == 4 and token.isalnum()):
        return receptor

    from urllib.request import Request, urlopen

    out = Path("data/receptors") / f"{token.lower()}.pdb"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        return str(out)

    url = f"https://files.rcsb.org/download/{token.upper()}.pdb"
    print(f"[receptor] '{token}' looks like a PDB ID; fetching "
          f"{url}", flush=True)
    try:
        req = Request(url, headers={"User-Agent": "CellSim/1.0"})
        with urlopen(req, timeout=60) as resp:
            data = resp.read()
        out.write_bytes(data)
        print(f"[receptor] wrote {out}  "
              f"({len(data) / 1024:.1f} KB)", flush=True)
        return str(out)
    except Exception as e:
        print(f"[receptor] PDB fetch failed: {e}. Provide a local "
              "--receptor path.",
              file=sys.stderr, flush=True)
        return receptor
