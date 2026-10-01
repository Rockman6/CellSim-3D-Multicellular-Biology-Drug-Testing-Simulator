"""CellSim physics-prior cache.

SQLite + HDF5 memoisation of everything expensive in Campaign 1:
  - AM1-BCC partial charges (antechamber/sqm, ~20 s/compound)
  - xTB GFN2 results (single-point, reactive fragments)
  - Vina docking poses + ΔG
  - alchemical FEP ΔΔG (perses — pending)
  - classical MD equilibrated snapshots

Cache key: (content_hash, method, method_version, cellsim_version).
Content hashes are deterministic over the *physical meaning*
(canonical InChI Key for compounds, protein atom-hash for receptors)
so the same compound under a renamed SMILES still hits the cache.

Non-AI. No learned model; just memoisation of deterministic physics
calls.
"""

from .store import Cache, CacheEntry
from .hashing import compound_hash, receptor_hash

__all__ = [
    "Cache",
    "CacheEntry",
    "compound_hash",
    "receptor_hash",
]
