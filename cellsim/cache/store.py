"""SQLite-backed physics-prior cache store.

One `cache.sqlite` file per CellSim workspace. Entries are keyed by
`(content_hash, method)`; the value is an arbitrary JSON-serialisable
dict produced by a physics calculation (ParametrizeResult, XtbResult,
DockingResult top pose, …).

Design decisions (non-AI, non-fancy):

  - **SQLite only for this MVP.** HDF5 for large binary tensors can
    be added later; all current CellSim outputs fit in JSON (< 1 MB
    per entry even for docking with 9 poses).
  - **Content-addressed.** Keys are hashed physical inputs, not
    user-chosen names, so a renamed compound still hits the cache.
  - **Method version baked into the key.** Upgrading the force
    field (e.g. Sage 2.1 → 2.2) creates new keys; old entries are
    never served for the new method. No cache poisoning.
  - **Serialisation is plain JSON.** Readable by any tool;
    diffable in git if a user wants to commit a calibration cache.
  - **Concurrency: single-writer, multi-reader.** SQLite's default
    is fine; each worker opens its own connection.
"""

from __future__ import annotations

import json
import logging
import sqlite3
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterable, Optional

logger = logging.getLogger(__name__)


_SCHEMA = """
CREATE TABLE IF NOT EXISTS entries (
    content_hash   TEXT NOT NULL,
    method         TEXT NOT NULL,
    value_json     TEXT NOT NULL,
    cellsim_commit TEXT,
    created_at     REAL NOT NULL,
    PRIMARY KEY (content_hash, method)
);
CREATE INDEX IF NOT EXISTS idx_entries_created ON entries (created_at);
CREATE INDEX IF NOT EXISTS idx_entries_method  ON entries (method);
"""


@dataclass
class CacheEntry:
    """One cache row."""

    content_hash: str
    method: str
    value: Any = field(default_factory=dict)
    cellsim_commit: Optional[str] = None
    created_at: Optional[float] = None

    def as_dict(self) -> dict:
        return asdict(self)


class Cache:
    """Thin SQLite wrapper for content-addressed physics memoisation."""

    def __init__(self, path: str | Path,
                 cellsim_commit: Optional[str] = None):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.path))
        self._conn.executescript(_SCHEMA)
        self._conn.commit()
        self.cellsim_commit = cellsim_commit

    def close(self) -> None:
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    # ---------- read / write ----------

    def get(self, content_hash: str, method: str) -> Optional[CacheEntry]:
        row = self._conn.execute(
            "SELECT value_json, cellsim_commit, created_at "
            "FROM entries WHERE content_hash = ? AND method = ?",
            (content_hash, method)).fetchone()
        if row is None:
            return None
        try:
            value = json.loads(row[0])
        except json.JSONDecodeError:
            logger.warning("cache row had malformed JSON; ignoring")
            return None
        return CacheEntry(
            content_hash=content_hash, method=method,
            value=value, cellsim_commit=row[1], created_at=row[2])

    def put(self, content_hash: str, method: str,
             value: Any, *, cellsim_commit: Optional[str] = None
             ) -> None:
        commit = cellsim_commit or self.cellsim_commit
        value_json = json.dumps(value, default=_json_default)
        self._conn.execute(
            "INSERT OR REPLACE INTO entries "
            "(content_hash, method, value_json, cellsim_commit, "
            " created_at) VALUES (?, ?, ?, ?, ?)",
            (content_hash, method, value_json, commit, time.time()))
        self._conn.commit()

    def delete(self, content_hash: str, method: str) -> bool:
        cur = self._conn.execute(
            "DELETE FROM entries WHERE content_hash = ? AND method = ?",
            (content_hash, method))
        self._conn.commit()
        return cur.rowcount > 0

    def has(self, content_hash: str, method: str) -> bool:
        row = self._conn.execute(
            "SELECT 1 FROM entries WHERE content_hash = ? AND method = ?",
            (content_hash, method)).fetchone()
        return row is not None

    # ---------- bulk ops / diagnostics ----------

    def iter_entries(self, method: Optional[str] = None
                      ) -> Iterable[CacheEntry]:
        if method:
            cur = self._conn.execute(
                "SELECT content_hash, method, value_json, "
                "cellsim_commit, created_at FROM entries "
                "WHERE method = ? ORDER BY created_at",
                (method,))
        else:
            cur = self._conn.execute(
                "SELECT content_hash, method, value_json, "
                "cellsim_commit, created_at FROM entries "
                "ORDER BY created_at")
        for row in cur:
            try:
                value = json.loads(row[2])
            except json.JSONDecodeError:
                continue
            yield CacheEntry(
                content_hash=row[0], method=row[1],
                value=value, cellsim_commit=row[3],
                created_at=row[4])

    def stats(self) -> dict:
        n = self._conn.execute(
            "SELECT COUNT(*) FROM entries").fetchone()[0]
        by_method = dict(self._conn.execute(
            "SELECT method, COUNT(*) FROM entries GROUP BY method"
        ).fetchall())
        size_bytes = self.path.stat().st_size if self.path.exists() else 0
        return {
            "path": str(self.path),
            "n_entries": n,
            "entries_by_method": by_method,
            "size_bytes": size_bytes,
        }

    def invalidate_method(self, method_prefix: str) -> int:
        """Drop every entry whose method string starts with `method_
        prefix` (used on FF / tool-version bumps to clear stale rows
        in one call)."""
        cur = self._conn.execute(
            "DELETE FROM entries WHERE method LIKE ?",
            (f"{method_prefix}%",))
        self._conn.commit()
        return cur.rowcount


def _json_default(o: Any) -> Any:
    """JSON fallback for common numpy / pathlib types."""
    try:
        import numpy as np
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
    except ImportError:
        pass
    if isinstance(o, Path):
        return str(o)
    if hasattr(o, "as_dict"):
        return o.as_dict()
    raise TypeError(f"cannot serialise {type(o).__name__}")
