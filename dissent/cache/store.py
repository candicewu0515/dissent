"""SQLite-backed HTTP response cache.

Two reasons this exists rather than being an optimisation:

1. NCBI rate-limits unkeyed clients to ~3 requests/second. A cold interpretation
   touches several NCBI endpoints; without a cache, repeated use trips the limit.
2. A live demo must not depend on network conditions. Pre-warming the cache for
   the demo variants makes the recording deterministic.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import sqlite3
import time
from pathlib import Path
from typing import Any

DEFAULT_TTL_SECONDS = 60 * 60 * 24 * 7  # one week; these sources change slowly


class ResponseCache:
    def __init__(self, path: str | Path = ".cache/dissent.sqlite", ttl: int = DEFAULT_TTL_SECONDS):
        self.path = Path(path)
        self.ttl = ttl
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    @contextlib.contextmanager
    def _connect(self):
        """Yield a connection and always close it.

        `with sqlite3.connect(...)` commits but does not close. Each cache read
        and write opened a fresh connection, so a long-running server leaked a
        file descriptor per operation -- around ten per interpretation.
        """
        conn = sqlite3.connect(self.path, timeout=10)
        try:
            conn.execute("PRAGMA journal_mode=WAL")
            with conn:
                yield conn
        finally:
            conn.close()

    def _init_db(self) -> None:
        with self._connect() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS responses (
                    key        TEXT PRIMARY KEY,
                    source     TEXT NOT NULL,
                    payload    TEXT NOT NULL,
                    stored_at  REAL NOT NULL
                )
                """
            )
            conn.execute("CREATE INDEX IF NOT EXISTS idx_source ON responses(source)")

    @staticmethod
    def make_key(source: str, request: Any) -> str:
        blob = json.dumps(request, sort_keys=True, default=str)
        digest = hashlib.sha256(blob.encode()).hexdigest()[:32]
        return f"{source}:{digest}"

    def get(self, key: str) -> Any | None:
        with self._connect() as conn:
            row = conn.execute(
                "SELECT payload, stored_at FROM responses WHERE key = ?", (key,)
            ).fetchone()
        if row is None:
            return None
        payload, stored_at = row
        if time.time() - stored_at > self.ttl:
            return None
        try:
            return json.loads(payload)
        except json.JSONDecodeError:
            return None

    def set(self, key: str, source: str, value: Any) -> None:
        with self._connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO responses (key, source, payload, stored_at) "
                "VALUES (?, ?, ?, ?)",
                (key, source, json.dumps(value), time.time()),
            )

    def clear(self, source: str | None = None) -> int:
        with self._connect() as conn:
            if source:
                cur = conn.execute("DELETE FROM responses WHERE source = ?", (source,))
            else:
                cur = conn.execute("DELETE FROM responses")
            return cur.rowcount

    def stats(self) -> dict[str, int]:
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT source, COUNT(*) FROM responses GROUP BY source"
            ).fetchall()
        return {src: n for src, n in rows}
