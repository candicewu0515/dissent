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
import logging
import os
import sqlite3
import time
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

DEFAULT_TTL_SECONDS = 60 * 60 * 24 * 7  # one week; these sources change slowly

# Relative to the working directory by default, which is fine for a checkout
# and wrong for a container with a read-only root. The environment variable is
# how a deployment points it at a writable volume.
DEFAULT_PATH = os.environ.get("DISSENT_CACHE_PATH", ".cache/dissent.sqlite")


class ResponseCache:
    def __init__(self, path: str | Path | None = None, ttl: int = DEFAULT_TTL_SECONDS):
        self.path = Path(path if path is not None else DEFAULT_PATH)
        self.ttl = ttl
        # A cache is an optimisation. It must never be able to fail a request,
        # so every failure below turns it off instead of raising: a read-only
        # filesystem, a full disk, a directory that cannot be created. The
        # server then runs slower and hits the public APIs more often, which is
        # a far better outcome than refusing to interpret a variant.
        self.enabled = True
        self.disabled_reason: str | None = None
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self._init_db()
        except (OSError, sqlite3.Error) as exc:
            self._disable(exc)

    def _disable(self, exc: Exception) -> None:
        if self.enabled:
            self.enabled = False
            self.disabled_reason = f"{type(exc).__name__}: {exc}"
            log.warning(
                "Response cache disabled (%s at %s); continuing without it.",
                self.disabled_reason,
                self.path,
            )

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
        if not self.enabled:
            return None
        try:
            with self._connect() as conn:
                row = conn.execute(
                    "SELECT payload, stored_at FROM responses WHERE key = ?", (key,)
                ).fetchone()
        except sqlite3.Error as exc:
            # A failed read is a miss, not an error. The caller fetches.
            self._disable(exc)
            return None
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
        if not self.enabled:
            return
        try:
            with self._connect() as conn:
                conn.execute(
                    "INSERT OR REPLACE INTO responses (key, source, payload, stored_at) "
                    "VALUES (?, ?, ?, ?)",
                    (key, source, json.dumps(value), time.time()),
                )
        except sqlite3.Error as exc:
            # The disk will not take it. Turn the cache off; the response we
            # already have is still good.
            self._disable(exc)
        except (TypeError, ValueError) as exc:
            # This one payload will not serialise. That is a fault in one
            # response, not in the cache, so the cache stays on.
            log.warning("Could not cache a %s response: %s", source, exc)

    def clear(self, source: str | None = None) -> int:
        if not self.enabled:
            return 0
        with self._connect() as conn:
            if source:
                cur = conn.execute("DELETE FROM responses WHERE source = ?", (source,))
            else:
                cur = conn.execute("DELETE FROM responses")
            return cur.rowcount

    def stats(self) -> dict[str, Any]:
        """Row counts per source, or why there are none.

        Surfaced on /api/health: a cache that silently turned itself off should
        be visible, because the only other symptom is the server being slow.
        """
        if not self.enabled:
            return {"enabled": False, "reason": self.disabled_reason or "unknown"}
        try:
            with self._connect() as conn:
                rows = conn.execute(
                    "SELECT source, COUNT(*) FROM responses GROUP BY source"
                ).fetchall()
        except sqlite3.Error as exc:
            self._disable(exc)
            return {"enabled": False, "reason": self.disabled_reason or str(exc)}
        return {src: n for src, n in rows}
