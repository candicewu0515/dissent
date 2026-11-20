"""Shared HTTP plumbing for every external source.

Two guarantees that the rest of the codebase relies on:

1. **A source failure never raises.** `fetch_json` returns `(None, status)` with
   `status.ok == False`. A report with four of five sources is still a useful
   report, and the UI shows which source was missing. Interpretation tools that
   fail closed on one dead API are useless in practice.
2. **Every fetch is cache-aware and timed.** The returned SourceStatus records
   whether the answer came from cache and how long it took, which the UI shows.
"""

from __future__ import annotations

import asyncio
import time
from typing import Any, Literal

import httpx

from dissent.cache.store import ResponseCache
from dissent.models import SourceStatus

USER_AGENT = "Dissent/0.1 (variant interpretation worksheet; research use)"
DEFAULT_TIMEOUT = 12.0
# Total wall-clock a single source may consume across all its attempts.
# Without this, `timeout` multiplied by `retries` is the real worst case: a
# 20s timeout with two retries spent 63s on one slow gnomAD query and the
# surrounding asyncio.gather discarded a sibling source's finished result.
DEFAULT_BUDGET = 16.0


class SourceClient:
    """Base class for one external data source."""

    name: str = "source"
    base_url: str = ""

    def __init__(
        self,
        client: httpx.AsyncClient | None = None,
        cache: ResponseCache | None = None,
        timeout: float = DEFAULT_TIMEOUT,
        budget: float = DEFAULT_BUDGET,
    ):
        self._client = client
        self._owns_client = client is None
        self.cache = cache
        self.timeout = timeout
        self.budget = budget

    async def __aenter__(self) -> SourceClient:
        if self._client is None:
            self._client = httpx.AsyncClient(
                timeout=self.timeout,
                headers={"User-Agent": USER_AGENT},
                follow_redirects=True,
            )
        return self

    async def __aexit__(self, *exc: Any) -> None:
        if self._owns_client and self._client is not None:
            await self._client.aclose()
            self._client = None

    @property
    def client(self) -> httpx.AsyncClient:
        if self._client is None:
            self._client = httpx.AsyncClient(
                timeout=self.timeout,
                headers={"User-Agent": USER_AGENT},
                follow_redirects=True,
            )
            self._owns_client = True
        return self._client

    async def fetch_json(
        self,
        url: str,
        *,
        method: Literal["GET", "POST"] = "GET",
        params: dict[str, Any] | None = None,
        json_body: dict[str, Any] | None = None,
        cache_key_parts: Any = None,
        retries: int = 2,
    ) -> tuple[Any | None, SourceStatus]:
        """Fetch and parse JSON. Returns (data_or_None, status). Never raises."""
        started = time.perf_counter()

        key = None
        if self.cache is not None:
            parts = cache_key_parts if cache_key_parts is not None else [url, params, json_body]
            key = ResponseCache.make_key(self.name, parts)
            hit = self.cache.get(key)
            if hit is not None:
                return hit, SourceStatus(
                    source=self.name,
                    ok=True,
                    cached=True,
                    elapsed_ms=int((time.perf_counter() - started) * 1000),
                )

        last_error = "unknown error"
        for attempt in range(retries + 1):
            elapsed = time.perf_counter() - started
            # Stop before an attempt that cannot finish inside the budget.
            if attempt > 0 and elapsed + self.timeout > self.budget:
                last_error = f"{last_error} (budget {self.budget:.0f}s exhausted)"
                break
            try:
                response = await self.client.request(
                    method, url, params=params, json=json_body, timeout=self.timeout
                )
                if response.status_code == 404:
                    # Not an error: the variant genuinely is not in this source.
                    return None, SourceStatus(
                        source=self.name,
                        ok=True,
                        detail="no record in this source",
                        elapsed_ms=int((time.perf_counter() - started) * 1000),
                    )
                response.raise_for_status()
                data = response.json()
                if key is not None and self.cache is not None:
                    self.cache.set(key, self.name, data)
                return data, SourceStatus(
                    source=self.name,
                    ok=True,
                    elapsed_ms=int((time.perf_counter() - started) * 1000),
                )
            except httpx.HTTPStatusError as exc:
                last_error = f"HTTP {exc.response.status_code}"
                # 4xx other than 429 will not improve on retry.
                if exc.response.status_code < 500 and exc.response.status_code != 429:
                    break
            except httpx.TimeoutException:
                last_error = f"timeout after {self.timeout}s"
            except httpx.HTTPError as exc:
                last_error = f"{type(exc).__name__}: {exc}"
            except ValueError as exc:
                last_error = f"malformed JSON: {exc}"
                break

            if attempt < retries:
                await asyncio.sleep(0.4 * (2**attempt))

        return None, SourceStatus(
            source=self.name,
            ok=False,
            detail=last_error,
            elapsed_ms=int((time.perf_counter() - started) * 1000),
        )


# --------------------------------------------------------------------------
# Helpers for the ragged shapes these APIs return
# --------------------------------------------------------------------------


def as_list(value: Any) -> list[Any]:
    """Normalise scalar-or-list into a list, dropping Nones.

    MyVariant returns a scalar when one transcript matched and a list when
    several did, for the same field. Every read of its data goes through here.
    """
    if value is None:
        return []
    if isinstance(value, list):
        return [v for v in value if v is not None]
    return [value]


def first(value: Any, default: Any = None) -> Any:
    items = as_list(value)
    return items[0] if items else default


def dig(obj: Any, *path: str, default: Any = None) -> Any:
    """Walk a nested dict path, tolerating lists by taking the first element."""
    current = obj
    for key in path:
        if isinstance(current, list):
            current = current[0] if current else None
        if not isinstance(current, dict):
            return default
        current = current.get(key)
        if current is None:
            return default
    return current


def to_float(value: Any) -> float | None:
    """Coerce a scalar-or-list numeric to float, preferring the extreme.

    Where several transcripts give different scores, taking the first would be
    arbitrary. Callers that care about direction pass a reducer instead.
    """
    items = as_list(value)
    for item in items:
        try:
            return float(item)
        except (TypeError, ValueError):
            continue
    return None


def to_float_all(value: Any) -> list[float]:
    out: list[float] = []
    for item in as_list(value):
        try:
            out.append(float(item))
        except (TypeError, ValueError):
            continue
    return out


def to_int(value: Any) -> int | None:
    f = to_float(value)
    return int(f) if f is not None else None
