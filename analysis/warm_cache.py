#!/usr/bin/env python3
"""Fetch everything the demo touches, so a recording does not wait on the internet.

A live demo over six public APIs is a demo of someone else's uptime. NCBI
rate-limits unkeyed clients, gnomAD's GraphQL endpoint has had minute-long
responses, and a recording that stalls for forty seconds is a recording that
gets re-shot. This fills the response cache ahead of time, so every query in
the demo resolves from disk.

It is deliberately *not* a fixture dump. The cache holds the same responses the
live code would fetch, keyed the same way, and it expires on the same schedule
as any other entry -- so nothing shown on camera is a recording of something
other than what the tool does. Run it shortly before filming.

    python analysis/warm_cache.py                 # the demo set
    python analysis/warm_cache.py --check         # report coverage, fetch nothing
    python analysis/warm_cache.py --also "TP53 R175H" "rs334"

Exit status is 0 only if every query resolved, so it can gate a recording.
"""

from __future__ import annotations

import argparse
import asyncio
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import httpx  # noqa: E402

from dissent.api.app import EXAMPLES, PANEL_EXAMPLES  # noqa: E402
from dissent.cache.store import ResponseCache  # noqa: E402
from dissent.report import interpret  # noqa: E402
from dissent.sources.base import USER_AGENT  # noqa: E402

# Four at a time, matching the panel's own limit: these are unpaid public
# endpoints and a warm-up that trips a rate limit leaves the cache worse than
# it found it.
CONCURRENCY = 4


def demo_queries() -> list[str]:
    """Every query the demo can reach, in the order the script uses them."""
    seen: set[str] = set()
    out: list[str] = []
    for q in [e["query"] for e in EXAMPLES]:
        if q.upper() not in seen:
            seen.add(q.upper())
            out.append(q)
    for panel in PANEL_EXAMPLES:
        for q in panel["queries"]:
            if q.upper() not in seen:
                seen.add(q.upper())
                out.append(q)
    # The identifier forms the demo types by hand, which take different code
    # paths through the parser than a gene-and-protein-change query does.
    for q in ["rs113488022", "chr7:140753336A>T", "41815", "rs1800562"]:
        if q.upper() not in seen:
            seen.add(q.upper())
            out.append(q)
    return out


async def warm(queries: list[str], cache: ResponseCache, check_only: bool) -> int:
    if not cache.enabled:
        print(f"! cache is disabled: {cache.disabled_reason}", file=sys.stderr)
        print("  Nothing can be warmed. Check DISSENT_CACHE_PATH.", file=sys.stderr)
        return 2

    sem = asyncio.Semaphore(CONCURRENCY)
    failures: list[tuple[str, str]] = []
    slow: list[tuple[str, float]] = []

    async with httpx.AsyncClient(
        timeout=30.0, headers={"User-Agent": USER_AGENT}, follow_redirects=True
    ) as client:

        async def one(q: str) -> None:
            async with sem:
                started = time.monotonic()
                try:
                    report = await asyncio.wait_for(
                        interpret(q, cache=cache, client=client), timeout=90
                    )
                except Exception as exc:  # noqa: BLE001
                    failures.append((q, f"{type(exc).__name__}: {exc}"))
                    print(f"  FAIL  {q:24s} {type(exc).__name__}: {exc}")
                    return

                elapsed = time.monotonic() - started
                bundle = report.bundle
                resolved = bool(bundle and (bundle.clinvar or bundle.consequence)) or bool(
                    report.query.chrom and report.query.pos
                )
                if not resolved:
                    # Not a warm-up failure: some demo inputs do not resolve on
                    # purpose. Flagged rather than counted against the run.
                    print(f"  ----  {q:24s} {elapsed:5.1f}s  does not resolve (expected for some)")
                    return

                # After a successful pass the second read should be instant. If
                # it is not, the cache is not being hit and the demo will still
                # wait on the network.
                if not check_only and elapsed > 2.0:
                    slow.append((q, elapsed))
                mark = "warm" if check_only else "ok  "
                print(f"  {mark}  {q:24s} {elapsed:5.1f}s  {report.classification.call}")

        await asyncio.gather(*(one(q) for q in queries))

    print()
    print(f"{len(queries) - len(failures)}/{len(queries)} resolved or deliberately did not")
    stats = cache.stats()
    if isinstance(stats, dict) and "enabled" not in stats:
        total = sum(stats.values())
        print(f"cache now holds {total} responses: " + ", ".join(
            f"{src} {n}" for src, n in sorted(stats.items())
        ))

    if slow:
        print()
        print("Still slow on the second pass — the cache is not being hit for these:")
        for q, secs in slow:
            print(f"  {q:24s} {secs:.1f}s")

    if failures:
        print()
        print("Failed, and would fail on camera:", file=sys.stderr)
        for q, why in failures:
            print(f"  {q:24s} {why}", file=sys.stderr)
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Warm the response cache before recording a demo."
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="run the queries and report timing without treating slowness as a problem",
    )
    parser.add_argument("--also", nargs="*", default=[], help="extra queries to warm")
    args = parser.parse_args(argv)

    queries = demo_queries() + list(args.also)
    cache = ResponseCache()
    print(f"Warming {len(queries)} queries into {cache.path}")
    print()
    rc = asyncio.run(warm(queries, cache, args.check))
    if rc == 0:
        print()
        print("Every query resolves from cache. Safe to record.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
