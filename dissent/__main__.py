"""Command line entry point.

    python -m dissent serve                 start the web UI
    python -m dissent "BRAF V600E"          print a report to the terminal
"""

from __future__ import annotations

import argparse
import asyncio
import json
import sys

from dissent.cache.store import ResponseCache
from dissent.report import interpret

_SEVERITY_MARK = {"high": "!!", "medium": "! ", "low": "  "}
_STATUS_MARK = {"met": "[MET]", "not_met": "[ - ]", "insufficient": "[ ? ]", "manual": "[MAN]"}


def _print_report(report, width: int = 80) -> None:
    rule = "=" * width
    thin = "-" * width

    print(rule)
    print(report.display_name or report.query.raw)
    q = report.query
    coords = f"{q.chrom}:{q.pos}{q.ref}>{q.alt}" if q.chrom and q.pos else "(unresolved)"
    print(f"{report.gene_symbol or '?'}   {coords}   {q.rsid or ''}   [{q.assembly}]")
    print(rule)

    print(f"\nDRAFT CLASSIFICATION: {report.classification.call}")
    print(f"  {report.classification.reasoning}\n")

    if report.conflicts:
        print(thin)
        print(f"CONFLICTS ({len(report.conflicts)})")
        print(thin)
        for c in report.conflicts:
            print(f"{_SEVERITY_MARK.get(c.severity, '  ')} [{c.severity.upper()}] {c.headline}")
            print(f"     {c.detail}\n")
    else:
        print("No conflicts detected.\n")

    print(thin)
    counts: dict[str, int] = {}
    for c in report.criteria:
        counts[c.status] = counts.get(c.status, 0) + 1
    summary = "  ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"ACMG WORKSHEET   {summary}")
    print(thin)
    for c in report.criteria:
        if c.status in ("met", "insufficient"):
            strength = f" ({c.strength.value})" if c.strength else ""
            print(f"{_STATUS_MARK[c.status]} {c.code}{strength}  {c.name}")
            print(f"      {c.rationale}")
            if c.what_would_change:
                print(f"      -> needs: {c.what_would_change}")
    manual = [c.code for c in report.criteria if c.status == "manual"]
    if manual:
        print(f"\n[MAN] Requires human review: {', '.join(manual)}")

    if report.bundle:
        print(f"\n{thin}")
        print("SOURCES")
        for s in report.bundle.sources:
            state = "ok" if s.ok else "FAILED"
            cached = " (cached)" if s.cached else ""
            detail = f" -- {s.detail}" if s.detail else ""
            print(f"  {s.source:14} {state}{cached}{detail}")

    print(
        "\nDraft worksheet from public data. Not a clinical determination; "
        "criteria marked [MAN] were not assessed."
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="dissent",
        description="Show the argument behind a variant interpretation, not just the answer.",
    )
    parser.add_argument("variant", nargs="?", help="rsID, HGVS, coordinates, or 'GENE p.Change'")
    parser.add_argument("--assembly", default="hg38", choices=["hg19", "hg38"])
    parser.add_argument("--json", action="store_true", help="emit the full report as JSON")
    parser.add_argument("--serve", action="store_true", help="start the web UI instead")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args(argv)

    if args.serve or args.variant == "serve":
        import uvicorn

        print(f"Dissent running at http://{args.host}:{args.port}")
        uvicorn.run("dissent.api.app:app", host=args.host, port=args.port, log_level="info")
        return 0

    if not args.variant:
        parser.print_help()
        return 2

    report = asyncio.run(interpret(args.variant, assembly=args.assembly, cache=ResponseCache()))

    if args.json:
        print(report.model_dump_json(indent=2))
    else:
        _print_report(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
