"""Command line entry point.

    python -m dissent serve                 start the web UI
    python -m dissent "BRAF V600E"          print a report to the terminal
    python -m dissent --panel list.txt      triage a panel, ordered worklist
    python -m dissent --panel list.txt --csv > out.csv
"""

from __future__ import annotations

import argparse
import asyncio
import json
import sys

from dissent.cache.store import ResponseCache
from dissent.report import interpret

_SEVERITY_MARK = {"high": "!!", "medium": "! ", "low": "  "}
_STATUS_MARK = {
    "met": "[MET]",
    "not_met": "[ - ]",
    "not_applicable": "[ / ]",
    "insufficient": "[ ? ]",
    "manual": "[MAN]",
}


def _resolved(report) -> bool:
    """Whether the lookup actually found the variant.

    Coordinates alone are enough: a variant can legitimately be absent from
    ClinVar and still be a real, annotatable locus.
    """
    q = report.query
    if q.chrom and q.pos:
        return True
    bundle = report.bundle
    return bool(bundle and (bundle.clinvar or bundle.consequence))


def _wrap(text: str, width: int) -> list[str]:
    import textwrap

    return textwrap.wrap(text, width=width) or [""]


def _print_report(report, width: int = 80) -> None:
    rule = "=" * width
    thin = "-" * width

    print(rule)
    print(report.display_name or report.query.raw)
    q = report.query
    coords = f"{q.chrom}:{q.pos}{q.ref}>{q.alt}" if q.chrom and q.pos else "(unresolved)"
    print(f"{report.gene_symbol or '?'}   {coords}   {q.rsid or ''}   [{q.assembly}]")
    print(rule)

    if report.summary:
        print()
        for line in _wrap(report.summary, width - 2):
            print(f" {line}")

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
            print(f"{_STATUS_MARK[c.status]} {c.code}{strength}  {c.plain or c.name}")
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
    parser.add_argument(
        "--panel",
        metavar="FILE",
        help="a file of variants, one per line; prints a triaged worklist (- for stdin)",
    )
    parser.add_argument("--csv", action="store_true", help="with --panel, emit CSV")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args(argv)

    if args.serve or args.variant == "serve":
        import uvicorn

        print(f"Dissent running at http://{args.host}:{args.port}")
        uvicorn.run("dissent.api.app:app", host=args.host, port=args.port, log_level="info")
        return 0

    if args.panel:
        return _run_panel(args)

    if not args.variant:
        parser.print_help()
        return 2

    report = asyncio.run(interpret(args.variant, assembly=args.assembly, cache=ResponseCache()))

    # A variant that did not resolve must not be presented as a finished
    # worksheet. `dissent banana` previously printed all 28 criteria reading
    # "Uncertain significance / Pathogenic evidence: nothing" and exited 0.
    if not _resolved(report):
        print(
            f"Could not resolve {args.variant!r}.\n\n"
            "Accepted forms:\n"
            "  rs113488022                 dbSNP identifier\n"
            "  BRAF V600E                  gene and protein change\n"
            "  chr7:140753336A>T           genomic coordinates\n"
            "  NM_004333.6:c.1799T>A       transcript and coding change\n\n"
            "Protein-notation indels such as 'CFTR F508del' are not yet parsed; "
            "use the rsID or coordinates for those.\n\n"
            "If the protein change came from clinical literature, note that many "
            "genes are still written in legacy numbering that omits the signal "
            "peptide, while HGVS counts from the initiator methionine. Factor V "
            "Leiden is R506Q in the literature and p.Arg534Gln in HGVS; sickle "
            "haemoglobin is E6V and p.Glu7Val. The two do not resolve to each "
            "other, so search the rsID or the c. notation instead of guessing "
            "the offset.",
            file=sys.stderr,
        )
        return 3

    if args.json:
        print(report.model_dump_json(indent=2))
    else:
        _print_report(report)
    return 0


def _run_panel(args) -> int:
    """Triage a file of variants into an ordered worklist."""
    from dissent.panel import TIERS_BY_KEY, interpret_panel, parse_panel, rows_to_csv

    raw = sys.stdin.read() if args.panel == "-" else open(args.panel, encoding="utf-8").read()
    queries = parse_panel(raw)
    if not queries:
        print(f"No variants found in {args.panel!r}.", file=sys.stderr)
        return 2

    result = asyncio.run(
        interpret_panel(queries, assembly=args.assembly, cache=ResponseCache())
    )

    if args.csv:
        sys.stdout.write(rows_to_csv(result.rows))
        return 0
    if args.json:
        print(result.model_dump_json(indent=2))
        return 0

    width = 92
    print("=" * width)
    print(f"{len(result.rows)} variants, {result.elapsed_ms / 1000:.1f}s")
    for tier in TIERS_BY_KEY.values():
        n = result.tier_counts.get(tier.key, 0)
        if n:
            print(f"  {n:3d}  {tier.name}")
    if result.evidence_not_public:
        print(
            f"\n  {result.evidence_not_public} reported confidently while the computed "
            "call stays uncertain.\n  Expected: a laboratory's segregation and "
            "functional evidence is not in a public API."
        )

    last = None
    for row in result.rows:
        if row.tier != last:
            last = row.tier
            tier = TIERS_BY_KEY.get(row.tier)
            print()
            print("-" * width)
            print((tier.name if tier else row.tier).upper())
            print("-" * width)
        stars = "" if row.best_stars is None else f"{row.best_stars}*"
        print(f"  {row.display_name or row.query}")
        print(
            f"      computed {row.computed_call or '-'}   "
            f"reported {row.reported_call or 'no record'} {stars}   "
            f"conflicts {row.conflicts}"
        )
        for reason in row.reasons:
            for line in _wrap(reason, width - 10):
                print(f"      - {line}")
        if row.evidence_not_public:
            print("      - submitter's evidence is not public (expected)")
        if row.error:
            print(f"      ! {row.error}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
