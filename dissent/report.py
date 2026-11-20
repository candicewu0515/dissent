"""Assemble the full report: evidence -> criteria -> classification -> conflicts."""

from __future__ import annotations

from datetime import date

import httpx

from dissent.acmg.combine import classify
from dissent.acmg.criteria import evaluate_all
from dissent.acmg.thresholds import DEFAULT, Thresholds
from dissent.cache.store import ResponseCache
from dissent.collect import collect_evidence, utc_now_iso
from dissent.conflicts.detectors import detect_all
from dissent.models import EvidenceBundle, Report, VariantQuery
from dissent.normalize.query import parse_variant
from dissent.summary import plain_summary


def build_report(
    bundle: EvidenceBundle,
    thresholds: Thresholds = DEFAULT,
    today: date | None = None,
) -> Report:
    """Pure function: evidence in, report out. No I/O, so it is directly testable."""
    criteria = evaluate_all(bundle, thresholds)
    classification = classify(criteria)
    conflicts = detect_all(
        bundle, criteria=criteria, classification=classification, t=thresholds, today=today
    )

    display = None
    if bundle.clinvar and bundle.clinvar.preferred_name:
        display = bundle.clinvar.preferred_name
    elif bundle.consequence and bundle.consequence.hgvs_p:
        display = f"{bundle.gene_symbol} {bundle.consequence.hgvs_p}"
    elif bundle.query.hgvs_g:
        display = bundle.query.hgvs_g

    return Report(
        query=bundle.query,
        gene_symbol=bundle.gene_symbol,
        display_name=display,
        summary=plain_summary(bundle, classification, conflicts),
        criteria=criteria,
        classification=classification,
        conflicts=conflicts,
        bundle=bundle,
        generated_at=utc_now_iso(),
    )


async def interpret(
    raw: str,
    assembly: str = "hg38",
    cache: ResponseCache | None = None,
    client: httpx.AsyncClient | None = None,
    thresholds: Thresholds = DEFAULT,
) -> Report:
    """Parse, fetch, and build a report for a free-form variant string."""
    query: VariantQuery = parse_variant(raw, assembly=assembly)  # type: ignore[arg-type]
    bundle = await collect_evidence(query, cache=cache, client=client)
    return build_report(bundle, thresholds=thresholds)
