"""Fetch every source and assemble one EvidenceBundle.

Fetching runs in two phases because of a real dependency: the gene symbol is
not known until MyVariant or VEP answers, and gnomAD constraint plus Open
Targets are keyed on it.

    phase 1 (parallel): MyVariant | Ensembl VEP | gnomAD variant
    phase 2 (parallel): gnomAD constraint | Open Targets | LitVar

Within a phase, requests are concurrent, so total latency is the slower of the
two phases rather than the sum of six round trips.
"""

from __future__ import annotations

import asyncio
from datetime import datetime, timezone
from typing import Any

import httpx

from dissent.cache.store import ResponseCache
from dissent.models import EvidenceBundle, PopulationFrequency, SourceStatus, VariantQuery
from dissent.sources.base import USER_AGENT, dig, to_float
from dissent.sources.ensembl import EnsemblClient
from dissent.sources.gnomad import GnomadClient
from dissent.sources.litvar import LitVarClient
from dissent.sources.myvariant import (
    MyVariantClient,
    extract_clinvar,
    extract_frequencies,
    extract_predictors,
)
from dissent.normalize.consequence import from_snpeff, transcript_from_preferred_name
from dissent.sources.opentargets import OpenTargetsClient


# How many PMIDs to keep. The report shows a count and a handful of links.
PMID_LIMIT = 200


def _spliceai_max(doc: dict[str, Any]) -> float | None:
    """Highest SpliceAI delta score across the four acceptor/donor predictions.

    SpliceAI reports gain/loss for both acceptor and donor separately; the
    convention for interpretation is to use the maximum.
    """
    node = doc.get("spliceai")
    if node is None:
        return None

    # MyVariant returns a list when several genes overlap the position, and a
    # bare object otherwise. Requiring a dict silently produced "no SpliceAI
    # score was available" for exactly the overlapping-gene cases, leaving BP7
    # unresolved on a wrong premise.
    entries = node if isinstance(node, list) else [node]

    candidates: list[float] = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        # Some builds nest the deltas under `masked` / `unmasked`.
        for container in (entry, entry.get("masked"), entry.get("unmasked")):
            if not isinstance(container, dict):
                continue
            for key in ("ds_ag", "ds_al", "ds_dg", "ds_dl"):
                value = to_float(container.get(key))
                if value is not None:
                    candidates.append(value)
    return max(candidates) if candidates else None


def _merge_frequencies(
    primary: list[PopulationFrequency], secondary: list[PopulationFrequency]
) -> list[PopulationFrequency]:
    """gnomAD's own API wins over MyVariant's embedded copy.

    MyVariant currently carries gnomAD v2 exomes/genomes; the live API serves
    v4, which has roughly an order of magnitude more samples. Keeping both
    would double-count the same cohort when computing a maximum, so v4 records
    replace the embedded ones and the embedded set is only a fallback.
    """
    if primary:
        return primary
    return secondary


def _resolve_coordinates(
    query: VariantQuery, doc: dict[str, Any], clinvar: Any
) -> VariantQuery:
    """Fill in rsID, chrom, pos, ref and alt from whatever the lookup returned.

    Downstream sources are keyed on these: gnomAD needs `chrom-pos-ref-alt` and
    VEP needs an rsID or a region. A query entered as `BRAF V600E` has none of
    them until MyVariant answers.
    """
    resolved = query.model_copy()

    if not resolved.rsid:
        resolved.rsid = (clinvar.rsid if clinvar else None) or _dbsnp_rsid(doc)
    if not resolved.chrom and doc.get("chrom"):
        resolved.chrom = str(doc["chrom"])

    # `vcf` is the authoritative allele record in a MyVariant document; the
    # top level has no ref/alt of its own.
    vcf = doc.get("vcf")
    if isinstance(vcf, dict):
        if resolved.ref is None and vcf.get("ref"):
            resolved.ref = str(vcf["ref"])
        if resolved.alt is None and vcf.get("alt"):
            resolved.alt = str(vcf["alt"])
        if resolved.pos is None and vcf.get("position"):
            try:
                resolved.pos = int(vcf["position"])
            except (TypeError, ValueError):
                pass

    if resolved.pos is None:
        # Take a position only from the block matching the declared assembly.
        # Falling through from hg38 to hg19 produced an hg19 coordinate still
        # labelled hg38, which was then queried against GRCh38-only gnomAD and
        # written into a mislabelled `hgvs_g`.
        node = doc.get("hg38" if resolved.assembly == "hg38" else "hg19")
        if isinstance(node, dict) and node.get("start"):
            resolved.pos = node["start"]

    if not resolved.hgvs_g and resolved.chrom and resolved.pos and resolved.ref and resolved.alt:
        resolved.hgvs_g = (
            f"chr{resolved.chrom}:g.{resolved.pos}{resolved.ref}>{resolved.alt}"
        )
    return resolved


def _dbsnp_rsid(doc: dict[str, Any]) -> str | None:
    rsid = dig(doc, "dbsnp", "rsid")
    if rsid:
        text = str(rsid)
        return text if text.startswith("rs") else f"rs{text}"
    return None


def _dbsnp_citations(doc: dict[str, Any]) -> list[int]:
    out: list[int] = []
    for value in (dig(doc, "dbsnp", "citations") or []):
        try:
            out.append(int(value))
        except (TypeError, ValueError):
            continue
    return out


def _dedupe_statuses(statuses: list[SourceStatus]) -> list[SourceStatus]:
    """One row per source, preferring the attempt that actually returned data.

    A source retried after coordinate resolution appears twice: once saying it
    had nothing to query, once with the real result. Showing both in the UI
    would read as a failure.
    """
    best: dict[str, SourceStatus] = {}
    for status in statuses:
        existing = best.get(status.source)
        if existing is None:
            best[status.source] = status
            continue
        existing_blank = bool(existing.detail and "no " in existing.detail.lower())
        if existing_blank and not status.detail:
            best[status.source] = status
        elif not existing.ok and status.ok:
            best[status.source] = status
    return list(best.values())


async def collect_evidence(
    query: VariantQuery,
    cache: ResponseCache | None = None,
    client: httpx.AsyncClient | None = None,
) -> EvidenceBundle:
    owns_client = client is None
    if client is None:
        client = httpx.AsyncClient(
            timeout=20.0, headers={"User-Agent": USER_AGENT}, follow_redirects=True
        )

    statuses: list[SourceStatus] = []
    try:
        # Budgets are tuned per source from measured behaviour. gnomAD's
        # variant endpoint is markedly slow for common variants, so it gets a
        # single attempt: its frequencies are also carried inside the MyVariant
        # document, which `_merge_frequencies` falls back to.
        myvariant = MyVariantClient(client=client, cache=cache, timeout=12.0, budget=26.0)
        # Ensembl VEP is the slowest source by an order of magnitude and the
        # least reliable. Sampled over ten variants: median 4.1s, but 18s and
        # 28s in the tail and two that never answered inside 30s. At the old
        # 25s timeout one of those tail cases held the whole page, which took
        # a cold lookup to a median of 25s.
        #
        # 8s rather than something longer, because raising the cap does not
        # buy availability: measured over the same eight variants, 8s and 12s
        # both left VEP unavailable on five of them, so the slow ones are
        # failing rather than being slow, and the extra four seconds are spent
        # waiting for an answer that does not arrive. No room to retry for the
        # same reason. Cold median with this cap: 9.7s.
        #
        # Losing it costs the second opinion on the molecular consequence, not
        # the worksheet -- snpEff supplies the consequence through MyVariant,
        # so the classification, the criteria and the other detectors are
        # unchanged, the annotator-disagreement check simply has nothing to
        # compare, and the source panel says VEP did not answer.
        ensembl = EnsemblClient(client=client, cache=cache, timeout=8.0, budget=9.0)
        gnomad = GnomadClient(client=client, cache=cache, timeout=10.0, budget=11.0)
        opentargets = OpenTargetsClient(client=client, cache=cache, timeout=10.0, budget=21.0)
        litvar = LitVarClient(client=client, cache=cache, timeout=12.0, budget=13.0)

        # ---- phase 1 -----------------------------------------------------
        # VEP is queried by region+allele, never by rsID, whenever coordinates
        # are available. An rsID names a position, not an allele: querying VEP
        # with a multi-allelic rsID returns every colocated allele and taking
        # the first silently reports the wrong amino acid change. Observed with
        # rs113488022, where it yielded p.Val600Gly instead of p.Val600Glu.
        can_query_vep_now = bool(query.chrom and query.pos and query.ref and query.alt)

        phase1 = [myvariant.fetch(query), gnomad.fetch_variant(query)]
        if can_query_vep_now:
            phase1.append(ensembl.fetch(query))

        phase1_results = await asyncio.gather(*phase1, return_exceptions=True)
        mv_result = phase1_results[0]
        gnomad_freq_result = phase1_results[1]
        vep_result = phase1_results[2] if can_query_vep_now else None

        doc: dict[str, Any] | None = None
        if isinstance(mv_result, tuple):
            doc, status = mv_result
            statuses.append(status)
        else:
            statuses.append(
                SourceStatus(source="myvariant", ok=False, detail=str(mv_result))
            )

        consequence = None
        if isinstance(vep_result, tuple):
            consequence, status = vep_result
            statuses.append(status)
        elif vep_result is not None:
            statuses.append(
                SourceStatus(source="ensembl_vep", ok=False, detail=str(vep_result))
            )

        gnomad_freqs: list[PopulationFrequency] = []
        if isinstance(gnomad_freq_result, tuple):
            gnomad_freqs, status = gnomad_freq_result
            statuses.append(status)
        else:
            statuses.append(SourceStatus(source="gnomad", ok=False, detail=str(gnomad_freq_result)))

        doc = doc or {}
        clinvar = extract_clinvar(doc)
        predictors = extract_predictors(doc)
        embedded_freqs = extract_frequencies(doc)

        # Primary consequence: snpEff, already present in `doc`, numbered
        # against the transcript ClinVar used where that is known.
        preferred_tx = transcript_from_preferred_name(
            clinvar.preferred_name if clinvar else None
        )
        snpeff_consequence = from_snpeff(doc, preferred_tx)
        vep_consequence = consequence
        consequence = snpeff_consequence or consequence

        # Fill in coordinates discovered downstream, so phase 2 and any later
        # re-query can use them.
        resolved = _resolve_coordinates(query, doc, clinvar)

        gene_symbol = (
            (clinvar.gene if clinvar else None)
            or (consequence.gene_symbol if consequence else None)
            or resolved.gene
            or dig(doc, "dbnsfp", "genename")
        )
        if isinstance(gene_symbol, list):
            gene_symbol = gene_symbol[0] if gene_symbol else None

        # Phase 1 runs before the variant is resolved, so a query given as
        # `GENE p.X` reaches gnomAD and VEP with nothing to look up. Now that
        # MyVariant has supplied an rsID and coordinates, retry whichever of
        # them came back empty. Retried concurrently, not in sequence.
        retry_gnomad = (
            not gnomad_freqs
            and bool(resolved.gnomad_variant_id)
            and resolved.assembly == "hg38"
        )
        has_coords = bool(resolved.chrom and resolved.pos and resolved.ref and resolved.alt)
        # Retry VEP whenever it has not answered yet: it is the cross-check,
        # so it is worth fetching even when snpEff already gave a consequence.
        retry_vep = vep_consequence is None and (has_coords or bool(resolved.rsid))

        # Build only the coroutines that will actually be awaited. Creating one
        # unconditionally and dropping it raises "coroutine was never awaited",
        # which printed a Python warning above the CLI's own output.
        retries: list[tuple[str, Any]] = []
        if retry_gnomad:
            retries.append(("gnomad", gnomad.fetch_variant(resolved)))
        if retry_vep:
            retries.append(("vep", ensembl.fetch(resolved)))

        if retries:
            outcomes = await asyncio.gather(
                *(coro for _, coro in retries), return_exceptions=True
            )
            for (kind, _), outcome in zip(retries, outcomes):
                if not isinstance(outcome, tuple):
                    continue
                value, status = outcome
                if kind == "gnomad":
                    if value:
                        gnomad_freqs = value
                elif value is not None:
                    vep_consequence = value
                    # snpEff may have had nothing to say (non-coding, say).
                    if consequence is None:
                        consequence = value
                statuses.append(status)

        # ---- phase 2 -----------------------------------------------------
        constraint_result, target_result, litvar_result = await asyncio.gather(
            gnomad.fetch_constraint(gene_symbol or ""),
            opentargets.fetch_gene(gene_symbol or ""),
            litvar.fetch_pmids(resolved.rsid),
            return_exceptions=True,
        )

        constraint = None
        if isinstance(constraint_result, tuple):
            constraint, status = constraint_result
            statuses.append(status)

        gene_description = None
        if isinstance(target_result, tuple):
            target, status = target_result
            statuses.append(status)
            if target:
                gene_description = target.get("function")

        pmids: list[int] = []
        literature_total = 0
        if isinstance(litvar_result, tuple):
            pmids, literature_total, status = litvar_result
            statuses.append(status)
        if not pmids:
            # dbSNP ships its own citation list, which costs no extra request
            # because it arrives inside the MyVariant document.
            pmids = _dbsnp_citations(doc)
            literature_total = len(pmids)
        pmids = pmids[:PMID_LIMIT]

        return EvidenceBundle(
            query=resolved,
            gene_symbol=gene_symbol,
            clinvar=clinvar,
            predictors=predictors,
            frequencies=_merge_frequencies(gnomad_freqs, embedded_freqs),
            constraint=constraint,
            consequence=consequence,
            consequence_vep=vep_consequence,
            literature_pmids=pmids,
            literature_total=literature_total,
            gene_description=gene_description,
            alternative_alleles=[str(x) for x in (doc.get("_siblings") or [])],
            splice_ai_max=_spliceai_max(doc),
            sources=_dedupe_statuses(statuses),
        )
    finally:
        if owns_client:
            await client.aclose()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")
