"""Ensembl REST VEP client.

Supplies the molecular consequence, which decides which ACMG criteria are even
applicable: PVS1 needs a loss-of-function term, PP2/BP1 need missense, BP7
needs synonymous, PM4 needs a protein-length change.

The canonical transcript is preferred. Where VEP returns several transcripts,
picking an arbitrary one changes the consequence and therefore the
classification, so selection is explicit and the chosen transcript is reported.
"""

from __future__ import annotations

from typing import Any

from dissent.models import Consequence, SourceStatus, VariantQuery
from dissent.sources.base import SourceClient, to_int

API = "https://rest.ensembl.org"

# VEP's own severity ordering, most severe first. Used to choose between
# transcripts when none is flagged canonical.
_SEVERITY = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intergenic_variant",
)
_SEVERITY_RANK = {term: i for i, term in enumerate(_SEVERITY)}


def _rank(terms: list[str]) -> int:
    return min((_SEVERITY_RANK.get(t, 999) for t in terms), default=999)


class EnsemblClient(SourceClient):
    name = "ensembl_vep"

    async def fetch(self, query: VariantQuery) -> tuple[Consequence | None, SourceStatus]:
        # Coordinates take priority over the rsID. The rsID endpoint returns
        # every allele colocated at that position, and reading the first gives
        # the wrong amino acid change on a multi-allelic site: rs113488022
        # yielded p.Val600Gly rather than BRAF p.Val600Glu this way.
        if query.chrom and query.pos and query.ref and query.alt:
            region = f"{query.chrom}:{query.pos}-{query.pos + len(query.ref) - 1}"
            url = f"{API}/vep/human/region/{region}/{query.alt}"
            key = ["region", region, query.alt]
        elif query.rsid:
            url = f"{API}/vep/human/id/{query.rsid}"
            key = ["id", query.rsid]
        else:
            return None, SourceStatus(
                source=self.name, ok=True, detail="no rsID or coordinates to query"
            )

        # VEP omits HGVS strings and MANE flags unless asked. Without hgvs=1
        # the report cannot show `c.` / `p.` notation; without mane=1 transcript
        # selection falls back to Ensembl canonical, which is not always MANE.
        params = {
            "content-type": "application/json",
            "hgvs": "1",
            "canonical": "1",
            "mane": "1",
            "numbers": "1",
        }
        data, status = await self.fetch_json(
            url, params=params, cache_key_parts=[*key, "hgvs-mane"]
        )
        if not isinstance(data, list) or not data:
            return None, status if not status.ok else SourceStatus(
                source=self.name, ok=True, detail="no VEP annotation", cached=status.cached
            )

        record = data[0]
        if not isinstance(record, dict):
            return None, status

        consequences = [
            c for c in (record.get("transcript_consequences") or []) if isinstance(c, dict)
        ]
        if not consequences:
            return (
                Consequence(
                    most_severe=record.get("most_severe_consequence"),
                    terms=[record.get("most_severe_consequence")]
                    if record.get("most_severe_consequence")
                    else [],
                ),
                status,
            )

        chosen = self._choose_transcript(consequences, query)
        terms = [t for t in (chosen.get("consequence_terms") or []) if t]

        return (
            Consequence(
                most_severe=record.get("most_severe_consequence") or (terms[0] if terms else None),
                transcript_id=chosen.get("transcript_id"),
                gene_symbol=chosen.get("gene_symbol"),
                hgvs_c=chosen.get("hgvsc"),
                hgvs_p=chosen.get("hgvsp"),
                protein_position=to_int(chosen.get("protein_start")),
                amino_acids=chosen.get("amino_acids"),
                canonical=bool(chosen.get("canonical")),
                biotype=chosen.get("biotype"),
                terms=terms,
            ),
            status,
        )

    @staticmethod
    def _choose_transcript(
        consequences: list[dict[str, Any]], query: VariantQuery
    ) -> dict[str, Any]:
        """Pick one transcript, preferring: requested gene > MANE > canonical > severity."""
        pool = consequences
        if query.gene:
            same_gene = [
                c
                for c in pool
                if str(c.get("gene_symbol", "")).upper() == query.gene.upper()
            ]
            if same_gene:
                pool = same_gene

        protein_coding = [c for c in pool if c.get("biotype") == "protein_coding"]
        if protein_coding:
            pool = protein_coding

        for key in ("mane_select", "mane"):
            mane = [c for c in pool if c.get(key)]
            if mane:
                return mane[0]

        canonical = [c for c in pool if c.get("canonical")]
        if canonical:
            return canonical[0]

        return min(pool, key=lambda c: _rank(c.get("consequence_terms") or []))
