"""gnomAD GraphQL client (v4).

Supplies two things MyVariant's embedded frequencies do not:

* **Homozygote counts**, needed for BS2 (observed in a healthy adult in the
  homozygous state).
* **Gene-level constraint** (pLI, LOEUF, mis_z), needed for PVS1's "loss of
  function is a known mechanism of disease" test and for PP2/BP1.

gnomAD v4 is GRCh38-only. A query built from hg19 coordinates would silently
return the wrong locus, so `fetch_variant` refuses hg19 input rather than
guessing a liftover.
"""

from __future__ import annotations

from typing import Any

from dissent.models import GeneConstraint, PopulationFrequency, SourceStatus, VariantQuery
from dissent.sources.base import SourceClient, to_float, to_int

API = "https://gnomad.broadinstitute.org/api"
VARIANT_URL = "https://gnomad.broadinstitute.org/variant/{vid}?dataset=gnomad_r4"
GENE_URL = "https://gnomad.broadinstitute.org/gene/{sym}?dataset=gnomad_r4"

_VARIANT_QUERY = """
query Variant($variantId: String!) {
  variant(variantId: $variantId, dataset: gnomad_r4) {
    variant_id
    rsids
    exome { ac an af homozygote_count populations { id ac an } }
    genome { ac an af homozygote_count populations { id ac an } }
  }
}
"""

_GENE_QUERY = """
query Gene($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id
    symbol
    gnomad_constraint { pLI oe_lof oe_lof_upper oe_mis mis_z lof_z }
  }
}
"""

# gnomAD population ids -> readable labels. Sub-populations are excluded: the
# ACMG frequency thresholds are defined on continental groupings.
_POPS = {
    "afr": "African/African-American",
    "amr": "Admixed American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "nfe": "Non-Finnish European",
    "sas": "South Asian",
    "mid": "Middle Eastern",
    "remaining": "Remaining",
}


def _frequency_from(node: dict[str, Any] | None, label: str) -> PopulationFrequency | None:
    if not isinstance(node, dict):
        return None
    an = to_int(node.get("an"))
    if not an:
        return None

    per_pop: dict[str, float] = {}
    for pop in node.get("populations") or []:
        if not isinstance(pop, dict):
            continue
        pid = str(pop.get("id", "")).lower()
        if pid not in _POPS:
            continue
        pop_an, pop_ac = to_int(pop.get("an")), to_int(pop.get("ac"))
        # Small denominators give wild frequency estimates. gnomAD's own
        # guidance for BA1/BS1 is to ignore populations below ~2000 alleles.
        if pop_an and pop_an >= 2000 and pop_ac is not None:
            per_pop[_POPS[pid]] = pop_ac / pop_an

    popmax_label, popmax_af = (None, None)
    if per_pop:
        popmax_label, popmax_af = max(per_pop.items(), key=lambda kv: kv[1])

    return PopulationFrequency(
        source=label,
        af=to_float(node.get("af")),
        ac=to_int(node.get("ac")),
        an=an,
        homozygote_count=to_int(node.get("homozygote_count")),
        popmax_af=popmax_af,
        popmax_population=popmax_label,
        per_population=per_pop,
    )


class GnomadClient(SourceClient):
    name = "gnomad"

    async def fetch_variant(
        self, query: VariantQuery
    ) -> tuple[list[PopulationFrequency], SourceStatus]:
        vid = query.gnomad_variant_id
        if not vid:
            return [], SourceStatus(
                source=self.name, ok=True, detail="no genomic coordinates to query"
            )
        if query.assembly != "hg38":
            return [], SourceStatus(
                source=self.name,
                ok=True,
                detail="gnomAD v4 is GRCh38-only; hg19 coordinates not queried",
            )

        data, status = await self.fetch_json(
            API,
            method="POST",
            json_body={"query": _VARIANT_QUERY, "variables": {"variantId": vid}},
            cache_key_parts=["variant", vid],
        )
        if not isinstance(data, dict):
            return [], status

        variant = (data.get("data") or {}).get("variant")
        if not isinstance(variant, dict):
            return [], SourceStatus(
                source=self.name, ok=True, detail="not in gnomAD v4", cached=status.cached
            )

        out = []
        for key, label in (("exome", "gnomAD v4 exomes"), ("genome", "gnomAD v4 genomes")):
            freq = _frequency_from(variant.get(key), label)
            if freq:
                out.append(freq)
        return out, status

    async def fetch_constraint(
        self, gene_symbol: str
    ) -> tuple[GeneConstraint | None, SourceStatus]:
        if not gene_symbol:
            return None, SourceStatus(source=self.name, ok=True, detail="no gene symbol")

        data, status = await self.fetch_json(
            API,
            method="POST",
            json_body={"query": _GENE_QUERY, "variables": {"symbol": gene_symbol}},
            cache_key_parts=["gene", gene_symbol],
        )
        if not isinstance(data, dict):
            return None, status

        gene = (data.get("data") or {}).get("gene")
        if not isinstance(gene, dict):
            return None, SourceStatus(
                source=self.name, ok=True, detail="gene not found", cached=status.cached
            )

        c = gene.get("gnomad_constraint") or {}
        return (
            GeneConstraint(
                gene_symbol=gene.get("symbol"),
                gene_id=gene.get("gene_id"),
                pli=to_float(c.get("pLI")),
                oe_lof=to_float(c.get("oe_lof")),
                oe_lof_upper=to_float(c.get("oe_lof_upper")),
                oe_mis=to_float(c.get("oe_mis")),
                mis_z=to_float(c.get("mis_z")),
                lof_z=to_float(c.get("lof_z")),
            ),
            status,
        )
