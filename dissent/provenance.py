"""Links back to the records a claim rests on.

These were duplicated across the criteria and the detectors, with the gnomAD
one hard-coding a database build. One copy, so a provenance label cannot drift
away from the data it describes.
"""

from __future__ import annotations

from dissent.models import EvidenceBundle, Provenance

CLINVAR_VARIATION = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/"
GNOMAD_VARIANT = "https://gnomad.broadinstitute.org/variant/{vid}?dataset=gnomad_r4"
GNOMAD_GENE = "https://gnomad.broadinstitute.org/gene/{symbol}?dataset=gnomad_r4"
LITVAR = "https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?text={rsid}"


def clinvar(bundle: EvidenceBundle) -> list[Provenance]:
    cv = bundle.clinvar
    if cv is None or cv.variant_id is None:
        return []
    return [
        Provenance(
            source="ClinVar",
            label=f"ClinVar variation {cv.variant_id}",
            url=CLINVAR_VARIATION.format(vid=cv.variant_id),
        )
    ]


def gnomad_variant(bundle: EvidenceBundle) -> list[Provenance]:
    """Link to the gnomAD record, naming the build the figures came from.

    A frequency may have come from the live v4 API or from the v2 copy embedded
    in the MyVariant document; labelling both "v4" put the wrong database build
    on a claim. An hg19 query links to no gnomAD record at all, since the
    browser is GRCh38 only.
    """
    vid = bundle.query.gnomad_variant_id
    if not vid or bundle.query.assembly != "hg38":
        return []
    used_v4 = any("v4" in f.source for f in bundle.frequencies)
    return [
        Provenance(
            source="gnomAD",
            label=f"gnomAD {'v4' if used_v4 else 'browser'} variant {vid}",
            url=GNOMAD_VARIANT.format(vid=vid),
        )
    ]


def gnomad_gene(bundle: EvidenceBundle) -> list[Provenance]:
    if not bundle.gene_symbol:
        return []
    return [
        Provenance(
            source="gnomAD",
            label=f"gnomAD gene constraint for {bundle.gene_symbol}",
            url=GNOMAD_GENE.format(symbol=bundle.gene_symbol),
        )
    ]


def litvar(bundle: EvidenceBundle) -> list[Provenance]:
    rsid = bundle.query.rsid
    total = bundle.literature_total or len(bundle.literature_pmids)
    if not rsid or not total:
        return []
    return [
        Provenance(
            source="LitVar",
            label=f"{total:,} linked publications",
            url=LITVAR.format(rsid=rsid),
        )
    ]
