"""A one-sentence reading of the report, for someone with no genetics background.

The worksheet and the conflict panel are written for a curator. This is the
line above them: what the evidence actually says, in ordinary words, before any
vocabulary is needed. It is generated from the same evidence everything else
rests on, never written by hand per variant.
"""

from __future__ import annotations

from dissent.models import Classification, Conflict, EvidenceBundle


def _one_in(af: float) -> str:
    """Carrier frequency as '1 in N people', from an allele frequency.

    Carriers are heterozygotes, 2p(1-p), not 1/p. Quoting the reciprocal of the
    allele frequency would understate how common a carrier is by roughly half.
    """
    carrier = 2 * af * (1 - af)
    if carrier <= 0:
        return "almost nobody"
    n = round(1 / carrier)
    if n <= 1:
        return "most people"
    return f"about 1 in {n:,} people"


def _records_calling_pathogenic(bundle: EvidenceBundle) -> int:
    cv = bundle.clinvar
    if cv is None:
        return 0
    return sum(1 for a in cv.actionable_assertions if a.significance.is_pathogenic_side)


def plain_summary(
    bundle: EvidenceBundle,
    classification: Classification | None,
    conflicts: list[Conflict],
) -> str | None:
    """The headline sentence, chosen by which conflict matters most."""
    top = conflicts[0] if conflicts else None
    cv = bundle.clinvar
    record = bundle.max_af_record()

    if top is None:
        if cv and cv.actionable_assertions:
            labels = {a.significance.value.replace("_", " ") for a in cv.actionable_assertions}
            if len(labels) == 1:
                return (
                    f"Every laboratory that has looked at this variant reached the same "
                    f"conclusion: {next(iter(labels))}. The evidence here does not "
                    "contradict itself."
                )
        return "Nothing in the available evidence contradicts itself."

    if top.kind == "frequency_versus_claim" and record:
        af, _, population = record
        n_path = _records_calling_pathogenic(bundle)
        where = population or "the general population"
        return (
            f"{n_path} laboratories call this variant disease-causing. It is carried by "
            f"{_one_in(af)} of {where} descent. Those two facts are hard to reconcile, "
            "and nothing in either database mentions the other."
        )

    if top.kind == "record_discordance":
        return (
            "Laboratories that have examined this variant did not reach the same "
            "conclusion, and the single label the database shows you is a summary that "
            "hides the disagreement."
        )

    if top.kind == "internal_contradiction":
        return (
            "The evidence pulls both ways: how common this variant is argues it is "
            "harmless, while prediction software argues it is damaging. Only one of "
            "those is a measurement in actual people."
        )

    if top.kind == "computed_versus_reported":
        if classification and classification.call.lower().startswith("uncertain"):
            return (
                "The database is confident about this variant and the public evidence "
                "alone is not enough to reach that conclusion. The gap is a measure of "
                "how much of the reasoning is not published anywhere."
            )
        return (
            "This worksheet and the database reach opposite conclusions from the same "
            "public evidence, which means one of them is using something the other "
            "cannot see."
        )

    if top.kind == "somatic_germline_misuse":
        return (
            "This variant is recorded as arising in tumours rather than being "
            "inherited. The criteria below are for inherited disease, so applying them "
            "here asks the wrong question of the data."
        )

    if top.kind == "weak_evidence_base":
        return (
            "The conclusions on record for this variant come without any stated method, "
            "yet they are displayed with the same weight as fully reviewed ones."
        )

    if top.kind == "stale_assertion":
        return (
            "Some conclusions on record for this variant have not been revisited in "
            "years, and still count toward the label the database shows."
        )

    if top.kind == "predictor_split":
        return (
            "Prediction software does not agree about this variant. Most tools show you "
            "only the majority."
        )

    return top.headline
