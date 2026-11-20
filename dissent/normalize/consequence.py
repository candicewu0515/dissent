"""Build a Consequence from the snpEff annotations inside a MyVariant document.

Why not VEP as the primary source: `rest.ensembl.org` latency is erratic,
measured between 3 and 60 seconds for the same request, and a report that
sometimes takes a minute is not usable. snpEff annotations arrive inside the
MyVariant document that has already been fetched, at no extra cost, and carry
the same fields the criteria need.

VEP is still queried, as an independent cross-check. Where the two disagree on
the molecular consequence, that disagreement is reported rather than resolved:
which of them is right determines whether PVS1 or PP2 applies at all.

Transcript choice is the hard part. For BRCA1 c.181T>G, snpEff reports
p.Cys61Gly on NM_007300.3 and p.Cys14Gly on NM_007297.3 — the same variant,
numbered against different transcripts. ClinVar's `preferred_name` names the
MANE transcript it used, so that is the anchor when available.
"""

from __future__ import annotations

import re
from typing import Any

from dissent.models import Consequence
from dissent.sources.base import as_list, to_int

# snpEff's putative_impact ordering, most severe first.
_IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

# Anchored at the start of the string on purpose. snpEff also emits
# structural annotations whose feature_id embeds a transcript late in the
# string, e.g. `1GZH:A_175-A_237:NM_000546.5`. An unanchored search matched
# those and selected a PDB interaction record as the primary consequence for
# TP53 p.Arg175His.
_TX_ACCESSION_ANCHORED = re.compile(r"^((?:NM|NR|NP|ENST)_?\d+)", re.I)
# For ClinVar's preferred_name the accession is not at position 0, so that one
# is searched rather than anchored.
_TX_ACCESSION_SEARCH = re.compile(r"\b((?:NM|NR|NP|ENST)_?\d+)", re.I)

# snpEff effects that are not transcript-level molecular consequences. They
# must not be selected as the variant's consequence.
_NON_CONSEQUENCE_EFFECTS = frozenset(
    {
        "structural_interaction_variant",
        "protein_protein_contact",
        "sequence_feature",
    }
)


def transcript_from_preferred_name(preferred_name: str | None) -> str | None:
    """`NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)` -> `NM_007294`."""
    if not preferred_name:
        return None
    match = _TX_ACCESSION_SEARCH.search(preferred_name)
    return match.group(1).upper() if match else None


def _bare_accession(feature_id: Any) -> str | None:
    if not feature_id:
        return None
    match = _TX_ACCESSION_ANCHORED.match(str(feature_id).strip())
    return match.group(1).upper() if match else None


def from_snpeff(doc: dict[str, Any], preferred_transcript: str | None = None) -> Consequence | None:
    """Pick one snpEff annotation and convert it to a Consequence."""
    annotations = [
        a
        for a in as_list((doc.get("snpeff") or {}).get("ann") if isinstance(doc.get("snpeff"), dict) else None)
        if isinstance(a, dict)
    ]
    if not annotations:
        return None

    # Keep only transcript-level consequence records.
    usable = [
        a
        for a in annotations
        if a.get("feature_type") == "transcript"
        and not (_NON_CONSEQUENCE_EFFECTS & set(str(a.get("effect") or "").split("&")))
    ]
    if not usable:
        return None

    coding = [a for a in usable if a.get("transcript_biotype") == "protein_coding"]
    pool = coding or usable

    chosen: dict[str, Any] | None = None

    # 1. The transcript ClinVar numbers against, when we know it.
    if preferred_transcript:
        want = preferred_transcript.upper()
        for annotation in pool:
            if _bare_accession(annotation.get("feature_id")) == want:
                chosen = annotation
                break

    # 2. Otherwise most severe impact, then the longest protein, which is a
    #    reasonable stand-in for the principal isoform.
    if chosen is None:
        def rank(a: dict[str, Any]) -> tuple[int, int]:
            impact = _IMPACT_RANK.get(str(a.get("putative_impact", "")).upper(), 4)
            length = to_int((a.get("protein") or {}).get("length")) or 0
            return (impact, -length)

        chosen = min(pool, key=rank)

    effect = str(chosen.get("effect") or "")
    # snpEff joins compound effects with `&`, e.g. `splice_region_variant&intron_variant`.
    terms = [t for t in effect.split("&") if t]

    return Consequence(
        most_severe=terms[0] if terms else None,
        transcript_id=chosen.get("feature_id"),
        gene_symbol=chosen.get("genename") or chosen.get("gene_id"),
        hgvs_c=chosen.get("hgvs_c"),
        hgvs_p=chosen.get("hgvs_p"),
        protein_position=to_int((chosen.get("protein") or {}).get("position")),
        canonical=bool(preferred_transcript)
        and _bare_accession(chosen.get("feature_id")) == preferred_transcript.upper(),
        biotype=chosen.get("transcript_biotype"),
        terms=terms,
    )


def disagreement(a: Consequence | None, b: Consequence | None) -> str | None:
    """Describe a consequence-level disagreement between two annotators.

    Only differences that change which criteria apply are reported. A variant
    called `missense_variant` by one tool and `splice_region_variant` by the
    other is a materially different interpretation; differing transcript IDs
    with the same effect are not.
    """
    if a is None or b is None:
        return None

    def bucket(c: Consequence) -> str:
        if c.is_lof:
            return "loss of function"
        if c.is_missense:
            return "missense"
        if c.is_synonymous:
            return "synonymous"
        if c.is_inframe_indel:
            return "in-frame indel"
        return c.most_severe or "other"

    ba, bb = bucket(a), bucket(b)
    if ba == bb:
        return None
    return f"{ba} vs {bb}"
