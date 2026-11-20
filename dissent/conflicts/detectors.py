"""Detect and rank disagreement in the evidence.

This is the part of the report that existing aggregators do not produce. They
reconcile sources into one label; the reconciliation is where the judgement
lives, and hiding it means the reader cannot audit it.

Each detector is a pure function over the evidence bundle (plus, for two of
them, the criteria already evaluated) and returns zero or one Conflict. They
are deliberately conservative: a detector that fires on every variant trains
the reader to ignore the panel.
"""

from __future__ import annotations

from collections import Counter
from datetime import date
from typing import Callable

from dissent.acmg.combine import (
    BENIGN,
    CONFLICTING,
    LIKELY_BENIGN,
    LIKELY_PATHOGENIC,
    PATHOGENIC,
)
from dissent.acmg.thresholds import DEFAULT, Thresholds
from dissent.models import (
    Classification,
    Conflict,
    CriterionResult,
    EvidenceBundle,
    Provenance,
    Significance,
)
from dissent.normalize.consequence import disagreement as consequence_disagreement

Detector = Callable[..., Conflict | None]

_SEVERITY_ORDER = {"high": 0, "medium": 1, "low": 2}


def _clinvar_link(bundle: EvidenceBundle) -> list[Provenance]:
    cv = bundle.clinvar
    if cv is None or cv.variant_id is None:
        return []
    return [
        Provenance(
            source="ClinVar",
            label=f"ClinVar variation {cv.variant_id}",
            url=f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv.variant_id}/",
        )
    ]


def _gnomad_link(bundle: EvidenceBundle) -> list[Provenance]:
    vid = bundle.query.gnomad_variant_id
    if not vid:
        return []
    return [
        Provenance(
            source="gnomAD",
            label=f"gnomAD v4 variant {vid}",
            url=f"https://gnomad.broadinstitute.org/variant/{vid}?dataset=gnomad_r4",
        )
    ]


# ==========================================================================
# 1. Submitters disagree with each other
# ==========================================================================


def submitter_discordance(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    cv = bundle.clinvar
    if cv is None:
        return None
    assertions = cv.actionable_assertions
    if len(assertions) < 2:
        return None

    counts = Counter(a.significance for a in assertions)
    has_path = any(s.is_pathogenic_side for s in counts)
    has_benign = any(s.is_benign_side for s in counts)
    has_vus = Significance.UNCERTAIN in counts
    flagged_conflicting = Significance.CONFLICTING in counts

    if not (flagged_conflicting or (has_path and has_benign) or (has_path and has_vus)):
        return None

    spread = ", ".join(
        f"{n}x {sig.value.replace('_', ' ')}" for sig, n in counts.most_common()
    )

    if has_path and has_benign:
        severity = "high"
        headline = "Submitters disagree across the pathogenic/benign divide"
    elif flagged_conflicting:
        severity = "high"
        headline = "ClinVar itself flags conflicting interpretations"
    else:
        severity = "medium"
        headline = "Submitters disagree between pathogenic and uncertain"

    # A majority label is not a consensus when the dissenters carry equal or
    # better review status, so report the dissenters' star ratings too.
    majority = counts.most_common(1)[0][0]
    dissent_stars = sorted(
        {a.stars for a in assertions if a.significance != majority}, reverse=True
    )
    star_note = (
        f" Dissenting records carry up to {dissent_stars[0]} review star(s)."
        if dissent_stars
        else ""
    )

    return Conflict(
        kind="submitter_discordance",
        severity=severity,  # type: ignore[arg-type]
        headline=headline,
        detail=(
            f"{len(assertions)} ClinVar assertions: {spread}.{star_note} "
            "An aggregate label hides this; the disagreement is the finding."
        ),
        sources=_clinvar_link(bundle),
    )


# ==========================================================================
# 2. The evidence base is weak
# ==========================================================================


def weak_evidence_base(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    cv = bundle.clinvar
    if cv is None or not cv.assertions:
        return None

    total = len(cv.assertions)
    no_criteria = [a for a in cv.assertions if a.stars == 0]
    if not no_criteria:
        return None

    share = len(no_criteria) / total
    best_stars = max((a.stars for a in cv.assertions), default=0)

    if best_stars >= 3:
        # An expert panel or practice guideline record outweighs the noise.
        return None
    if share < 0.5 and best_stars >= 2:
        return None

    severity = "high" if best_stars == 0 else "medium"
    return Conflict(
        kind="weak_evidence_base",
        severity=severity,  # type: ignore[arg-type]
        headline=(
            "No ClinVar record states assertion criteria"
            if best_stars == 0
            else "Most ClinVar records state no assertion criteria"
        ),
        detail=(
            f"{len(no_criteria)} of {total} assertions are 'no assertion criteria "
            f"provided' (0 review stars); the best available is {best_stars} star(s). "
            "Such records are opinions without a stated method, and interfaces "
            "usually render them with the same weight as curated ones."
        ),
        sources=_clinvar_link(bundle),
    )


# ==========================================================================
# 3. The assertion is old
# ==========================================================================


def stale_assertion(
    bundle: EvidenceBundle, t: Thresholds = DEFAULT, today: date | None = None
) -> Conflict | None:
    cv = bundle.clinvar
    if cv is None:
        return None
    dated = [a for a in cv.actionable_assertions if a.last_evaluated]
    if not dated:
        return None

    today = today or date.today()
    newest = max(a.last_evaluated for a in dated)  # type: ignore[type-var]
    age_years = (today - newest).days / 365.25
    if age_years < t.stale_assertion_years:
        return None

    oldest = min(a.last_evaluated for a in dated)  # type: ignore[type-var]
    return Conflict(
        kind="stale_assertion",
        severity="medium",
        headline=f"No ClinVar assertion has been re-evaluated in {age_years:.0f} years",
        detail=(
            f"Most recent evaluation {newest.isoformat()}, oldest "
            f"{oldest.isoformat()}. Population databases and predictors have changed "
            "substantially since; an assertion from before gnomAD may rest on "
            "frequency data that no longer holds."
        ),
        sources=_clinvar_link(bundle),
    )


# ==========================================================================
# 4. Frequency contradicts a pathogenic claim
# ==========================================================================


def frequency_versus_claim(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    cv = bundle.clinvar
    af = bundle.max_af()
    freq = bundle.best_frequency()
    if cv is None or af is None or freq is None:
        return None
    if freq.an is not None and freq.an < t.min_allele_number:
        return None

    pathogenic = [a for a in cv.actionable_assertions if a.significance.is_pathogenic_side]
    if not pathogenic:
        return None

    if af >= t.ba1_af:
        severity, bound = "high", f"the {t.ba1_af * 100:g}% BA1 stand-alone benign threshold"
    elif af >= t.bs1_af:
        severity, bound = "medium", f"the {t.bs1_af * 100:g}% BS1 working threshold"
    else:
        return None

    where = freq.popmax_population or freq.source
    return Conflict(
        kind="frequency_versus_claim",
        severity=severity,  # type: ignore[arg-type]
        headline="Reported pathogenic, but too common in the population",
        detail=(
            f"{len(pathogenic)} ClinVar assertion(s) call this pathogenic or likely "
            f"pathogenic, yet it reaches {af * 100:.3g}% in {where} "
            f"(AN {freq.an:,} in {freq.source}) -- above {bound}. Either the disorder "
            "is far more common or less penetrant than assumed, or the assertion "
            "predates this frequency data."
        ),
        sources=_clinvar_link(bundle) + _gnomad_link(bundle),
    )


# ==========================================================================
# 5. Predictors split
# ==========================================================================


def predictor_split(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    called = [p for p in bundle.predictors if p.call in ("damaging", "benign")]
    if len(called) < t.predictor_split_min_tools:
        return None

    counts = Counter(p.call for p in called)
    damaging, benign = counts.get("damaging", 0), counts.get("benign", 0)
    minority = min(damaging, benign)
    if minority == 0:
        return None

    fraction = minority / len(called)
    if fraction < t.predictor_split_min_minority_fraction:
        return None

    minority_call = "benign" if benign < damaging else "damaging"
    dissenters = ", ".join(
        sorted(p.name for p in called if p.call == minority_call)
    )

    return Conflict(
        kind="predictor_split",
        severity="medium" if fraction < 0.4 else "high",
        headline=f"Computational predictors split {damaging}-{benign}",
        detail=(
            f"Of {len(called)} predictors with a call, {damaging} say damaging and "
            f"{benign} say benign ({fraction * 100:.0f}% minority). Dissenting "
            f"{minority_call}: {dissenters}. Tools sharing training data are not "
            "independent, so agreement among them is weaker evidence than its count "
            "suggests -- which is why PP3 here rests on one calibrated predictor "
            "rather than a tally."
        ),
    )


# ==========================================================================
# 6. Our worksheet disagrees with ClinVar
# ==========================================================================

_DIRECTION = {
    PATHOGENIC: "pathogenic",
    LIKELY_PATHOGENIC: "pathogenic",
    BENIGN: "benign",
    LIKELY_BENIGN: "benign",
}


def computed_versus_reported(
    bundle: EvidenceBundle,
    classification: Classification | None = None,
    t: Thresholds = DEFAULT,
) -> Conflict | None:
    cv = bundle.clinvar
    if cv is None or classification is None:
        return None

    assertions = cv.actionable_assertions
    if not assertions:
        return None

    reported = Counter(a.significance for a in assertions).most_common(1)[0][0]
    reported_direction = (
        "pathogenic" if reported.is_pathogenic_side
        else "benign" if reported.is_benign_side
        else None
    )
    computed_direction = _DIRECTION.get(classification.call)

    if reported_direction is None or computed_direction is None:
        return None
    if reported_direction == computed_direction:
        return None

    return Conflict(
        kind="computed_versus_reported",
        severity="high",
        headline=f"This worksheet reaches {classification.call}; ClinVar's majority is "
        f"{reported.value.replace('_', ' ')}",
        detail=(
            f"{classification.reasoning} The difference is expected in one direction: "
            "this worksheet sees only public data and cannot apply the family, "
            "functional and case-level criteria a submitting lab had. Where it lands "
            "on the *opposite* side, one of the two is using evidence the other lacks, "
            "and that is worth finding out which."
        ),
        sources=_clinvar_link(bundle),
    )


# ==========================================================================
# 7. Germline criteria applied to a somatic variant
# ==========================================================================


def somatic_germline_misuse(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    cv = bundle.clinvar
    if cv is None:
        return None

    origins = cv.origins
    somatic_only = bool(origins) and origins <= {"somatic", "unknown", "not applicable", "not provided"}
    has_somatic = "somatic" in origins
    if not has_somatic:
        return None

    germline = {"germline", "maternal", "paternal", "inherited", "de novo", "biparental"}
    has_germline = bool(origins & germline)

    if somatic_only:
        severity = "high"
        headline = "Somatic variant: the germline ACMG framework does not apply"
        detail = (
            f"Every ClinVar assertion with a stated origin is somatic or unspecified "
            f"({', '.join(sorted(origins))}). The ACMG/AMP 2015 criteria scored above "
            "are for germline variants. Somatic variants are classified by clinical "
            "actionability under the AMP/ASCO/CAP tiers, which ask a different "
            "question: not whether the variant causes inherited disease, but whether "
            "it is therapeutically or prognostically meaningful in this tumour type. "
            "Reading the worksheet above as a germline result is a category error."
        )
    else:
        severity = "medium"
        headline = "Both somatic and germline assertions exist for this variant"
        detail = (
            f"ClinVar assertion origins include {', '.join(sorted(origins))}. The same "
            "nucleotide change can be a recurrent somatic driver and a rare germline "
            "finding; the two are classified under different frameworks and their "
            "evidence should not be pooled. Check which context applies before using "
            "the worksheet above."
        )
        if not has_germline:
            return None

    return Conflict(
        kind="somatic_germline_misuse",
        severity=severity,  # type: ignore[arg-type]
        headline=headline,
        detail=detail,
        sources=_clinvar_link(bundle),
    )


# ==========================================================================
# 8. The input itself was ambiguous
# ==========================================================================


def ambiguous_input(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    siblings = bundle.alternative_alleles
    if not siblings:
        return None
    if bundle.query.kind not in ("rsid", "gene_protein", "unknown"):
        return None

    shown = ", ".join(siblings[:4])
    more = f" (+{len(siblings) - 4} more)" if len(siblings) > 4 else ""
    return Conflict(
        kind="ambiguous_input",
        severity="medium",
        headline="This position carries more than one alternate allele",
        detail=(
            f"'{bundle.query.raw}' resolved to {bundle.query.ref}>{bundle.query.alt}, "
            f"but the same position also has: {shown}{more}. An rsID names a position, "
            "not an allele, so a multi-allelic site hides several variants behind one "
            "identifier -- and they can have opposite classifications. The allele "
            "reported here is the one carrying the most ClinVar evidence; confirm it "
            "is the one you meant."
        ),
    )


# ==========================================================================
# 9. Annotators disagree about the consequence
# ==========================================================================


def annotator_disagreement(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    difference = consequence_disagreement(bundle.consequence, bundle.consequence_vep)
    if not difference:
        return None

    return Conflict(
        kind="annotator_disagreement",
        severity="high",
        headline=f"snpEff and VEP disagree on the consequence: {difference}",
        detail=(
            f"snpEff: {bundle.consequence.most_severe if bundle.consequence else 'none'} "
            f"on {bundle.consequence.transcript_id if bundle.consequence else '-'}. "
            f"VEP: {bundle.consequence_vep.most_severe if bundle.consequence_vep else 'none'} "
            f"on {bundle.consequence_vep.transcript_id if bundle.consequence_vep else '-'}. "
            "Which is right decides which criteria are even applicable: PVS1 needs a "
            "null variant, PP2 and PP3 need missense, BP7 needs synonymous. The "
            "worksheet above follows snpEff on the transcript ClinVar numbers against."
        ),
    )


# ==========================================================================
# 10. Literature exists that the records do not reflect
# ==========================================================================


def unreviewed_literature(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    n_papers = bundle.literature_total or len(bundle.literature_pmids)
    if n_papers < 10:
        return None
    cv = bundle.clinvar
    if cv is None:
        return None
    best_stars = max((a.stars for a in cv.assertions), default=0)
    if best_stars >= 2:
        return None

    return Conflict(
        kind="unreviewed_literature",
        severity="low",
        headline=f"{n_papers} publications mention this variant, but no record above "
        f"{best_stars} review star(s)",
        detail=(
            "A substantial literature exists while the strongest ClinVar record still "
            "states weak or no assertion criteria. The published evidence may simply "
            "not have been curated yet, which is where a manual review has the most to "
            "add -- particularly for PS3 and BS3."
        ),
        sources=_clinvar_link(bundle),
    )


# ==========================================================================
# 11. The worksheet contradicts itself
# ==========================================================================

_STRENGTH_WEIGHT = {
    "stand_alone": 8,
    "very_strong": 8,
    "strong": 4,
    "moderate": 2,
    "supporting": 1,
}


def internal_contradiction(
    bundle: EvidenceBundle,
    criteria: list[CriterionResult] | None = None,
    t: Thresholds = DEFAULT,
) -> Conflict | None:
    """Both a pathogenic and a benign criterion are met on the same variant.

    This is not an error in the framework. It is the single most informative
    thing a worksheet can show, because it localises the disagreement to
    specific lines of evidence rather than leaving it in an aggregate label.
    The textbook case is a common polymorphism that missense predictors call
    damaging: BA1 and PP3 are met at once, and only one of them can be guiding
    a correct conclusion.
    """
    if not criteria:
        return None

    met_path = [c for c in criteria if c.status == "met" and c.direction == "pathogenic"]
    met_benign = [c for c in criteria if c.status == "met" and c.direction == "benign"]
    if not met_path or not met_benign:
        return None

    def weight(items: list[CriterionResult]) -> int:
        return sum(_STRENGTH_WEIGHT.get(c.strength.value if c.strength else "", 0) for c in items)

    path_weight, benign_weight = weight(met_path), weight(met_benign)
    path_codes = ", ".join(f"{c.code}" for c in met_path)
    benign_codes = ", ".join(f"{c.code}" for c in met_benign)

    # Both directions carrying real weight is a genuine standoff; one side
    # clearly dominating is worth noting but not alarming.
    severity = "high" if min(path_weight, benign_weight) >= 4 else "medium"

    stronger = "benign" if benign_weight > path_weight else "pathogenic"
    margin = (
        "the benign evidence is stronger"
        if benign_weight > path_weight
        else "the pathogenic evidence is stronger"
        if path_weight > benign_weight
        else "neither side is stronger"
    )

    return Conflict(
        kind="internal_contradiction",
        severity=severity,  # type: ignore[arg-type]
        headline=f"Pathogenic and benign criteria are both met ({path_codes} vs {benign_codes})",
        detail=(
            f"Met toward pathogenic: {path_codes}. Met toward benign: {benign_codes}. "
            f"By ACMG strength weighting {margin} ({stronger} side). "
            "Where a frequency criterion opposes a computational one, the frequency "
            "observation is a direct measurement in people and the predictor output is "
            "a model inference, so they do not deserve equal credence despite both "
            "appearing as met rows above."
        ),
    )


# ==========================================================================
# Runner
# ==========================================================================


def detect_all(
    bundle: EvidenceBundle,
    criteria: list[CriterionResult] | None = None,
    classification: Classification | None = None,
    t: Thresholds = DEFAULT,
    today: date | None = None,
) -> list[Conflict]:
    """Run every detector and return the conflicts, most severe first."""
    found: list[Conflict | None] = [
        somatic_germline_misuse(bundle, t),
        submitter_discordance(bundle, t),
        frequency_versus_claim(bundle, t),
        computed_versus_reported(bundle, classification, t),
        annotator_disagreement(bundle, t),
        weak_evidence_base(bundle, t),
        stale_assertion(bundle, t, today=today),
        internal_contradiction(bundle, criteria, t),
        predictor_split(bundle, t),
        ambiguous_input(bundle, t),
        unreviewed_literature(bundle, t),
    ]
    conflicts = [c for c in found if c is not None]
    conflicts.sort(key=lambda c: _SEVERITY_ORDER.get(c.severity, 9))
    return conflicts
