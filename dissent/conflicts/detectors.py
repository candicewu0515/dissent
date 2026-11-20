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
from dissent.baseline.rates import load_baseline
from dissent.models import (
    UNKNOWN_REVIEW_STARS,
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


def record_discordance(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    """Disagreement between ClinVar records for one variant.

    The unit here is the RCV, which is a variant-plus-condition aggregate, not
    a submission. That distinction matters and an earlier version of this
    detector got it wrong, reporting "16 submitters disagree" for BRAF V600E
    when the records were:

        RCV000014992  Pathogenic           somatic   Carcinoma of colon
        RCV000067669  Pathogenic           somatic   Melanoma
        RCV004018627  Likely pathogenic    germline  Cardiovascular phenotype
        RCV005089260  Uncertain            germline  RASopathy

    Those are four different questions, not four answers to one. Genuine
    submitter-level discordance lives in SCV records, which are not carried in
    this payload. So the detector now says what it can see: records differ, and
    whether they differ *across conditions* or *within one*. Disagreement
    inside a single condition is the serious kind.
    """
    cv = bundle.clinvar
    if cv is None:
        return None
    records = cv.actionable_assertions
    if len(records) < 2:
        return None

    counts = Counter(a.significance for a in records)
    has_path = any(s.is_pathogenic_side for s in counts)
    has_benign = any(s.is_benign_side for s in counts)
    has_vus = Significance.UNCERTAIN in counts
    flagged = Significance.CONFLICTING in counts

    if not (flagged or (has_path and has_benign) or (has_path and has_vus)):
        return None

    # Group by condition. Disagreement within one condition cannot be explained
    # away as two different clinical questions.
    by_condition: dict[str, set[Significance]] = {}
    for a in records:
        key = (a.condition or "unspecified condition").strip().lower()
        by_condition.setdefault(key, set()).add(a.significance)
    within = {c: sigs for c, sigs in by_condition.items() if len(sigs) > 1}

    spread = ", ".join(
        f"{n}x {sig.value.replace('_', ' ')}" for sig, n in counts.most_common()
    )
    n_submitters = sum(a.number_submitters for a in records) or None
    submitter_note = (
        f" across {n_submitters} submissions" if n_submitters and n_submitters > len(records) else ""
    )

    if within:
        severity = "high"
        condition = next(iter(within))
        headline = "ClinVar records disagree about the same condition"
        context = (
            f"For '{condition}' alone the records carry "
            f"{', '.join(sorted(s.value.replace('_', ' ') for s in within[condition]))}. "
            "That is a real disagreement, not two different clinical questions."
        )
    elif has_path and has_benign:
        severity = "high"
        headline = "ClinVar records span the pathogenic/benign divide"
        context = (
            f"The records cover {len(by_condition)} distinct conditions, so part of "
            "the spread may reflect different clinical questions rather than "
            "disagreement. Check which condition applies to your case."
        )
    elif flagged:
        severity = "high"
        headline = "ClinVar itself flags conflicting classifications"
        context = "ClinVar's own aggregation could not reconcile the submitted records."
    else:
        severity = "medium"
        headline = "ClinVar records range from pathogenic to uncertain"
        context = (
            f"Spread across {len(by_condition)} distinct condition(s)."
            if len(by_condition) > 1
            else "All for the same condition."
        )

    majority = counts.most_common(1)[0][0]
    dissent_stars = sorted(
        {a.stars for a in records if a.significance != majority}, reverse=True
    )
    star_note = (
        f" Dissenting records carry up to {dissent_stars[0]} of 4 review stars."
        if dissent_stars
        else ""
    )

    return Conflict(
        kind="record_discordance",
        severity=severity,  # type: ignore[arg-type]
        headline=headline,
        detail=(
            f"{len(records)} ClinVar records{submitter_note}: {spread}.{star_note} "
            f"{context} An aggregate label hides this; the disagreement is the finding."
        ),
        sources=_clinvar_link(bundle),
    )


# Kept so existing imports and tests continue to resolve.
submitter_discordance = record_discordance


# ==========================================================================
# 2. The evidence base is weak
# ==========================================================================


def weak_evidence_base(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> Conflict | None:
    cv = bundle.clinvar
    if cv is None or not cv.assertions:
        return None

    known = [a for a in cv.assertions if a.stars != UNKNOWN_REVIEW_STARS]
    if not known:
        return None

    total = len(known)
    no_criteria = [a for a in known if a.stars == 0]
    if not no_criteria:
        return None

    share = len(no_criteria) / total
    best_stars = max(a.stars for a in known)

    if best_stars >= 3:
        # An expert panel or practice guideline record outweighs the noise.
        return None
    if share < 0.5 and best_stars >= 2:
        return None

    # The headline must match the data: "most" is only true above half.
    if best_stars == 0:
        severity, headline = "high", "No ClinVar record states assertion criteria"
    elif share > 0.5:
        severity, headline = "medium", "Most ClinVar records state no assertion criteria"
    else:
        severity, headline = "medium", "Some ClinVar records state no assertion criteria"

    return Conflict(
        kind="weak_evidence_base",
        severity=severity,  # type: ignore[arg-type]
        headline=headline,
        detail=(
            f"{len(no_criteria)} of {total} records are 'no assertion criteria "
            f"provided' (0 review stars); the best available is {best_stars} of 4. "
            "Such records are opinions without a stated method, and interfaces "
            "usually render them with the same weight as curated ones."
        ),
        sources=_clinvar_link(bundle),
    )


def unrecognised_review_status(
    bundle: EvidenceBundle, t: Thresholds = DEFAULT
) -> Conflict | None:
    """ClinVar used a review status this build does not know.

    Surfaced rather than scored as zero. ClinVar has reworded these labels
    before; silently treating an unknown one as "no assertion criteria" would
    turn every well-curated record in the database into a high-severity
    evidence-quality warning, and the warning would look entirely plausible.
    """
    cv = bundle.clinvar
    if cv is None:
        return None
    unknown = {
        a.review_status for a in cv.assertions if a.stars == UNKNOWN_REVIEW_STARS
    }
    unknown.discard(None)
    if not unknown:
        return None

    return Conflict(
        kind="unrecognised_review_status",
        severity="low",
        headline="ClinVar used a review status this build does not recognise",
        detail=(
            f"Unmapped: {', '.join(sorted(str(u) for u in unknown))}. Those records "
            "are excluded from the evidence-quality assessment rather than scored "
            "as unreviewed. ClinVar renames these labels periodically; the star "
            "ratings shown may understate the curation behind this variant."
        ),
        sources=_clinvar_link(bundle),
    )


# ==========================================================================
# 3. The assertion is old
# ==========================================================================


def stale_assertion(
    bundle: EvidenceBundle, t: Thresholds = DEFAULT, today: date | None = None
) -> Conflict | None:
    """Old classifications sitting underneath a recent one.

    An earlier version compared only the *newest* evaluation date, so a single
    fresh submission marked the whole variant current. That is the opposite of
    useful: the risk is precisely an assertion from before gnomAD existed,
    still contributing to the aggregate label, with a recent record above it
    drawing the eye. BRAF V600E carries a record last evaluated in 2009 and was
    reported as fresh.
    """
    cv = bundle.clinvar
    if cv is None:
        return None
    dated = [a for a in cv.actionable_assertions if a.last_evaluated]
    if not dated:
        return None

    today = today or date.today()

    def years(d: date) -> float:
        return (today - d).days / 365.25

    stale = [a for a in dated if years(a.last_evaluated) >= t.stale_assertion_years]
    if not stale:
        return None

    newest = max(a.last_evaluated for a in dated)
    oldest = min(a.last_evaluated for a in stale)
    undated = len(cv.actionable_assertions) - len(dated)

    everything_old = len(stale) == len(dated)
    severity = "medium" if everything_old else "low"

    if everything_old:
        headline = (
            f"No ClinVar record has been re-evaluated in "
            f"{years(newest):.1f} years"
        )
    else:
        headline = (
            f"{len(stale)} of {len(dated)} ClinVar records have not been "
            f"re-evaluated in {t.stale_assertion_years}+ years"
        )

    undated_note = (
        f" A further {undated} record(s) carry no evaluation date at all."
        if undated
        else ""
    )

    # The authoritative record being the oldest is worth saying out loud: the
    # reader is otherwise drawn to recent submissions that carry less weight.
    best_stars = max((a.stars for a in dated), default=0)
    if best_stars >= 2:
        best_dates = [a.last_evaluated for a in dated if a.stars == best_stars]
        newest_authoritative = max(best_dates)
        if newest_authoritative < newest:
            undated_note += (
                f" Note the best-reviewed record ({best_stars} of 4 stars) is itself "
                f"from {newest_authoritative.isoformat()}, older than the most recent "
                "submission, which carries less review weight."
            )

    return Conflict(
        kind="stale_assertion",
        severity=severity,  # type: ignore[arg-type]
        headline=headline,
        detail=(
            f"Oldest stale evaluation {oldest.isoformat()} "
            f"({years(oldest):.1f} years ago); most recent across all records "
            f"{newest.isoformat()}.{undated_note} Population databases and "
            "predictors have changed substantially since: an assertion made "
            "before gnomAD may rest on frequency data that no longer holds, and "
            "it still contributes to the aggregate label."
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
    t_criteria: list[CriterionResult] | None = None,
) -> Conflict | None:
    """Compare this worksheet's call with ClinVar's.

    Two cases matter, and an earlier version reported only the first.

    *Opposite directions* -- one side says pathogenic, the other benign -- is
    rare and alarming.

    *This worksheet says uncertain while a well-reviewed ClinVar record is
    confident* is far more common and was previously silent, because the call
    had no direction to compare. That silence hid the most useful case there
    is: BRCA1 c.181T>G is classified Pathogenic by an expert panel at three
    review stars, this engine reaches Uncertain significance from public data
    alone, and the conflict panel said "no conflicts detected". The gap is not
    an error in either party -- it is the measure of how much of the evidence
    lives outside public databases, which is exactly what a reader needs to
    know before trusting a worksheet.
    """
    cv = bundle.clinvar
    if cv is None or classification is None:
        return None

    records = cv.actionable_assertions
    if not records:
        return None

    reported = Counter(a.significance for a in records).most_common(1)[0][0]
    reported_direction = (
        "pathogenic" if reported.is_pathogenic_side
        else "benign" if reported.is_benign_side
        else None
    )
    computed_direction = _DIRECTION.get(classification.call)
    best_stars = max((a.stars for a in cv.assertions), default=0)
    reported_text = reported.value.replace("_", " ")

    if reported_direction is None:
        return None

    # Case 1: the two land on opposite sides.
    if computed_direction and computed_direction != reported_direction:
        return Conflict(
            kind="computed_versus_reported",
            severity="high",
            headline=(
                f"This worksheet reaches {classification.call}; ClinVar's majority "
                f"is {reported_text}"
            ),
            detail=(
                f"{classification.reasoning} The difference is expected in one "
                "direction: this worksheet sees only public data and cannot apply "
                "the family, functional and case-level criteria a submitting lab "
                "had. Where it lands on the *opposite* side, one of the two is "
                "using evidence the other lacks, and that is worth finding out which."
            ),
            sources=_clinvar_link(bundle),
        )

    # Case 2: we are uncertain and ClinVar is not.
    if computed_direction is None:
        manual = [c.code for c in (t_criteria or []) if c.status == "manual"]
        confident = best_stars >= 2
        return Conflict(
            kind="computed_versus_reported",
            severity="medium" if confident else "low",
            headline=(
                f"ClinVar reports {reported_text}; this worksheet cannot get there "
                "from public data"
            ),
            detail=(
                f"The strongest ClinVar record carries {best_stars} of 4 review "
                f"stars. {classification.reasoning} "
                + (
                    "A review at this level rests on evidence that is not public: "
                    "segregation in families, functional assays, case series and "
                    "phenotype. "
                    if confident
                    else "The supporting records are weakly reviewed, so the gap may "
                    "be in their evidence rather than in this one's. "
                )
                + (
                    f"The criteria that would close the gap are the ones marked for "
                    f"human review: {', '.join(manual[:8])}"
                    + (f" and {len(manual) - 8} more." if len(manual) > 8 else ".")
                    if manual
                    else ""
                )
            ),
            sources=_clinvar_link(bundle),
        )

    return None


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

    if benign_weight > path_weight:
        margin = "the benign evidence is stronger"
    elif path_weight > benign_weight:
        margin = "the pathogenic evidence is stronger"
    else:
        margin = "neither side outweighs the other"

    return Conflict(
        kind="internal_contradiction",
        severity=severity,  # type: ignore[arg-type]
        headline=f"Pathogenic and benign criteria are both met ({path_codes} vs {benign_codes})",
        detail=(
            f"Met toward pathogenic: {path_codes}. Met toward benign: {benign_codes}. "
            f"By ACMG strength weighting, {margin}. "
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
        record_discordance(bundle, t),
        frequency_versus_claim(bundle, t),
        computed_versus_reported(bundle, classification, t, criteria),
        annotator_disagreement(bundle, t),
        weak_evidence_base(bundle, t),
        unrecognised_review_status(bundle, t),
        stale_assertion(bundle, t, today=today),
        internal_contradiction(bundle, criteria, t),
        predictor_split(bundle, t),
        ambiguous_input(bundle, t),
        unreviewed_literature(bundle, t),
    ]
    conflicts = [c for c in found if c is not None]
    _attach_baselines(conflicts, bundle.gene_symbol)
    conflicts.sort(key=lambda c: _SEVERITY_ORDER.get(c.severity, 9))
    return conflicts


# Which measured rate answers "is this unusual?" for each detector.
_BASELINE_FOR = {
    "record_discordance": "any_label_disagreement_among_multi_submitter",
    "weak_evidence_base": "every_record_lacks_criteria",
    "stale_assertion": "stale_5y",
    "frequency_versus_claim": "pathogenic_vs_benign",
    "internal_contradiction": None,
}

_BASELINE_PHRASING = {
    "any_label_disagreement_among_multi_submitter": (
        "{rate}% of variants with more than one record show this"
    ),
    "every_record_lacks_criteria": (
        "{rate}% of variants have no record stating assertion criteria"
    ),
    "stale_5y": "{rate}% of variants have not been re-evaluated in five years",
    "pathogenic_vs_benign": (
        "only {rate}% of variants carry an outright pathogenic-versus-benign conflict"
    ),
}


def _attach_baselines(conflicts: list[Conflict], gene_symbol: str | None) -> None:
    """Add the measured frequency of each conflict type, where one exists.

    "Submitters disagree about this variant" is an anecdote until the reader
    knows whether that is routine. These rates come from ClinVar's own
    submission release, so they are a measurement rather than an impression.
    """
    baseline = load_baseline()
    if baseline is None:
        return

    gene = baseline.gene(gene_symbol)
    gene_rank = baseline.gene_rank(gene_symbol)

    for conflict in conflicts:
        key = _BASELINE_FOR.get(conflict.kind)
        if not key:
            continue
        rate = baseline.rates_pct.get(key)
        if rate is None:
            continue

        phrasing = _BASELINE_PHRASING.get(key)
        if not phrasing:
            continue
        note = (
            f"Measured across all {baseline.variants:,} ClinVar variants: "
            + phrasing.format(rate=f"{rate:g}")
            + "."
        )

        # Gene context, but only for the detector it actually describes.
        if conflict.kind == "record_discordance" and gene and gene_rank:
            note += (
                f" In {gene.gene} specifically the rate is {gene.rate_pct:.1f}% "
                f"({gene.discordant:,} of {gene.variants:,} variants), ranking it "
                f"{gene_rank[0]:,} of {gene_rank[1]:,} genes."
            )
        conflict.baseline = note
