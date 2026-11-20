"""Panel triage: the tier rule, the ordering, the parser and the export.

The rule decides which variants a reader looks at first, so a wrong tier is
worse than a wrong number -- it is a variant quietly moved to the bottom of a
worklist. Every branch is walked here with its boundary, and the recorded
bundles are used as the end-to-end check that the rule survives contact with
real records.
"""

from datetime import date, timedelta

import pytest

from dissent.acmg.combine import (
    BENIGN,
    CONFLICTING,
    LIKELY_BENIGN,
    LIKELY_PATHOGENIC,
    PATHOGENIC,
    UNCERTAIN,
)
from dissent.models import (
    Classification,
    Conflict,
    ClinVarAssertion,
    ClinVarEvidence,
    Consequence,
    EvidenceBundle,
    Report,
    Significance,
    VariantQuery,
)
from dissent.panel import (
    CONCORDANT,
    CONTRADICTED,
    CSV_COLUMNS,
    MAX_PANEL_SIZE,
    SPLIT,
    THIN,
    TIERS,
    UNRESOLVED,
    PanelRow,
    parse_panel,
    row_from_report,
    rows_to_csv,
    sort_rows,
    tier_counts,
    triage,
)
from dissent.report import build_report
from tests.conftest import available_bundles, load_bundle

TODAY = date(2026, 11, 20)


# --------------------------------------------------------------------------
# Builders
# --------------------------------------------------------------------------


def assertion(
    significance: Significance,
    stars: int = 2,
    review: str = "criteria provided, multiple submitters, no conflicts",
    evaluated: date | None = date(2025, 1, 1),
    condition: str = "Some condition",
) -> ClinVarAssertion:
    return ClinVarAssertion(
        accession="RCV000000001",
        significance_raw=significance.value,
        significance=significance,
        review_status=review,
        stars=stars,
        last_evaluated=evaluated,
        condition=condition,
    )


def report_with(call: str, assertions: list[ClinVarAssertion] | None, resolved: bool = True) -> Report:
    """A minimal report. Only the fields the tier rule reads are populated, so
    a test that passes cannot be passing because of something else."""
    bundle = EvidenceBundle(
        query=VariantQuery(
            raw="TEST p.Ala1Val", kind="gene_protein", gene="TEST", protein_change="p.Ala1Val"
        ),
        clinvar=ClinVarEvidence(variant_id=1, gene="TEST", assertions=assertions or [])
        if assertions is not None
        else None,
        consequence=Consequence(most_severe="missense_variant") if resolved else None,
    )
    return Report(
        query=bundle.query,
        gene_symbol="TEST",
        display_name="TEST p.Ala1Val",
        classification=Classification(call=call, reasoning="test", met_codes=["PM2"]),
        bundle=bundle,
    )


# --------------------------------------------------------------------------
# Tier 5: unresolved
# --------------------------------------------------------------------------


def test_unresolved_when_nothing_was_found():
    report = report_with(UNCERTAIN, assertions=None, resolved=False)
    tier, reasons = triage(report, today=TODAY)
    assert tier is UNRESOLVED
    assert reasons


def test_an_annotated_coordinate_is_resolved_without_a_clinvar_record():
    """A variant can legitimately be absent from ClinVar and still be real --
    as long as something recognised the locus."""
    bundle = EvidenceBundle(
        query=VariantQuery(raw="chr7:140753336A>T", kind="gcoord", chrom="7", pos=140753336),
        consequence=Consequence(most_severe="missense_variant"),
    )
    report = Report(
        query=bundle.query,
        classification=Classification(call=UNCERTAIN, reasoning="t"),
        bundle=bundle,
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is THIN
    assert any("No ClinVar record" in r for r in reasons)
    # The protein-numbering advice is for protein-change queries. A coordinate
    # carries no protein numbering to be wrong about.
    assert not any("protein number" in r for r in reasons)


def test_a_coordinate_nothing_recognises_is_unresolved_not_thin():
    """Coordinates parse into a chromosome and a position whether or not a
    variant exists there. `chr9:999999999A>T` used to arrive as "Thin
    evidence" and sort above a fully curated variant."""
    bundle = EvidenceBundle(
        query=VariantQuery(raw="chr9:999999999A>T", kind="gcoord", chrom="9", pos=999999999)
    )
    report = Report(
        query=bundle.query,
        classification=Classification(call=UNCERTAIN, reasoning="t"),
        bundle=bundle,
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is UNRESOLVED
    assert any("no annotation source recognised" in r for r in reasons)


def test_garbage_coordinates_sort_below_a_curated_variant():
    rows = sort_rows([
        PanelRow(query="chr9:999999999A>T", tier=UNRESOLVED.key, tier_rank=UNRESOLVED.rank),
        PanelRow(query="CFTR G542X", tier=CONCORDANT.key, tier_rank=CONCORDANT.rank),
    ])
    assert rows[0].query == "CFTR G542X"


# --------------------------------------------------------------------------
# Tier 1: contradicted
# --------------------------------------------------------------------------


def test_records_opposed_pathogenic_versus_benign_is_tier_one():
    report = report_with(
        UNCERTAIN,
        [assertion(Significance.PATHOGENIC), assertion(Significance.BENIGN)],
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is CONTRADICTED
    assert any("pathogenic vs benign" in r for r in reasons)


@pytest.mark.parametrize(
    "computed,reported",
    [
        (BENIGN, Significance.PATHOGENIC),
        (BENIGN, Significance.LIKELY_PATHOGENIC),
        (LIKELY_BENIGN, Significance.PATHOGENIC),
        (PATHOGENIC, Significance.BENIGN),
        (LIKELY_PATHOGENIC, Significance.BENIGN),
        (PATHOGENIC, Significance.LIKELY_BENIGN),
    ],
)
def test_computed_call_crossing_the_line_is_tier_one(computed, reported):
    tier, _ = triage(report_with(computed, [assertion(reported)]), today=TODAY)
    assert tier is CONTRADICTED


@pytest.mark.parametrize("computed", [BENIGN, PATHOGENIC, UNCERTAIN])
def test_uncertain_reported_never_crosses_the_line(computed):
    """Uncertain is not a side. Calling it one would put most of ClinVar in
    tier 1 and make the tier meaningless."""
    tier, _ = triage(report_with(computed, [assertion(Significance.UNCERTAIN)]), today=TODAY)
    assert tier is not CONTRADICTED


def test_low_penetrance_is_not_treated_as_pathogenic():
    """'Pathogenic, low penetrance' is a distinct ClinVar label. Flattening it
    into Pathogenic manufactures a contradiction against a benign computed
    call -- exactly the bug this project reports in other tools."""
    tier, _ = triage(
        report_with(BENIGN, [assertion(Significance.PATHOGENIC_LOW_PENETRANCE)]), today=TODAY
    )
    assert tier is not CONTRADICTED


# --------------------------------------------------------------------------
# Tier 2: submitters split
# --------------------------------------------------------------------------


def test_pathogenic_against_uncertain_is_tier_two():
    report = report_with(
        UNCERTAIN,
        [assertion(Significance.PATHOGENIC), assertion(Significance.UNCERTAIN)],
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is SPLIT
    assert any("disagree" in r for r in reasons)


def test_split_outranks_thin():
    """A record that is both split and 0-star belongs in the higher tier: the
    disagreement is the finding, the thin evidence is why it happened."""
    report = report_with(
        UNCERTAIN,
        [
            assertion(Significance.PATHOGENIC, stars=0, review="no assertion criteria provided"),
            assertion(Significance.UNCERTAIN, stars=0, review="no assertion criteria provided"),
        ],
    )
    tier, _ = triage(report, today=TODAY)
    assert tier is SPLIT


@pytest.mark.parametrize(
    "off_axis",
    [
        Significance.DRUG_RESPONSE,
        Significance.RISK_FACTOR,
        Significance.PATHOGENIC_LOW_PENETRANCE,
    ],
)
def test_labels_answering_a_different_question_do_not_create_a_split(off_axis):
    """Drug response, risk factor and low penetrance are not claims on the
    pathogenic/benign axis, so they cannot disagree with one. Counting them
    made MTHFR A222V read as '8 records disagree' when the records agreed
    about pathogenicity and some were about warfarin dosing -- which is the
    error this tool exists to find, committed by this tool."""
    report = report_with(
        BENIGN,
        [assertion(Significance.BENIGN), assertion(off_axis)],
    )
    tier, _ = triage(report, today=TODAY)
    assert tier is not SPLIT


def test_a_real_split_still_counts_only_the_axis_records():
    report = report_with(
        UNCERTAIN,
        [
            assertion(Significance.PATHOGENIC),
            assertion(Significance.UNCERTAIN),
            assertion(Significance.DRUG_RESPONSE),
        ],
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is SPLIT
    assert any("2 records disagree" in r for r in reasons)
    assert any("not a pathogenic/benign claim" in r for r in reasons)


def test_the_reported_label_prefers_a_pathogenicity_claim():
    """A 'drug response' headline beside a computed benign call reads as a
    contradiction and is not one, so the label comes from the best-reviewed
    record that is actually on the axis."""
    report = report_with(
        BENIGN,
        [
            assertion(Significance.DRUG_RESPONSE, stars=3, review="reviewed by expert panel"),
            assertion(Significance.BENIGN, stars=2),
        ],
    )
    row = row_from_report("TEST", report, today=TODAY)
    assert row.reported_call == Significance.BENIGN.value


def test_the_reported_label_falls_back_when_nothing_is_on_the_axis():
    report = report_with(UNCERTAIN, [assertion(Significance.DRUG_RESPONSE, stars=2)])
    row = row_from_report("TEST", report, today=TODAY)
    assert row.reported_call == Significance.DRUG_RESPONSE.value


def test_unresolved_rows_report_no_computed_call():
    """The engine still returns 'Uncertain significance' over an empty bundle.
    Printing that next to a failed lookup reads as a result rather than as the
    absence of one."""
    row = row_from_report("HBB E6V", report_with(UNCERTAIN, None, resolved=False), today=TODAY)
    assert row.computed_call is None
    assert not row.resolved


def test_placeholder_labels_do_not_create_a_split():
    """'not provided' is the absence of a classification, not a third opinion."""
    report = report_with(
        PATHOGENIC,
        [assertion(Significance.PATHOGENIC), assertion(Significance.NOT_PROVIDED)],
    )
    tier, _ = triage(report, today=TODAY)
    assert tier is not SPLIT


def test_disagreement_inside_one_condition_is_named_as_such():
    report = report_with(
        UNCERTAIN,
        [
            assertion(Significance.PATHOGENIC, condition="Long QT syndrome"),
            assertion(Significance.UNCERTAIN, condition="Long QT syndrome"),
        ],
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is SPLIT
    assert any("same question" in r for r in reasons)
    assert not any("distinct conditions" in r for r in reasons)


def test_disagreement_across_conditions_is_qualified_rather_than_asserted():
    """Two records about two diseases are not necessarily two opinions."""
    report = report_with(
        UNCERTAIN,
        [
            assertion(Significance.PATHOGENIC, condition="Melanoma"),
            assertion(Significance.UNCERTAIN, condition="Cardiofaciocutaneous syndrome"),
        ],
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is SPLIT
    assert any("distinct conditions" in r for r in reasons)
    assert not any("same question" in r for r in reasons)


def test_a_within_condition_split_outranks_the_across_condition_note():
    """If any one condition has both sides, that is real disagreement and the
    row should say so rather than explain it away."""
    report = report_with(
        UNCERTAIN,
        [
            assertion(Significance.PATHOGENIC, condition="Long QT syndrome"),
            assertion(Significance.UNCERTAIN, condition="Long QT syndrome"),
            assertion(Significance.PATHOGENIC, condition="Something else"),
        ],
    )
    _, reasons = triage(report, today=TODAY)
    assert any("same question" in r for r in reasons)


def test_records_with_no_stated_condition_are_grouped_together():
    """Two records that both omit the condition are about the same unstated
    thing as far as anything here can tell, so they must not read as two
    distinct conditions."""
    report = report_with(
        UNCERTAIN,
        [
            assertion(Significance.PATHOGENIC, condition=""),
            assertion(Significance.UNCERTAIN, condition=""),
        ],
    )
    _, reasons = triage(report, today=TODAY)
    assert any("same question" in r for r in reasons)


# --------------------------------------------------------------------------
# Tier 3: thin evidence
# --------------------------------------------------------------------------


def test_no_clinvar_record_is_thin_and_says_why():
    report = report_with(UNCERTAIN, [])
    tier, reasons = triage(report, today=TODAY)
    assert tier is THIN
    # For a GENE p.Change query the reason has to carry the protein-numbering
    # warning: resolving to a real variant with no record is most often the
    # legacy-numbering trap, and silence there reads as "nothing known".
    assert any("protein number" in r for r in reasons)


def test_zero_stars_is_thin():
    report = report_with(
        PATHOGENIC,
        [assertion(Significance.PATHOGENIC, stars=0, review="no assertion criteria provided")],
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is THIN
    assert any("0 of 4" in r for r in reasons)


def test_one_star_is_thin_regardless_of_submitter_count():
    """ClinVar's one star means a single submitter or conflicting
    classifications. Both are thin, and `number_submitters` comes off an RCV
    aggregate that counts submissions per condition, which is a different
    question -- so the star rating alone decides."""
    a = assertion(Significance.PATHOGENIC, stars=1, review="criteria provided, single submitter")
    a.number_submitters = 7
    tier, _ = triage(report_with(PATHOGENIC, [a]), today=TODAY)
    assert tier is THIN


def test_two_stars_is_not_thin():
    tier, _ = triage(report_with(PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=2)]), today=TODAY)
    assert tier is CONCORDANT


@pytest.mark.parametrize(
    "years,expect_thin",
    [(4.0, False), (4.9, False), (5.01, True), (9.0, True)],
)
def test_staleness_boundary(years, expect_thin):
    """Walk the five-year threshold rather than asserting one side of it."""
    evaluated = TODAY - timedelta(days=int(years * 365.25))
    report = report_with(
        PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=2, evaluated=evaluated)]
    )
    tier, reasons = triage(report, today=TODAY)
    assert (tier is THIN) is expect_thin
    assert any("years ago" in r for r in reasons) is expect_thin


def test_expert_panel_record_outweighs_a_zero_star_one():
    report = report_with(
        PATHOGENIC,
        [
            assertion(Significance.PATHOGENIC, stars=3, review="reviewed by expert panel"),
            assertion(Significance.PATHOGENIC, stars=0, review="no assertion criteria provided"),
        ],
    )
    tier, _ = triage(report, today=TODAY)
    assert tier is CONCORDANT


def test_unknown_review_status_is_reported_as_unknown_not_as_zero_stars():
    """An unrecognised status is unknown, not worthless. Scoring it zero would
    turn every well-curated record into 'no criteria provided' the next time
    ClinVar rewords a label; calling it four stars would vouch for a record
    nobody has checked. It goes in the tier for evidence that cannot be
    assessed, and the reason has to say which it is."""
    report = report_with(
        PATHOGENIC,
        [assertion(Significance.PATHOGENIC, stars=-1, review="reviewed by a committee of cats")],
    )
    tier, reasons = triage(report, today=TODAY)
    assert tier is THIN
    assert not any("0 of 4" in r for r in reasons)
    assert any("not recognised" in r for r in reasons)


def test_a_known_status_still_reports_its_star_count():
    _, reasons = triage(
        report_with(PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=2)]), today=TODAY
    )
    assert any("2 of 4 stars" in r for r in reasons)


def test_one_unrecognised_record_does_not_hide_a_good_one():
    """The star rating is taken from the records that have one, so a single
    reworded label must not drag a well-reviewed variant into tier 3."""
    report = report_with(
        PATHOGENIC,
        [
            assertion(Significance.PATHOGENIC, stars=-1, review="something new"),
            assertion(Significance.PATHOGENIC, stars=3, review="reviewed by expert panel"),
        ],
    )
    tier, _ = triage(report, today=TODAY)
    assert tier is CONCORDANT


# --------------------------------------------------------------------------
# "No disagreement found" has to mean it
# --------------------------------------------------------------------------


def conflict(kind: str, severity: str = "high") -> Conflict:
    return Conflict(kind=kind, severity=severity, headline=f"{kind} headline", detail="d")


def test_a_high_severity_finding_the_rule_cannot_see_lifts_the_row():
    """The tier facts read ClinVar's labels and nothing else, so a variant
    where snpEff and Ensembl VEP disagree about the molecular consequence --
    meaning the worksheet may be classifying the wrong change -- could sit in
    the bottom tier with nothing said about it."""
    report = report_with(PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=3)])
    assert triage(report, today=TODAY)[0] is CONCORDANT

    report.conflicts = [conflict("annotator_disagreement")]
    tier, reasons = triage(report, today=TODAY)
    assert tier is THIN
    assert any("molecular consequence" in r for r in reasons)


def test_a_finding_the_rule_already_accounts_for_does_not_lift_the_row():
    """Otherwise every variant with more than one record lands in tier 3,
    because the detectors restate what the tier already measured."""
    report = report_with(PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=3)])
    report.conflicts = [conflict("record_discordance"), conflict("computed_versus_reported")]
    assert triage(report, today=TODAY)[0] is CONCORDANT


def test_a_medium_severity_finding_does_not_lift_the_row():
    """Medium findings are worth reading and not worth reordering a worklist
    for; otherwise the ranking tracks the detectors rather than the evidence."""
    report = report_with(PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=3)])
    report.conflicts = [conflict("predictor_split", severity="medium")]
    assert triage(report, today=TODAY)[0] is CONCORDANT


def test_the_rule_reads_the_kind_not_the_headline():
    """A rule keyed on prose would change meaning when someone rewords a
    message. The headline here says the opposite of the kind."""
    report = report_with(PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=3)])
    c = conflict("internal_contradiction")
    c.headline = "everything is fine"
    report.conflicts = [c]
    tier, reasons = triage(report, today=TODAY)
    assert tier is THIN
    assert any("both met" in r for r in reasons)
    assert not any("everything is fine" in r for r in reasons)


def test_an_unknown_conflict_kind_is_still_named_readably():
    """A detector added later must not print a bare slug into a worklist."""
    report = report_with(PATHOGENIC, [assertion(Significance.PATHOGENIC, stars=3)])
    report.conflicts = [conflict("some_new_detector")]
    _, reasons = triage(report, today=TODAY)
    assert any("some new detector" in r for r in reasons)


def test_every_labelled_conflict_kind_is_one_the_detectors_can_raise():
    """A label for a kind nothing emits is dead text that looks like coverage."""
    import re
    from pathlib import Path

    from dissent.conflicts import detectors as module
    from dissent.panel import CONFLICT_LABELS, TIER_ACCOUNTS_FOR

    emitted = set(re.findall(r'kind="([a-z_]+)"', Path(module.__file__).read_text()))
    assert set(CONFLICT_LABELS) <= emitted
    assert TIER_ACCOUNTS_FOR <= emitted


def test_every_detector_kind_is_either_accounted_for_or_labelled():
    """A kind in neither set would reach a worklist as a slug, and nobody
    would have decided whether it should affect the order."""
    import re
    from pathlib import Path

    from dissent.conflicts import detectors as module
    from dissent.panel import CONFLICT_LABELS, TIER_ACCOUNTS_FOR

    emitted = set(re.findall(r'kind="([a-z_]+)"', Path(module.__file__).read_text()))
    uncovered = emitted - set(CONFLICT_LABELS) - TIER_ACCOUNTS_FOR
    assert not uncovered, f"detector kinds nobody classified: {sorted(uncovered)}"


# --------------------------------------------------------------------------
# The flag that is deliberately not a tier
# --------------------------------------------------------------------------


def test_confident_reported_uncertain_computed_is_flagged_not_ranked():
    """This is the expected state for most pathogenic variants: segregation and
    functional evidence is not in a public API. It was a tier in the first
    version of this module and swallowed 21 of 30 rows in a test panel."""
    row = row_from_report(
        "TEST", report_with(UNCERTAIN, [assertion(Significance.PATHOGENIC, stars=2)]), today=TODAY
    )
    assert row.evidence_not_public
    assert row.tier == CONCORDANT.key


def test_flag_is_not_set_when_the_reported_call_is_also_uncertain():
    row = row_from_report(
        "TEST", report_with(UNCERTAIN, [assertion(Significance.UNCERTAIN, stars=2)]), today=TODAY
    )
    assert not row.evidence_not_public


def test_flag_is_not_set_when_there_is_no_record():
    row = row_from_report("TEST", report_with(UNCERTAIN, []), today=TODAY)
    assert not row.evidence_not_public


# --------------------------------------------------------------------------
# Records that are not a pathogenicity claim must not answer for one
# --------------------------------------------------------------------------


def test_an_off_axis_record_does_not_lend_its_stars_to_a_pathogenic_claim():
    """A three-star expert-panel record about drug response was donating its
    stars to a one-star pathogenic assertion, moving the row out of thin
    evidence entirely. The row then read "reported pathogenic, 3 of 4 stars,
    no disagreement found" for a claim one submitter made with no
    corroboration."""
    alone = report_with(
        PATHOGENIC,
        [assertion(Significance.PATHOGENIC, stars=1, review="criteria provided, single submitter")],
    )
    assert triage(alone, today=TODAY)[0] is THIN

    with_drug_response = report_with(
        PATHOGENIC,
        [
            assertion(Significance.PATHOGENIC, stars=1, review="criteria provided, single submitter"),
            assertion(Significance.DRUG_RESPONSE, stars=3, review="reviewed by expert panel"),
        ],
    )
    tier, _ = triage(with_drug_response, today=TODAY)
    assert tier is THIN
    assert row_from_report("t", with_drug_response, today=TODAY).best_stars == 1


def test_a_recent_off_axis_record_does_not_make_an_old_claim_current():
    """"When was this last looked at" is a question about the pathogenicity
    claim. A warfarin-dosing record re-evaluated this year was answering it
    for an assertion last examined in 2005."""
    old = assertion(Significance.PATHOGENIC, stars=2, evaluated=date(2005, 1, 1))
    recent_unrelated = assertion(Significance.DRUG_RESPONSE, stars=2, evaluated=date(2026, 1, 1))

    tier, reasons = triage(report_with(PATHOGENIC, [old, recent_unrelated]), today=TODAY)
    assert tier is THIN
    assert any("years ago" in r for r in reasons)


def test_the_off_axis_note_names_the_records_instead_of_guessing():
    """The wording hard-coded "(drug response, risk factor or low penetrance)"
    and printed it over two MTHFR records that both read "Conflicting
    interpretations of pathogenicity" -- telling the reader they answer a
    different question when they say the opposite."""
    report = report_with(
        UNCERTAIN,
        [
            assertion(Significance.PATHOGENIC),
            assertion(Significance.UNCERTAIN),
            assertion(Significance.CONFLICTING),
        ],
    )
    _, reasons = triage(report, today=TODAY)
    joined = " ".join(reasons)
    assert "conflicting" in joined.lower()
    assert "drug response" not in joined


def test_the_no_criteria_denominator_counts_only_the_rated_records():
    """"1 of 3 records state no criteria" was printed where one of the three
    had an unrecognised status and another was not a pathogenicity claim."""
    report = report_with(
        PATHOGENIC,
        [
            assertion(Significance.PATHOGENIC, stars=0, review="no assertion criteria provided"),
            assertion(Significance.PATHOGENIC, stars=2),
            assertion(Significance.DRUG_RESPONSE, stars=2),
        ],
    )
    _, reasons = triage(report, today=TODAY)
    assert any("1 of 2 records state no criteria" in r for r in reasons)


def test_when_nothing_is_on_the_axis_the_records_rate_themselves():
    """A variant whose only records are drug-response ones still has records.
    Reporting it as having none would be worse than rating what is there."""
    report = report_with(UNCERTAIN, [assertion(Significance.DRUG_RESPONSE, stars=2)])
    row = row_from_report("t", report, today=TODAY)
    assert row.records == 1
    assert row.best_stars == 2
    assert row.reported_call == Significance.DRUG_RESPONSE.value


# --------------------------------------------------------------------------
# The engine's own contradiction
# --------------------------------------------------------------------------


def test_the_engine_meeting_both_sides_reaches_tier_one():
    """combine.py returns "Uncertain significance (conflicting evidence)" when
    pathogenic and benign criteria are both satisfied -- the strongest
    internal signal the tool has. _call_side bucketed it with every ordinary
    VUS, so the row could land in "No disagreement found"."""
    report = report_with(CONFLICTING, [assertion(Significance.PATHOGENIC, stars=2)])
    tier, reasons = triage(report, today=TODAY)
    assert tier is CONTRADICTED
    assert any("both met" in r for r in reasons)


def test_an_ordinary_vus_does_not_reach_tier_one():
    """The branch above must key on the conflicting verdict, not on the word
    'uncertain', or every VUS in ClinVar lands in tier 1."""
    report = report_with(UNCERTAIN, [assertion(Significance.PATHOGENIC, stars=2)])
    assert triage(report, today=TODAY)[0] is not CONTRADICTED


def test_the_conflicting_verdict_is_not_flagged_as_private_evidence():
    row = row_from_report(
        "t", report_with(CONFLICTING, [assertion(Significance.PATHOGENIC, stars=2)]), today=TODAY
    )
    assert not row.evidence_not_public


# --------------------------------------------------------------------------
# The flag must not fire where the records themselves disagree
# --------------------------------------------------------------------------


def test_a_contradicted_row_is_not_flagged_as_private_evidence():
    """The row printed "the submitter's evidence is not public -- expected,
    not a flag" directly underneath "these records disagree pathogenic vs
    benign". Two submitters on opposite sides is not one submitter holding
    private evidence."""
    row = row_from_report(
        "t",
        report_with(UNCERTAIN, [assertion(Significance.PATHOGENIC), assertion(Significance.BENIGN)]),
        today=TODAY,
    )
    assert row.tier == CONTRADICTED.key
    assert not row.evidence_not_public


def test_a_split_row_is_not_flagged_as_private_evidence():
    row = row_from_report(
        "t",
        report_with(
            UNCERTAIN, [assertion(Significance.PATHOGENIC), assertion(Significance.UNCERTAIN)]
        ),
        today=TODAY,
    )
    assert not row.evidence_not_public


# --------------------------------------------------------------------------
# ClinVar's own "submitters disagree" label
# --------------------------------------------------------------------------


def test_a_conflicting_aggregate_is_not_read_as_an_absent_record():
    """ClinVar's current wording fell through to OTHER, was dropped as not
    actionable, and the variant was reported as having no ClinVar record --
    while the record present says, literally, that submitters disagree."""
    from dissent.models import normalize_significance

    sig = normalize_significance("Conflicting classifications of pathogenicity")
    report = report_with(UNCERTAIN, [assertion(sig, stars=1)])
    _, reasons = triage(report, today=TODAY)
    assert not any("No ClinVar record" in r for r in reasons)


def test_a_contradicted_row_always_names_what_was_reported():
    """Where every record carries an unreadable review status, reported_label
    stays None while the sides are known, and the row printed "reported: no
    record" beside a reason saying "reported pathogenic"."""
    report = report_with(
        BENIGN, [assertion(Significance.PATHOGENIC, stars=-1, review="brand new wording")]
    )
    row = row_from_report("t", report, today=TODAY)
    assert row.tier == CONTRADICTED.key
    assert row.records == 1
    assert row.reported_call, "a row with a record must not display 'no record'"


# --------------------------------------------------------------------------
# Ordering
# --------------------------------------------------------------------------


def row(tier_rank: int, high: int = 0, conflicts: int = 0, stars: int | None = 2, name: str = "a") -> PanelRow:
    return PanelRow(
        query=name,
        display_name=name,
        tier=TIERS[tier_rank - 1].key,
        tier_rank=tier_rank,
        high_severity=high,
        conflicts=conflicts,
        best_stars=stars,
    )


def test_tier_dominates_the_order():
    out = sort_rows([row(4, high=9, conflicts=9), row(1)])
    assert out[0].tier_rank == 1


def test_within_a_tier_high_severity_comes_first():
    out = sort_rows([row(3, high=0, conflicts=5), row(3, high=1, conflicts=1)])
    assert out[0].high_severity == 1


def test_within_a_tier_fewer_stars_comes_first():
    """A one-star record with a problem is a worse place to be than a
    three-star one with the same problem: nobody has checked it."""
    out = sort_rows([row(3, stars=3, name="b"), row(3, stars=1, name="a")])
    assert out[0].best_stars == 1


def test_missing_star_rating_sorts_last_not_first():
    out = sort_rows([row(3, stars=None, name="b"), row(3, stars=4, name="a")])
    assert out[0].best_stars == 4


def test_tier_counts_cover_every_tier():
    counts = tier_counts([row(1), row(1), row(3)])
    assert set(counts) == {t.key for t in TIERS}
    assert counts[CONTRADICTED.key] == 2
    assert counts[THIN.key] == 1
    assert counts[CONCORDANT.key] == 0


# --------------------------------------------------------------------------
# Input parsing
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("BRAF V600E\nrs334", ["BRAF V600E", "rs334"]),
        ("BRAF V600E, rs334", ["BRAF V600E", "rs334"]),
        ('"BRAF V600E"\n', ["BRAF V600E"]),
        ("  BRAF V600E  \n\n\n", ["BRAF V600E"]),
        ("# a comment\nrs334", ["rs334"]),
        ("rs334\nRS334", ["rs334"]),
        ("", []),
        (",,,\n\n", []),
    ],
)
def test_parse_panel(raw, expected):
    assert parse_panel(raw) == expected


def test_parse_panel_honours_the_limit():
    assert len(parse_panel("\n".join(f"rs{i}" for i in range(200)))) == MAX_PANEL_SIZE


def test_parse_panel_keeps_the_first_of_a_duplicate_pair():
    assert parse_panel("BRAF V600E\nbraf v600e") == ["BRAF V600E"]


# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------


def test_csv_has_a_header_and_one_row_each():
    body = rows_to_csv([row(1, name="x"), row(3, name="y")])
    lines = body.strip().split("\n")
    assert lines[0] == ",".join(CSV_COLUMNS)
    assert len(lines) == 3


def test_csv_quotes_a_reason_containing_a_comma():
    r = row(1, name="x")
    r.reasons = ["two records disagree, badly"]
    body = rows_to_csv([r])
    assert '"two records disagree, badly"' in body


def test_csv_names_the_tier_in_words():
    assert CONTRADICTED.name in rows_to_csv([row(1)])


# --------------------------------------------------------------------------
# End to end, against the recorded bundles
# --------------------------------------------------------------------------


@pytest.mark.parametrize("name", available_bundles())
def test_every_recorded_bundle_triages_with_a_substantive_reason(name):
    """Whatever tier a real variant lands in, it must say why, in words that
    carry information.

    The earlier version asserted `row.tier in TIERS` and `row.reasons` --
    both true by construction, so the test could not fail without a type
    error. These assertions can: a reason has to name a number or a named
    state, not merely exist.
    """
    report = build_report(load_bundle(name), today=TODAY)
    row = row_from_report(name, report, today=TODAY)
    assert row.reasons, f"{name} got tier {row.tier} with no reason"

    joined = " ".join(row.reasons)
    assert len(joined) > 20, f"{name}: reason too thin to audit: {joined!r}"
    informative = any(ch.isdigit() for ch in joined) or any(
        word in joined.lower()
        for word in ("computed", "both met", "not recognised", "no clinvar record", "did not resolve")
    )
    assert informative, f"{name}: reason states nothing checkable: {joined!r}"


@pytest.mark.parametrize("name", available_bundles())
def test_a_recorded_bundle_never_lands_in_the_bottom_tier_with_high_severity_conflicts(name):
    """"No disagreement found" has to mean it. A row carrying a high-severity
    conflict is a row with something to read."""
    report = build_report(load_bundle(name), today=TODAY)
    row = row_from_report(name, report, today=TODAY)
    if row.tier == CONCORDANT.key:
        assert row.high_severity == 0, (
            f"{name}: tier 4 while carrying {row.high_severity} high-severity "
            f"conflicts: {row.conflict_headlines}"
        )


def test_hfe_c282y_is_contradicted():
    """Asserted pathogenic at ~7% allele frequency; the computed call is
    benign. The two are on opposite sides of the clinical line."""
    report = build_report(load_bundle("hfe_c282y"), today=TODAY)
    row = row_from_report("HFE C282Y", report, today=TODAY)
    assert row.tier == CONTRADICTED.key
    assert row.computed_call == BENIGN


def test_braf_v600e_is_split_and_says_the_records_cover_many_conditions():
    """A ClinVar record is an RCV: one aggregate per variant-and-condition,
    not one submitter. BRAF V600E has 16, almost all somatic tumour types
    beside one germline RASopathy, and the row used to assert flatly that
    "submitters" disagreed. The earlier test only checked it was not tier 4,
    which passed for four of the five tiers."""
    report = build_report(load_bundle("braf_v600e"), today=TODAY)
    row = row_from_report("BRAF V600E", report, today=TODAY)
    assert row.tier == SPLIT.key
    joined = " ".join(row.reasons)
    assert "distinct conditions" in joined
    assert "submitters" not in joined.lower()


def test_tier_ranks_are_unique_and_contiguous_from_one():
    """The page sorts on tier_rank and groups on tier.

    The earlier version asserted `TIERS[rank - 1].key == tier` on rows built
    from those same Tier objects -- a tautology. What actually has to hold is
    a property of the TIERS table itself.
    """
    ranks = [t.rank for t in TIERS]
    assert ranks == sorted(ranks)
    assert ranks == list(range(1, len(TIERS) + 1))
    assert len({t.key for t in TIERS}) == len(TIERS)
    assert len({t.name for t in TIERS}) == len(TIERS)


def test_unresolved_is_the_last_tier():
    """sort_rows puts failures at the bottom by rank alone, so a new tier
    appended after UNRESOLVED would sort above it."""
    assert TIERS[-1] is UNRESOLVED
