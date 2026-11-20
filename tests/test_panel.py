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
    LIKELY_BENIGN,
    LIKELY_PATHOGENIC,
    PATHOGENIC,
    UNCERTAIN,
    classify,
)
from dissent.acmg.criteria import evaluate_all
from dissent.acmg.thresholds import DEFAULT
from dissent.models import (
    Classification,
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
        query=VariantQuery(raw="TEST p.Ala1Val", kind="gene_protein", gene="TEST"),
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


def test_coordinates_alone_count_as_resolved():
    """A variant can legitimately be absent from ClinVar and still be real."""
    bundle = EvidenceBundle(
        query=VariantQuery(raw="chr7:140753336A>T", kind="gcoord", chrom="7", pos=140753336)
    )
    report = Report(
        query=bundle.query,
        classification=Classification(call=UNCERTAIN, reasoning="t"),
        bundle=bundle,
    )
    tier, _ = triage(report, today=TODAY)
    assert tier is not UNRESOLVED


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
    assert any("different question" in r for r in reasons)


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


# --------------------------------------------------------------------------
# Tier 3: thin evidence
# --------------------------------------------------------------------------


def test_no_clinvar_record_is_thin_and_says_why():
    report = report_with(UNCERTAIN, [])
    tier, reasons = triage(report, today=TODAY)
    assert tier is THIN
    # The reason has to carry the protein-numbering warning: a GENE p.Change
    # query that resolves to a real variant with no record is most often the
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
def test_every_recorded_bundle_triages_with_a_reason(name):
    """Whatever tier a real variant lands in, it must say why. A row with an
    empty justification is a ranking the reader cannot audit."""
    bundle = load_bundle(name)
    report = build_report(bundle, today=TODAY)
    row = row_from_report(name, report, today=TODAY)
    assert row.tier in {t.key for t in TIERS}
    assert row.reasons, f"{name} got tier {row.tier} with no reason"


def test_hfe_c282y_is_contradicted():
    """Asserted pathogenic at ~7% allele frequency; the computed call is
    benign. The two are on opposite sides of the clinical line."""
    report = build_report(load_bundle("hfe_c282y"), today=TODAY)
    row = row_from_report("HFE C282Y", report, today=TODAY)
    assert row.tier == CONTRADICTED.key
    assert row.computed_call == BENIGN


def test_braf_v600e_is_not_concordant():
    """Somatic assertions mixed into a germline record, most without criteria.
    Whatever tier it lands in, 'nothing to see' is not an option."""
    report = build_report(load_bundle("braf_v600e"), today=TODAY)
    row = row_from_report("BRAF V600E", report, today=TODAY)
    assert row.tier != CONCORDANT.key


def test_tier_rank_matches_the_tier_key():
    """The page sorts on tier_rank and groups on tier. They must not drift."""
    for name in available_bundles():
        report = build_report(load_bundle(name), today=TODAY)
        r = row_from_report(name, report, today=TODAY)
        assert TIERS[r.tier_rank - 1].key == r.tier
