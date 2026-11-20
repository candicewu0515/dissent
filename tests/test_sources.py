"""Tests for the normalisation layer: raw API payloads in, typed models out.

This is where the bugs live. The code here handles 21 predictors with
incompatible letter codes, fields that are a scalar on one variant and a list
on the next, transcript selection, and three different notions of "which
allele". None of it was covered by the bundle-level fixtures, which are
already-normalised output.

Payloads under `fixtures/raw/` are verbatim responses recorded from the live
APIs; nothing here touches the network.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pytest

from tests.conftest import load_bundle as conftest_load_bundle

from dissent.collect import _dbsnp_citations, _dedupe_statuses, _resolve_coordinates
from dissent.models import (
    UNKNOWN_REVIEW_STARS,
    Significance,
    SourceStatus,
    VariantQuery,
    normalize_significance,
    review_stars,
)
from dissent.normalize.consequence import from_snpeff, transcript_from_preferred_name
from dissent.sources.base import as_list, dig, first, to_float, to_float_all, to_int
from dissent.sources.gnomad import _frequency_from, _graphql_errors
from dissent.sources.myvariant import (
    MyVariantClient,
    _same_position_siblings,
    extract_clinvar,
    extract_frequencies,
    extract_predictors,
)

RAW = Path(__file__).parent / "fixtures" / "raw"


def raw(name: str):
    return json.loads((RAW / f"{name}.json").read_text())


def raw_names() -> list[str]:
    return sorted(p.stem for p in RAW.glob("myvariant_*.json"))


# ==========================================================================
# ClinVar extraction
# ==========================================================================


@pytest.mark.parametrize("name", raw_names())
def test_clinvar_extraction_survives_every_recorded_payload(name):
    evidence = extract_clinvar(raw(name))
    assert evidence is not None
    assert evidence.assertions
    for a in evidence.assertions:
        assert isinstance(a.significance, Significance)
        assert a.stars >= UNKNOWN_REVIEW_STARS


def test_clinvar_records_carry_condition_and_origin():
    """The discordance detector groups by condition, so it must survive extraction."""
    evidence = extract_clinvar(raw("myvariant_braf"))
    assert any(a.condition for a in evidence.assertions)
    assert any(a.origin for a in evidence.assertions)
    assert "somatic" in evidence.origins


def test_clinvar_dates_parse_from_the_recorded_format():
    evidence = extract_clinvar(raw("myvariant_braf"))
    dated = [a.last_evaluated for a in evidence.assertions if a.last_evaluated]
    assert dated
    assert all(isinstance(d, date) for d in dated)
    assert min(dated).year < 2015, "BRAF carries a long-stale assertion"


def test_preferred_name_is_recovered_for_display():
    evidence = extract_clinvar(raw("myvariant_brca1"))
    assert evidence.preferred_name
    assert "BRCA1" in evidence.preferred_name


# ==========================================================================
# Significance and review status vocabulary
# ==========================================================================


@pytest.mark.parametrize(
    "raw_label,expected",
    [
        ("Pathogenic", Significance.PATHOGENIC),
        ("pathogenic", Significance.PATHOGENIC),
        ("Likely pathogenic", Significance.LIKELY_PATHOGENIC),
        ("Uncertain significance", Significance.UNCERTAIN),
        ("Benign", Significance.BENIGN),
        ("Likely benign", Significance.LIKELY_BENIGN),
        ("Conflicting interpretations of pathogenicity", Significance.CONFLICTING),
        # ClinVar's current wording after the 2023/24 relabelling. The star
        # table was updated for it and this one was not, so it fell through to
        # OTHER and was dropped as not actionable -- a variant whose only
        # record says "submitters disagree" was reported as having no ClinVar
        # record at all.
        ("Conflicting classifications of pathogenicity", Significance.CONFLICTING),
        ("not provided", Significance.NOT_PROVIDED),
    ],
)
def test_significance_vocabulary(raw_label, expected):
    assert normalize_significance(raw_label) is expected


@pytest.mark.parametrize(
    "raw_label,expected",
    [
        ("Benign; other", Significance.BENIGN),
        ("Pathogenic; risk factor", Significance.PATHOGENIC),
        ("Pathogenic; drug response", Significance.PATHOGENIC),
        ("Likely pathogenic; other", Significance.LIKELY_PATHOGENIC),
        ("Benign/Likely benign; other", Significance.BENIGN),
        ("Uncertain significance; risk factor", Significance.UNCERTAIN),
    ],
)
def test_a_compound_label_keeps_its_primary_classification(raw_label, expected):
    """ClinVar joins the primary classification to its additional ones with a
    semicolon, primary first. Splitting on the comma alone left these falling
    through to OTHER, where is_actionable drops them -- so a pathogenic record
    vanished from the evidence, taking with it any opposition it was half of.
    """
    sig = normalize_significance(raw_label)
    assert sig is expected
    assert sig.is_actionable


def test_a_semicolon_label_does_not_break_the_comma_ones():
    """"Pathogenic, low penetrance" is a classification in its own right, not
    a primary plus a qualifier, so the whole string has to be tried first."""
    assert (
        normalize_significance("Pathogenic, low penetrance")
        is Significance.PATHOGENIC_LOW_PENETRANCE
    )
    assert (
        normalize_significance("Likely pathogenic, low penetrance")
        is Significance.PATHOGENIC_LOW_PENETRANCE
    )


def test_a_label_nobody_recognises_is_still_other():
    """The fallbacks must not become a way to guess."""
    assert normalize_significance("banana; mango") is Significance.OTHER
    assert normalize_significance("") is Significance.NOT_PROVIDED
    assert normalize_significance(None) is Significance.NOT_PROVIDED


def test_low_penetrance_is_not_flattened_into_pathogenic():
    """Low penetrance is a submitter's answer to 'why is this common?'.

    Collapsing it to Pathogenic manufactured the frequency-versus-claim
    conflict the tool then reported as a finding.
    """
    sig = normalize_significance("Pathogenic, low penetrance")
    assert sig is Significance.PATHOGENIC_LOW_PENETRANCE
    assert not sig.is_pathogenic_side
    assert sig.is_reduced_penetrance


@pytest.mark.parametrize(
    "status,stars",
    [
        ("practice guideline", 4),
        ("reviewed by expert panel", 3),
        ("criteria provided, multiple submitters, no conflicts", 2),
        ("criteria provided, single submitter", 1),
        ("no assertion criteria provided", 0),
        # Absent and empty are unknown, not zero. Zero is ClinVar stating that
        # a submitter gave no assertion criteria; a missing field says nothing
        # about the submitter, and reading it as zero turns a gap in the
        # payload into a claim about someone's work.
        (None, UNKNOWN_REVIEW_STARS),
        ("", UNKNOWN_REVIEW_STARS),
        ("   ", UNKNOWN_REVIEW_STARS),
        ("a status ClinVar has not used before", UNKNOWN_REVIEW_STARS),
        # The 2023/24 relabelling. Both spellings rate the same.
        ("criteria provided, conflicting classifications", 1),
        ("criteria provided, conflicting interpretations", 1),
    ],
)
def test_review_status_star_ratings(status, stars):
    assert review_stars(status) == stars


def test_an_unrecognised_review_status_is_unknown_not_zero():
    """ClinVar renames these labels; a rename must not fake an evidence problem."""
    assert review_stars("criteria provided, some new wording") == UNKNOWN_REVIEW_STARS
    assert UNKNOWN_REVIEW_STARS < 0


# ==========================================================================
# Predictors — the part with 21 incompatible code systems
# ==========================================================================


@pytest.mark.parametrize("name", raw_names())
def test_predictor_extraction_survives_every_recorded_payload(name):
    for p in extract_predictors(raw(name)):
        assert p.name
        assert p.call in (None, "damaging", "benign", "uncertain")


def test_the_letter_p_is_read_per_predictor_not_globally():
    """`P` means three different things, and a shared map would invert two.

    PolyPhen-2: possibly damaging. AlphaMissense: Pathogenic.
    MutationTaster: polymorphism, i.e. benign.
    """
    doc = {
        "dbnsfp": {
            "polyphen2": {"hvar": {"pred": ["P"], "score": [0.7]}},
            "alphamissense": {"pred": ["P"], "score": [0.9]},
            "mutationtaster": {"pred": ["P"], "score": [0.9]},
        }
    }
    calls = {p.name: p.call for p in extract_predictors(doc)}
    assert calls["PolyPhen-2 HVAR"] == "uncertain"
    assert calls["AlphaMissense"] == "damaging"
    assert calls["MutationTaster"] == "benign"


def test_a_scalar_and_a_single_element_list_give_the_same_result():
    """MyVariant returns a scalar for one transcript and a list for several."""
    scalar = extract_predictors({"dbnsfp": {"metasvm": {"pred": "D", "score": 0.85}}})
    listed = extract_predictors({"dbnsfp": {"metasvm": {"pred": ["D"], "score": [0.85]}}})
    assert [(p.name, p.call, p.score) for p in scalar] == [
        (p.name, p.call, p.score) for p in listed
    ]


def test_a_damaging_call_on_any_transcript_makes_the_predictor_damaging():
    doc = {"dbnsfp": {"sift": {"pred": ["T", "D", "T"], "score": [0.4, 0.0, 0.5]}}}
    sift = next(p for p in extract_predictors(doc) if p.name == "SIFT")
    assert sift.call == "damaging"


def test_score_only_predictors_are_marked_as_threshold_derived():
    """REVEL has no categorical call, so its call is inferred and flagged."""
    doc = {"dbnsfp": {"revel": {"score": [0.9]}}}
    revel = next(p for p in extract_predictors(doc) if p.name == "REVEL")
    assert revel.call == "damaging"
    assert revel.call_from_threshold is True


def test_a_predictor_whose_damaging_side_is_low_scores_uses_the_right_direction():
    """ESM1b is damaging below its cutoff, not above."""
    low = extract_predictors({"dbnsfp": {"esm1b": {"score": [-12.0]}}})
    high = extract_predictors({"dbnsfp": {"esm1b": {"score": [-1.0]}}})
    assert next(p for p in low if p.name == "ESM1b").call == "damaging"
    assert next(p for p in high if p.name == "ESM1b").call == "benign"


def test_absent_predictors_are_omitted_not_invented():
    assert extract_predictors({}) == []
    assert extract_predictors({"dbnsfp": {}}) == []


def test_a_real_payload_yields_a_substantial_panel():
    predictors = extract_predictors(raw("myvariant_brca1"))
    assert len(predictors) >= 15
    assert {"REVEL", "AlphaMissense", "SIFT"} <= {p.name for p in predictors}


# ==========================================================================
# Frequencies
# ==========================================================================


@pytest.mark.parametrize("name", raw_names())
def test_frequency_extraction_survives_every_recorded_payload(name):
    for f in extract_frequencies(raw(name)):
        assert f.source
        if f.af is not None:
            assert 0.0 <= f.af <= 1.0
        for value in f.per_population.values():
            assert 0.0 <= value <= 1.0


def test_the_myvariant_fallback_derives_no_popmax():
    """Its per-population fields carry a frequency with no denominator.

    An unguarded maximum over them let a handful of alleles in `af_oth` drive
    a BA1 decision, which is the false-benign direction.
    """
    frequencies = extract_frequencies(raw("myvariant_hfe"))
    assert frequencies
    for f in frequencies:
        assert f.popmax_af is None
        assert f.popmax_population is None
        # The figures are still carried, for display.
        assert f.per_population
        # And the global frequency, which does have an allele number, remains.
        assert f.an


def test_the_myvariant_fallback_names_the_build_it_actually_is():
    """It is gnomAD v2, not v4; the provenance must not claim otherwise."""
    for f in extract_frequencies(raw("myvariant_hfe")):
        assert "v2" in f.source and "MyVariant" in f.source


def test_popmax_from_the_live_api_excludes_bottlenecked_populations():
    """gnomAD excludes these from its own grpmax, and so must this.

    Founder effects inflate specific pathogenic variants in exactly these
    groups, which is the situation where a naive maximum declares a founder
    variant benign.
    """
    node = {
        "ac": 100, "an": 200000, "af": 0.0005, "homozygote_count": 0,
        "populations": [
            {"id": "nfe", "ac": 200, "an": 100000},   # 0.2%
            {"id": "fin", "ac": 4000, "an": 20000},   # 20%, founder effect
            {"id": "asj", "ac": 3000, "an": 10000},   # 30%, founder effect
            {"id": "remaining", "ac": 900, "an": 3000},
        ],
    }
    freq = _frequency_from(node, "test")
    assert freq is not None
    # Reported, because a reader interpreting a Finnish patient needs it.
    assert freq.per_population["Finnish"] == pytest.approx(0.2)
    assert freq.per_population["Ashkenazi Jewish"] == pytest.approx(0.3)
    # But not driving the criteria.
    assert freq.popmax_population == "Non-Finnish European"
    assert freq.popmax_af == pytest.approx(0.002)


def test_gnomad_api_frequencies_ignore_underpowered_populations():
    """A handful of alleles cannot support a BA1 or BS1 decision."""
    node = {
        "ac": 100, "an": 100000, "af": 0.001, "homozygote_count": 0,
        "populations": [
            {"id": "nfe", "ac": 90, "an": 90000},   # powered
            {"id": "asj", "ac": 10, "an": 50},      # not powered: 20% of 50 alleles
        ],
    }
    freq = _frequency_from(node, "test")
    assert freq is not None
    assert "Ashkenazi Jewish" not in freq.per_population
    assert freq.popmax_population == "Non-Finnish European"


def test_a_graphql_payload_with_no_allele_number_is_discarded():
    assert _frequency_from({"ac": 1, "an": 0}, "test") is None
    assert _frequency_from(None, "test") is None


# ==========================================================================
# GraphQL error detection
# ==========================================================================


def test_graphql_errors_are_detected_in_a_200_response():
    """GraphQL reports schema failures with HTTP 200 and an errors array."""
    payload = {"errors": [{"message": 'Cannot query field "populations"'}], "data": None}
    detected = _graphql_errors(payload)
    assert detected and "populations" in detected


def test_a_healthy_graphql_payload_reports_no_error():
    assert _graphql_errors({"data": {"variant": {}}}) is None
    assert _graphql_errors(raw("gnomad_variant")) is None


# ==========================================================================
# Transcript selection
# ==========================================================================


@pytest.mark.parametrize("name", raw_names())
def test_consequence_selection_survives_every_recorded_payload(name):
    doc = raw(name)
    evidence = extract_clinvar(doc)
    preferred = transcript_from_preferred_name(evidence.preferred_name if evidence else None)
    cons = from_snpeff(doc, preferred)
    if cons is None:
        return
    assert cons.most_severe
    assert cons.transcript_id


def test_selection_follows_the_transcript_clinvar_numbers_against():
    doc = raw("myvariant_brca1")
    evidence = extract_clinvar(doc)
    preferred = transcript_from_preferred_name(evidence.preferred_name)
    cons = from_snpeff(doc, preferred)
    assert cons is not None
    assert cons.hgvs_p and "61" in cons.hgvs_p, cons.hgvs_p


def test_a_structural_annotation_is_never_selected():
    """snpEff emits PDB interaction records with HIGH impact."""
    doc = raw("myvariant_tp53")
    cons = from_snpeff(doc, "NM_000546")
    assert cons is not None
    assert cons.most_severe != "structural_interaction_variant"
    assert not cons.transcript_id.startswith("1")


# ==========================================================================
# Hit selection
# ==========================================================================


def test_siblings_are_only_alleles_at_the_same_position():
    """A phrase search matches the same residue change in unrelated genes."""
    best = {"_id": "chr7:g.140753336A>T"}
    hits = [
        best,
        {"_id": "chr7:g.140753336A>C"},
        {"_id": "chr7:g.140753336A>G"},
        {"_id": "chr10:g.66280555A>T"},
    ]
    siblings = _same_position_siblings(best, hits)
    assert set(siblings) == {"chr7:g.140753336A>C", "chr7:g.140753336A>G"}


def test_a_hit_in_the_wrong_gene_is_refused_rather_than_returned():
    query = VariantQuery(raw="TP53 R273H", kind="gene_protein", gene="TP53",
                         protein_change="p.R273H")
    hits = [{"_id": "chr19:g.1622058C>T", "dbnsfp": {"genename": ["TCF3"]}, "_score": 90}]
    assert MyVariantClient._best_hit(hits, query) is None


def test_the_gene_constrained_query_is_tried_before_the_phrase_search():
    """A bare phrase for p.Arg273His ranks TP53 below unrelated genes."""
    query = VariantQuery(raw="TP53 R273H", kind="gene_protein", gene="TP53",
                         protein_change="p.R273H")
    terms = MyVariantClient._search_terms(query)
    assert terms
    assert "dbnsfp.genename:TP53" in terms[0]
    assert any(t == '"p.Arg273His"' for t in terms), "phrase search kept as a fallback"


# ==========================================================================
# Coordinate resolution
# ==========================================================================


def test_coordinates_come_from_the_vcf_block():
    """The top level of a MyVariant document carries no ref or alt."""
    doc = raw("myvariant_braf")
    resolved = _resolve_coordinates(
        VariantQuery(raw="BRAF V600E", kind="gene_protein", gene="BRAF"),
        doc,
        extract_clinvar(doc),
    )
    assert resolved.chrom == "7"
    assert resolved.ref and resolved.alt
    assert resolved.gnomad_variant_id == f"7-{resolved.pos}-{resolved.ref}-{resolved.alt}"


def test_resolution_does_not_overwrite_what_the_user_supplied():
    doc = raw("myvariant_braf")
    query = VariantQuery(raw="x", kind="gcoord", chrom="7", pos=1, ref="A", alt="T")
    resolved = _resolve_coordinates(query, doc, extract_clinvar(doc))
    assert resolved.pos == 1, "an explicit position must win over the document"


def test_dbsnp_citations_are_read_as_integers():
    pmids = _dbsnp_citations(raw("myvariant_braf"))
    assert pmids
    assert all(isinstance(p, int) for p in pmids)


# ==========================================================================
# Source status bookkeeping
# ==========================================================================


def test_a_populated_attempt_replaces_an_empty_one():
    """A source retried after resolution must show one row, not two."""
    statuses = [
        SourceStatus(source="gnomad", ok=True, detail="no genomic coordinates to query"),
        SourceStatus(source="gnomad", ok=True, elapsed_ms=120),
    ]
    merged = _dedupe_statuses(statuses)
    assert len(merged) == 1
    assert merged[0].detail is None


def test_a_success_replaces_a_failure_for_the_same_source():
    statuses = [
        SourceStatus(source="ensembl_vep", ok=False, detail="timeout"),
        SourceStatus(source="ensembl_vep", ok=True, elapsed_ms=90),
    ]
    assert _dedupe_statuses(statuses)[0].ok is True


# ==========================================================================
# Ragged-shape helpers
# ==========================================================================


def test_shape_helpers():
    assert as_list(None) == []
    assert as_list(3) == [3]
    assert as_list([1, None, 2]) == [1, 2]
    assert first([5, 6]) == 5
    assert first(None, "fallback") == "fallback"
    assert dig({"a": {"b": {"c": 1}}}, "a", "b", "c") == 1
    assert dig({"a": [{"b": 2}]}, "a", "b") == 2
    assert dig({"a": 1}, "a", "b") is None
    assert to_float(["0.83"]) == 0.83
    assert to_float(["not a number"]) is None
    assert to_float_all([0.1, "x", 0.2]) == [0.1, 0.2]
    assert to_int([7.9]) == 7


# ==========================================================================
# Degrading when a source does not answer
# ==========================================================================


def test_a_report_is_complete_without_the_second_annotator():
    """Ensembl VEP is capped at 8s and fails outright on a fifth of lookups,
    so running without it is an expected path rather than an edge case.

    It supplies a second opinion on the molecular consequence, not the
    worksheet: snpEff provides the consequence through MyVariant. If losing
    VEP ever changed the classification or the criteria, the cap would be
    trading correctness for speed, which is not the trade being made.
    """
    from dissent.acmg.combine import classify
    from dissent.acmg.criteria import evaluate_all

    bundle = conftest_load_bundle("hfe_c282y")

    with_vep = classify(evaluate_all(bundle))
    assert bundle.consequence is not None, "fixture should carry a consequence"

    bundle.consequence_vep = None
    without_vep = classify(evaluate_all(bundle))

    assert without_vep.call == with_vep.call
    assert without_vep.met_codes == with_vep.met_codes


def test_losing_the_second_annotator_silences_the_comparison_rather_than_faking_it():
    """With one annotator there is nothing to compare, so the detector must
    stay quiet instead of reporting agreement."""
    from dissent.conflicts.detectors import annotator_disagreement

    bundle = conftest_load_bundle("braf_v600e")
    bundle.consequence_vep = None
    assert annotator_disagreement(bundle) is None
