"""The conflict detectors, against recorded real-world evidence.

Detectors must be conservative. One that fires on every variant trains the
reader to ignore the panel, which defeats the point of having one.
"""

from datetime import date, timedelta

import pytest

from dissent.acmg.combine import classify
from dissent.acmg.criteria import evaluate_all
from dissent.conflicts.detectors import (
    detect_all,
    frequency_versus_claim,
    internal_contradiction,
    predictor_split,
    somatic_germline_misuse,
    stale_assertion,
    submitter_discordance,
    weak_evidence_base,
)
from dissent.models import EvidenceBundle, VariantQuery
from tests.conftest import available_bundles, load_bundle


def conflicts_for(name: str, today: date | None = None):
    bundle = load_bundle(name)
    criteria = evaluate_all(bundle)
    classification = classify(criteria)
    return detect_all(bundle, criteria=criteria, classification=classification, today=today)


def kinds(conflicts) -> set[str]:
    return {c.kind for c in conflicts}


# ---- hygiene --------------------------------------------------------------


def test_empty_evidence_produces_no_conflicts():
    """No data must not be reported as disagreement."""
    assert detect_all(EvidenceBundle(query=VariantQuery(raw="x"))) == []


@pytest.mark.parametrize("name", available_bundles())
def test_every_conflict_has_a_headline_and_detail(name):
    for c in conflicts_for(name):
        assert c.headline.strip()
        assert c.detail.strip()
        assert c.severity in ("high", "medium", "low")


@pytest.mark.parametrize("name", available_bundles())
def test_conflicts_are_ordered_by_severity(name):
    order = {"high": 0, "medium": 1, "low": 2}
    severities = [order[c.severity] for c in conflicts_for(name)]
    assert severities == sorted(severities)


@pytest.mark.parametrize("name", available_bundles())
def test_detectors_are_pure(name):
    bundle = load_bundle(name)
    before = bundle.model_dump_json()
    detect_all(bundle, criteria=evaluate_all(bundle))
    assert bundle.model_dump_json() == before


# ---- somatic / germline ---------------------------------------------------


def test_somatic_germline_mix_is_flagged_for_braf(braf):
    """BRAF V600E carries both somatic and germline assertions in ClinVar."""
    assert "somatic" in braf.clinvar.origins
    conflict = somatic_germline_misuse(braf)
    assert conflict is not None
    assert "somatic" in conflict.detail.lower()


def test_germline_only_variant_is_not_flagged_as_somatic(brca1):
    assert "somatic" not in brca1.clinvar.origins
    assert somatic_germline_misuse(brca1) is None


# ---- submitter discordance ------------------------------------------------


def test_discordance_detected_where_submitters_disagree(braf):
    conflict = submitter_discordance(braf)
    assert conflict is not None
    # The spread itself must be shown, not just asserted.
    assert "x pathogenic" in conflict.detail or "x uncertain" in conflict.detail


def test_discordance_reports_dissenting_review_status(braf):
    """A majority is not a consensus if dissenters are equally well reviewed."""
    conflict = submitter_discordance(braf)
    assert conflict is not None and "review star" in conflict.detail


def test_no_discordance_where_submitters_agree(brca1):
    sigs = {a.significance for a in brca1.clinvar.actionable_assertions}
    assert len(sigs) == 1, "fixture should be unanimous"
    assert submitter_discordance(brca1) is None


# ---- frequency versus claim ----------------------------------------------


def test_pathogenic_claim_above_ba1_threshold_is_flagged():
    """HFE C282Y: pathogenic assertions at ~7% allele frequency."""
    hfe = load_bundle("hfe_c282y")
    conflict = frequency_versus_claim(hfe)
    assert conflict is not None
    assert conflict.severity == "high"
    assert "too common" in conflict.headline.lower()


def test_rare_pathogenic_variant_is_not_flagged_on_frequency(brca1):
    assert frequency_versus_claim(brca1) is None


# ---- predictor split ------------------------------------------------------


def test_predictor_split_detected_where_tools_disagree(mthfr):
    conflict = predictor_split(mthfr)
    assert conflict is not None
    assert "split" in conflict.headline.lower()
    # Dissenting tools are named so the reader can weigh them.
    assert "dissenting" in conflict.detail.lower()


def test_predictor_split_needs_enough_tools_to_be_meaningful():
    """Two tools disagreeing is not a finding."""
    cftr = load_bundle("cftr_g542x")
    called = [p for p in cftr.predictors if p.call in ("damaging", "benign")]
    assert len(called) < 6, "fixture should have few predictor calls"
    assert predictor_split(cftr) is None


# ---- weak evidence base ---------------------------------------------------


def test_zero_star_records_are_flagged(braf):
    conflict = weak_evidence_base(braf)
    assert conflict is not None
    assert "assertion criteria" in conflict.detail


def test_expert_panel_review_suppresses_the_weak_evidence_flag():
    """A 3-star record outweighs surrounding unsourced opinions."""
    bundle = load_bundle("braf_v600e")
    for assertion in bundle.clinvar.assertions:
        assertion.stars = 3
    assert weak_evidence_base(bundle) is None


# ---- staleness ------------------------------------------------------------


def test_stale_records_are_flagged_even_under_a_recent_one(braf):
    """A fresh submission must not mark decade-old records current.

    This previously compared only the newest evaluation date across the whole
    variant, so one recent record made everything beneath it look current --
    the opposite of useful, since the risk is an assertion predating gnomAD
    still contributing to the aggregate label. BRAF V600E carries a record
    last evaluated in 2009 and was reported as fresh.
    """
    dated = [a for a in braf.clinvar.actionable_assertions if a.last_evaluated]
    newest = max(a.last_evaluated for a in dated)
    oldest = min(a.last_evaluated for a in dated)
    assert (newest - oldest).days > 365 * 10, "fixture should span over a decade"

    conflict = stale_assertion(braf, today=newest)
    assert conflict is not None, "old records under a recent one must still be flagged"
    assert oldest.isoformat() in conflict.detail
    assert "of" in conflict.headline  # "N of M records ..."


def test_staleness_is_quiet_when_every_record_is_recent():
    """Built synthetically: no real fixture has uniformly recent records.

    BRCA1 C61G, the best-curated variant in the set, still carries a 2015
    evaluation -- which is itself the 3-star expert-panel record.
    """
    from dissent.models import ClinVarAssertion, ClinVarEvidence, Significance

    today = date(2026, 1, 1)
    bundle = EvidenceBundle(
        query=VariantQuery(raw="x"),
        clinvar=ClinVarEvidence(
            variant_id=1,
            assertions=[
                ClinVarAssertion(
                    significance=Significance.PATHOGENIC,
                    stars=2,
                    last_evaluated=date(2025, 6, 1),
                ),
                ClinVarAssertion(
                    significance=Significance.PATHOGENIC,
                    stars=1,
                    last_evaluated=date(2024, 3, 1),
                ),
            ],
        ),
    )
    assert stale_assertion(bundle, today=today) is None


def test_an_outdated_authoritative_record_is_called_out(brca1):
    """The expert-panel review being the oldest is the point, not a detail."""
    dated = [a for a in brca1.clinvar.actionable_assertions if a.last_evaluated]
    best = max(a.stars for a in dated)
    assert best >= 2
    newest_overall = max(a.last_evaluated for a in dated)
    newest_best = max(a.last_evaluated for a in dated if a.stars == best)
    assert newest_best < newest_overall, "fixture premise: best review is not the newest"

    conflict = stale_assertion(brca1, today=newest_overall)
    assert conflict is not None
    assert "best-reviewed record" in conflict.detail


def test_all_records_stale_is_more_severe_than_some(braf):
    dated = [a for a in braf.clinvar.actionable_assertions if a.last_evaluated]
    newest = max(a.last_evaluated for a in dated)
    some = stale_assertion(braf, today=newest)
    everything = stale_assertion(braf, today=date(newest.year + 20, 1, 1))
    assert some is not None and everything is not None
    assert some.severity == "low"
    assert everything.severity == "medium"


# ---- internal contradiction ----------------------------------------------


def test_contradiction_flagged_when_frequency_opposes_predictors(mthfr):
    """MTHFR A222V: 49% allele frequency, yet REVEL calls it damaging."""
    criteria = evaluate_all(mthfr)
    conflict = internal_contradiction(mthfr, criteria)
    assert conflict is not None
    assert "PP3" in conflict.headline and "BA1" in conflict.headline


def test_contradiction_explains_why_the_two_are_not_equal_evidence(mthfr):
    conflict = internal_contradiction(mthfr, evaluate_all(mthfr))
    assert conflict is not None
    assert "direct measurement" in conflict.detail


def test_no_contradiction_when_evidence_points_one_way(brca1):
    assert internal_contradiction(brca1, evaluate_all(brca1)) is None


# ---- the panel as a whole -------------------------------------------------


def test_mthfr_surfaces_the_expected_set(mthfr):
    found = kinds(conflicts_for("mthfr_a222v"))
    assert "internal_contradiction" in found
    assert "predictor_split" in found


def test_hfe_surfaces_the_frequency_contradiction():
    assert "frequency_versus_claim" in kinds(conflicts_for("hfe_c282y"))


def test_well_curated_variant_has_a_quiet_panel(brca1):
    """Unanimous, well-reviewed, rare, predictors agreeing: little to report."""
    found = conflicts_for("brca1_c61g")
    high = [c for c in found if c.severity == "high"]
    assert not high, f"unexpected high-severity conflicts: {[c.kind for c in high]}"
