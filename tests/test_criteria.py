"""Criteria evaluated against recorded real-world evidence.

Each fixture was chosen because it exercises a distinct branch:

    cftr_g542x   stop_gained in a constrained gene       -> PVS1
    hfe_c282y    pathogenic claims at 7% allele frequency -> BA1 vs ClinVar
    mthfr_a222v  common polymorphism, predictors damaging -> BA1 and PP3 at once
    brca1_c61g   well-curated pathogenic missense         -> PM2, PP3
    braf_v600e   somatic/germline mix                     -> origin handling
"""

import pytest

from dissent.acmg.criteria import ALL_CRITERIA, evaluate_all
from dissent.models import EvidenceBundle, VariantQuery
from tests.conftest import available_bundles, load_bundle

ALL_CODES = {
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
}


def by_code(results):
    return {r.code: r for r in results}


def test_all_28_criteria_are_always_returned():
    """Omitting what cannot be evaluated invites absence to be read as evidence."""
    results = evaluate_all(EvidenceBundle(query=VariantQuery(raw="x")))
    assert len(results) == 28
    assert {r.code for r in results} == ALL_CODES


@pytest.mark.parametrize("name", available_bundles())
def test_every_fixture_evaluates_all_criteria_without_error(name):
    results = evaluate_all(load_bundle(name))
    assert len(results) == 28
    assert {r.code for r in results} == ALL_CODES


@pytest.mark.parametrize("name", available_bundles())
def test_every_result_carries_a_rationale(name):
    """A status with no stated reason is not auditable."""
    for r in evaluate_all(load_bundle(name)):
        assert r.rationale.strip(), f"{name}/{r.code} has an empty rationale"


@pytest.mark.parametrize("name", available_bundles())
def test_unresolved_criteria_say_what_would_settle_them(name):
    for r in evaluate_all(load_bundle(name)):
        if r.status in ("insufficient", "manual"):
            assert r.what_would_change, f"{name}/{r.code} gives no path to resolution"


# ---- PVS1 -----------------------------------------------------------------


def test_pvs1_met_for_nonsense_in_a_constrained_gene():
    nf1 = load_bundle("nf1_r1947ter")
    assert nf1.consequence is not None and nf1.consequence.is_lof
    pvs1 = by_code(evaluate_all(nf1))["PVS1"]
    assert pvs1.status == "met"
    assert pvs1.strength.value == "very_strong"


def test_pvs1_met_branch_says_constraint_is_only_a_proxy():
    """High pLI supports haploinsufficiency; it does not establish mechanism."""
    pvs1 = by_code(evaluate_all(load_bundle("nf1_r1947ter")))["PVS1"]
    assert "proxy" in pvs1.rationale.lower()


def test_pvs1_low_pli_is_not_treated_as_counter_evidence():
    """Recessive genes score as unconstrained because carriers are healthy.

    CFTR has pLI 3.8e-65 and LOEUF 1.19 while loss of function is certainly the
    mechanism of cystic fibrosis. Gating PVS1 on constraint would discard every
    recessive disease gene, so a nonsense variant there stays open, not declined.
    """
    cftr = load_bundle("cftr_g542x")
    assert cftr.consequence is not None and cftr.consequence.is_lof
    assert cftr.constraint is not None and (cftr.constraint.pli or 0) < 0.1

    pvs1 = by_code(evaluate_all(cftr))["PVS1"]
    assert pvs1.status == "insufficient", "must not be declined on constraint alone"
    assert "recessive" in pvs1.rationale.lower()
    assert "carriers are healthy" in pvs1.rationale.lower()


def test_pvs1_not_met_for_missense():
    pvs1 = by_code(evaluate_all(load_bundle("brca1_c61g")))["PVS1"]
    assert pvs1.status == "not_met"


# ---- frequency criteria ---------------------------------------------------


def test_ba1_met_for_a_variant_above_five_percent():
    hfe = load_bundle("hfe_c282y")
    assert hfe.max_af() is not None and hfe.max_af() > 0.05
    ba1 = by_code(evaluate_all(hfe))["BA1"]
    assert ba1.status == "met"
    assert ba1.strength.value == "stand_alone"


def test_ba1_not_met_for_a_rare_variant():
    assert by_code(evaluate_all(load_bundle("brca1_c61g")))["BA1"].status == "not_met"


def test_pm2_met_for_a_rare_variant_and_applied_at_supporting():
    """ClinGen SVI downgrades PM2 from Moderate to Supporting."""
    pm2 = by_code(evaluate_all(load_bundle("brca1_c61g")))["PM2"]
    assert pm2.status == "met"
    assert pm2.strength.value == "supporting"


def test_pm2_not_met_for_a_common_variant():
    assert by_code(evaluate_all(load_bundle("mthfr_a222v")))["PM2"].status == "not_met"


def test_bs1_threshold_is_declared_as_a_default_not_gene_specific():
    bs1 = by_code(evaluate_all(load_bundle("mthfr_a222v")))["BS1"]
    assert bs1.status == "met"
    assert "default" in bs1.rationale.lower()
    assert bs1.what_would_change


def test_bs2_met_where_homozygotes_are_abundant():
    mthfr = load_bundle("mthfr_a222v")
    bs2 = by_code(evaluate_all(mthfr))["BS2"]
    assert bs2.status == "met"
    # The caveat matters: gnomAD is not a healthy cohort and inheritance is unknown.
    assert "inheritance" in bs2.rationale.lower()


# ---- computational criteria ----------------------------------------------


def test_pp3_strength_follows_the_calibrated_predictor():
    """Pejaver 2022 calibration, not a tally of agreeing tools."""
    pp3 = by_code(evaluate_all(load_bundle("brca1_c61g")))["PP3"]
    assert pp3.status == "met"
    assert pp3.strength.value in ("strong", "moderate", "supporting")
    assert "revel" in pp3.rationale.lower()


def test_pp3_and_bp4_are_never_both_met():
    """They are mirror images; satisfying both would be incoherent."""
    for name in available_bundles():
        results = by_code(evaluate_all(load_bundle(name)))
        both = results["PP3"].status == "met" and results["BP4"].status == "met"
        assert not both, f"{name} met both PP3 and BP4"


def test_pp3_not_applied_to_a_nonsense_variant():
    """REVEL is calibrated for missense; applying it elsewhere is unjustified."""
    pp3 = by_code(evaluate_all(load_bundle("cftr_g542x")))["PP3"]
    assert pp3.status == "not_met"
    assert "missense" in pp3.rationale.lower()


# ---- retired criteria -----------------------------------------------------


@pytest.mark.parametrize("code", ["PP5", "BP6"])
@pytest.mark.parametrize("name", available_bundles())
def test_retired_criteria_are_never_applied(code, name):
    """ClinGen retired PP5/BP6: using another lab's conclusion double-counts."""
    result = by_code(evaluate_all(load_bundle(name)))[code]
    assert result.status != "met"
    assert "retired" in result.rationale.lower()


# ---- determinism ----------------------------------------------------------


@pytest.mark.parametrize("name", available_bundles())
def test_evaluation_is_deterministic(name):
    bundle = load_bundle(name)
    first = [(r.code, r.status, r.strength) for r in evaluate_all(bundle)]
    second = [(r.code, r.status, r.strength) for r in evaluate_all(bundle)]
    assert first == second


def test_criteria_are_pure_and_do_not_mutate_the_bundle():
    bundle = load_bundle("braf_v600e")
    before = bundle.model_dump_json()
    for fn in ALL_CRITERIA:
        fn(bundle)
    assert bundle.model_dump_json() == before
