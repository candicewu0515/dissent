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


@pytest.mark.parametrize(
    "term,strength",
    [
        ("stop_gained", "very_strong"),
        ("frameshift_variant", "very_strong"),
        ("splice_donor_variant", "very_strong"),
        ("splice_acceptor_variant", "very_strong"),
        # Abou Tayoun 2018: translation may reinitiate at a downstream ATG, so
        # initiation-codon loss is not reliably a null allele.
        ("start_lost", "moderate"),
    ],
)
def test_pvs1_strength_depends_on_the_kind_of_null_variant(term, strength):
    from dissent.models import Consequence, EvidenceBundle, GeneConstraint, VariantQuery

    bundle = EvidenceBundle(
        query=VariantQuery(raw="x"),
        gene_symbol="TESTGENE",
        consequence=Consequence(most_severe=term, terms=[term]),
        constraint=GeneConstraint(pli=0.99, oe_lof_upper=0.2),
    )
    pvs1 = by_code(evaluate_all(bundle))["PVS1"]
    assert pvs1.status == "met"
    assert pvs1.strength.value == strength


def test_pvs1_names_the_rest_of_the_decision_tree_it_did_not_run():
    pvs1 = by_code(evaluate_all(load_bundle("nf1_r1947ter")))["PVS1"]
    assert "nonsense-mediated decay" in pvs1.what_would_change


def test_pvs1_not_met_for_missense():
    pvs1 = by_code(evaluate_all(load_bundle("brca1_c61g")))["PVS1"]
    assert pvs1.status == "not_met"


# ---- frequency criteria ---------------------------------------------------


def _bundle_with_frequency(af: float, an: int = 800_000, hom: int = 0):
    """A minimal bundle carrying one population-frequency record."""
    from dissent.models import EvidenceBundle, PopulationFrequency, VariantQuery

    return EvidenceBundle(
        query=VariantQuery(raw="x"),
        frequencies=[
            PopulationFrequency(
                source="test cohort",
                af=af,
                an=an,
                ac=int(af * an),
                homozygote_count=hom,
                popmax_af=af,
                popmax_population="Test population",
                per_population={"Test population": af},
            )
        ],
    )


@pytest.mark.parametrize(
    "af,ba1,bs1,pm2",
    [
        (0.200, "met", "met", "not_met"),
        (0.050, "met", "met", "not_met"),      # BA1 boundary, inclusive
        (0.049, "not_met", "met", "not_met"),
        (0.010, "not_met", "met", "not_met"),  # BS1 boundary, inclusive
        (0.009, "not_met", "not_met", "not_met"),
        (0.0001, "not_met", "not_met", "met"), # PM2 boundary, inclusive
        (0.00009, "not_met", "not_met", "met"),
        (0.0, "not_met", "not_met", "met"),
    ],
)
def test_frequency_thresholds_at_and_around_each_boundary(af, ba1, bs1, pm2):
    results = by_code(evaluate_all(_bundle_with_frequency(af)))
    assert results["BA1"].status == ba1, f"BA1 at {af}"
    assert results["BS1"].status == bs1, f"BS1 at {af}"
    assert results["PM2"].status == pm2, f"PM2 at {af}"


def test_frequency_criteria_refuse_an_underpowered_cohort():
    """An 8% estimate from 50 alleles must not reach stand-alone BA1."""
    results = by_code(evaluate_all(_bundle_with_frequency(0.08, an=50)))
    for code in ("BA1", "BS1", "PM2"):
        assert results[code].status == "insufficient", code


def test_ba1_met_for_a_variant_above_five_percent():
    hfe = load_bundle("hfe_c282y")
    ba1 = by_code(evaluate_all(hfe))["BA1"]
    assert ba1.status == "met"
    assert ba1.strength.value == "stand_alone"
    # The cohort the number came from must be named in the rationale.
    assert "AN" in ba1.rationale


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


def test_bs2_is_never_scored_only_reported():
    """BS2 needs inheritance mode, age of onset and phenotype; gnomAD has none.

    Scoring it from a homozygote count alone inverts the answer on HFE C282Y,
    whose gnomAD homozygotes are the genotype that causes haemochromatosis.
    """
    for name in available_bundles():
        bs2 = by_code(evaluate_all(load_bundle(name)))["BS2"]
        assert bs2.status == "manual", f"{name}: BS2 must not be scored"
        assert bs2.what_would_change


def test_bs2_still_reports_the_homozygote_count_as_context():
    mthfr = load_bundle("mthfr_a222v")
    assert any(f.homozygote_count for f in mthfr.frequencies)
    bs2 = by_code(evaluate_all(mthfr))["BS2"]
    assert "homozygotes" in bs2.rationale


def test_frequency_data_alone_cannot_reach_benign_via_two_strong_criteria():
    """BS1 and BS2 both read gnomAD, so scoring both double-counted one row.

    Anything above 1% allele frequency has homozygotes, so the pair fired
    together and '>=2 strong' reached the terminal Benign call from a single
    frequency observation.
    """
    from dissent.acmg.combine import classify

    hfe = load_bundle("hfe_c282y")
    criteria = evaluate_all(hfe)
    strong_benign = [
        c for c in criteria
        if c.status == "met" and c.direction == "benign" and c.strength.value == "strong"
    ]
    assert len(strong_benign) <= 1, [c.code for c in strong_benign]
    # BA1 is stand-alone and still carries the call on its own.
    assert classify(criteria).call == "Benign"


# ---- computational criteria ----------------------------------------------


def _bundle_with_revel(score: float):
    """A minimal missense bundle carrying one REVEL score."""
    from dissent.models import Consequence, EvidenceBundle, PredictorResult, VariantQuery

    return EvidenceBundle(
        query=VariantQuery(raw="x"),
        consequence=Consequence(
            most_severe="missense_variant", terms=["missense_variant"]
        ),
        predictors=[PredictorResult(name="REVEL", score=score, call="damaging")],
    )


@pytest.mark.parametrize(
    "revel,status,strength",
    [
        # Pejaver 2022 calibrated thresholds, tested at and either side of each step.
        (0.999, "met", "strong"),
        (0.932, "met", "strong"),
        (0.931, "met", "moderate"),
        (0.773, "met", "moderate"),
        (0.772, "met", "supporting"),
        (0.644, "met", "supporting"),
        (0.643, "not_met", None),
        (0.100, "not_met", None),
    ],
)
def test_pp3_strength_steps_at_the_calibrated_thresholds(revel, status, strength):
    """Each rung of the ladder, so inverting or shifting it fails the suite."""
    pp3 = by_code(evaluate_all(_bundle_with_revel(revel)))["PP3"]
    assert pp3.status == status, f"REVEL {revel}"
    assert (pp3.strength.value if pp3.status == "met" else None) == strength


@pytest.mark.parametrize(
    "revel,status,strength",
    [
        (0.016, "met", "strong"),
        (0.017, "met", "moderate"),
        (0.183, "met", "moderate"),
        (0.184, "met", "supporting"),
        (0.290, "met", "supporting"),
        (0.291, "not_met", None),
        (0.900, "not_met", None),
    ],
)
def test_bp4_strength_steps_at_the_calibrated_thresholds(revel, status, strength):
    bp4 = by_code(evaluate_all(_bundle_with_revel(revel)))["BP4"]
    assert bp4.status == status, f"REVEL {revel}"
    assert (bp4.strength.value if bp4.status == "met" else None) == strength


@pytest.mark.parametrize("revel", [0.0, 0.1, 0.29, 0.3, 0.5, 0.644, 0.9, 1.0])
def test_pp3_and_bp4_are_never_both_met_across_the_whole_range(revel):
    """Mirror criteria: satisfying both would be incoherent at any score."""
    results = by_code(evaluate_all(_bundle_with_revel(revel)))
    assert not (results["PP3"].status == "met" and results["BP4"].status == "met")


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
