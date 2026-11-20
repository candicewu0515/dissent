"""The measured ClinVar base rates and the lookups built on them.

These figures are quoted in the UI as measurements, so they must be internally
consistent: a rate that disagrees with its own numerator and denominator is
worse than no rate at all.
"""

from __future__ import annotations

import json

import pytest

from dissent.baseline.rates import DATA_PATH, Baseline, GeneRate, load_baseline


@pytest.fixture(scope="module")
def baseline() -> Baseline:
    loaded = load_baseline()
    assert loaded is not None, "baseline data should ship with the repository"
    return loaded


def test_the_data_file_is_present_and_parses():
    assert DATA_PATH.exists()
    json.loads(DATA_PATH.read_text())


def test_absent_data_is_handled_rather_than_fatal(tmp_path):
    """The tool works without the file; it simply omits the comparison."""
    load_baseline.cache_clear()
    try:
        assert load_baseline(tmp_path / "nothing.json") is None
    finally:
        load_baseline.cache_clear()


def test_headline_counts_are_consistent_with_their_rates(baseline):
    counts, rates = baseline.counts, baseline.rates_pct
    pairs = [
        ("any_label_disagreement", "any_label_disagreement"),
        ("pathogenic_vs_benign", "pathogenic_vs_benign"),
        ("every_record_lacks_criteria", "every_record_lacks_criteria"),
    ]
    for count_key, rate_key in pairs:
        derived = 100.0 * counts[count_key] / baseline.variants
        assert derived == pytest.approx(rates[rate_key], abs=0.02), rate_key


def test_every_rate_is_a_percentage(baseline):
    for key, value in baseline.rates_pct.items():
        assert 0.0 <= value <= 100.0, key


def test_the_multi_record_rate_exceeds_the_whole_database_rate(baseline):
    """Variants with one record cannot disagree, so the conditional rate is higher."""
    overall = baseline.rates_pct["any_label_disagreement"]
    among_multi = baseline.rates_pct["any_label_disagreement_among_multi_submitter"]
    assert among_multi > overall


@pytest.mark.parametrize("gene", ["BRCA1", "BRCA2", "TP53", "CFTR", "LDLR"])
def test_well_known_genes_have_a_quotable_rate(gene):
    rate = load_baseline().gene(gene)
    assert rate is not None
    assert rate.variants >= 30
    assert 0 <= rate.discordant <= rate.variants


def test_gene_lookup_is_case_insensitive(baseline):
    assert baseline.gene("brca1") == baseline.gene("BRCA1")


def test_an_unmeasured_gene_returns_nothing_rather_than_zero(baseline):
    """A gene below the variant floor has no rate, which is not a rate of 0%."""
    assert baseline.gene("NOTAREALGENE") is None
    assert baseline.gene(None) is None
    assert baseline.gene_rank("NOTAREALGENE") is None


def test_the_two_variant_floors_are_distinct_and_recorded(baseline):
    """A rate is quotable at a lower bar than a ranking is trustworthy at.

    At 30 variants the standard error on a rate near 30% is roughly nine
    percentage points, enough for a small gene to top a table on noise. The
    leaderboard therefore uses a higher floor, and the two populations are
    never presented as the same one.
    """
    assert baseline.leaderboard_min_variants > baseline.gene_min_variants
    for g in baseline.most_discordant_genes:
        assert g["variants"] >= baseline.leaderboard_min_variants


def test_the_leaderboard_top_need_not_be_rank_one_overall(baseline):
    """It ranks within the larger, noisier population, and says so."""
    top = baseline.most_discordant_genes[0]
    rank = baseline.gene_rank(top["gene"])
    assert rank is not None
    assert rank[1] == len(baseline.gene_table)
    # Ahead of every other gene that clears the leaderboard floor...
    for g in baseline.most_discordant_genes[1:]:
        assert top["rate_pct"] >= g["rate_pct"]
    # ...but small genes below that floor can still rank above it.
    assert rank[0] >= 1


def test_rank_one_really_is_the_highest_rate_in_the_table(baseline):
    best = max(
        baseline.gene_table.items(),
        key=lambda kv: (kv[1][1] / kv[1][0]) if kv[1][0] else 0,
    )[0]
    assert baseline.gene_rank(best)[0] == 1


def test_the_leaderboard_is_ordered_by_rate(baseline):
    rates = [g["rate_pct"] for g in baseline.most_discordant_genes]
    assert rates == sorted(rates, reverse=True)


def test_leaderboard_rates_match_their_own_counts(baseline):
    for g in baseline.most_discordant_genes:
        derived = 100.0 * g["discordant"] / g["variants"]
        assert derived == pytest.approx(g["rate_pct"], abs=0.01), g["gene"]


def test_rank_tracks_rate(baseline):
    """A gene with a higher rate must not rank worse than one with a lower rate."""
    pairs = [("BRCA2", "BRCA1"), ("TP53", "CFTR")]
    for higher, lower in pairs:
        a, b = baseline.gene(higher), baseline.gene(lower)
        assert a and b and a.rate_pct > b.rate_pct, f"{higher} vs {lower}"
        ra, rb = baseline.gene_rank(higher), baseline.gene_rank(lower)
        assert ra[0] < rb[0], f"{higher} should rank ahead of {lower}"


def test_gene_rate_arithmetic():
    assert GeneRate(gene="X", variants=200, discordant=50).rate_pct == pytest.approx(25.0)
    assert GeneRate(gene="X", variants=0, discordant=0).rate_pct == 0.0


def test_the_clinically_consequential_genes_are_worse_than_average(baseline):
    """The finding that gives the number its weight, asserted so it stays true."""
    overall = baseline.rates_pct["any_label_disagreement"]
    for gene in ("BRCA1", "BRCA2", "TP53", "MSH2"):
        rate = baseline.gene(gene)
        assert rate is not None
        assert rate.rate_pct > overall * 2, gene


# ==========================================================================
# The pathogenic-versus-benign list
# ==========================================================================


@pytest.fixture(scope="module")
def opposed() -> list[dict]:
    return json.loads(DATA_PATH.read_text())["opposed_variants"]


def test_the_opposed_count_matches_the_headline_rate(baseline, opposed):
    assert len(opposed) == baseline.counts["pathogenic_vs_benign"]


def test_every_opposed_variant_really_opposes(opposed):
    """Each must carry at least one label from each side, or it does not belong."""
    pathogenic = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
    benign = {"benign", "likely benign", "benign/likely benign"}
    for v in opposed:
        labels = {label for label, _ in v["labels"]}
        assert labels & pathogenic, v["id"]
        assert labels & benign, v["id"]


def test_label_counts_are_consistent_with_the_record_count(opposed):
    for v in opposed:
        assert sum(n for _, n in v["labels"]) == v["records"]


def test_the_list_leads_with_the_best_reviewed_contradictions(opposed):
    """Two expert-panel records disagreeing is a different finding from two
    unsourced opinions doing so, so review status orders the list."""
    stars = [v["best_stars"] for v in opposed]
    assert stars == sorted(stars, reverse=True)
    assert stars[0] >= 3, "some conflicts should carry an expert-panel review"


def test_star_ratings_use_the_full_scale(opposed):
    """A binary has-criteria flag made the ordering meaningless."""
    assert len({v["best_stars"] for v in opposed}) > 2


def test_identifiers_are_usable_as_input(opposed):
    from dissent.normalize.query import parse_variant

    for v in opposed[:25]:
        q = parse_variant(str(v["id"]))
        assert q.kind == "clinvar_id"
        assert q.clinvar_id == v["id"]
