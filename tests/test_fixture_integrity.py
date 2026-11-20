"""The recorded bundles must still say what the current code would derive.

The fixtures are serialized `EvidenceBundle`s, so fields like `significance`
and `stars` are restored as stored values -- `normalize_significance` and
`review_stars` never run on them. That makes every end-to-end test in this
suite an assertion about a mapping that may no longer exist.

It had already happened. Four fields across two fixtures held values the
current mappers do not produce, and one of them mattered: HFE's
"Pathogenic; other" record was stored as `other`, which `is_actionable` drops,
so a pathogenic assertion was missing from the evidence of the variant this
project uses as its headline example.

These tests make a mapping change fail loudly here rather than silently
weaken every test that reads a fixture.
"""

import json
from pathlib import Path

import pytest

from dissent.models import normalize_significance, review_stars

FIXTURES = sorted((Path(__file__).parent / "fixtures").glob("*.json"))


def assertions_in(path: Path):
    data = json.loads(path.read_text())
    return (data.get("clinvar") or {}).get("assertions", [])


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_stored_significance_matches_the_current_mapping(path):
    for a in assertions_in(path):
        raw = a.get("significance_raw")
        if not raw:
            continue
        expected = normalize_significance(raw).value
        assert a.get("significance") == expected, (
            f"{path.name} {a.get('accession')}: raw {raw!r} now maps to "
            f"{expected!r}, fixture stores {a.get('significance')!r}. "
            "Re-derive the fixtures."
        )


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
def test_stored_stars_match_the_current_mapping(path):
    for a in assertions_in(path):
        expected = review_stars(a.get("review_status"))
        assert a.get("stars") == expected, (
            f"{path.name} {a.get('accession')}: review status "
            f"{a.get('review_status')!r} now rates {expected}, fixture stores "
            f"{a.get('stars')}. Re-derive the fixtures."
        )


def test_the_fixtures_exercise_the_compound_label_path():
    """At least one recorded record must carry a semicolon-joined label, or
    the mapping that handles them is only tested in isolation."""
    compound = [
        a.get("significance_raw")
        for p in FIXTURES
        for a in assertions_in(p)
        if ";" in (a.get("significance_raw") or "")
    ]
    assert compound, "no fixture carries a compound ClinVar label"


def test_no_fixture_silently_drops_a_record_that_states_a_classification():
    """The specific regression: a record whose label clearly states a
    classification, stored under a value that `is_actionable` filters out.

    Asserts actionability rather than a particular side. "Pathogenic/
    Pathogenic, low penetrance" correctly keeps its qualifier and is not on
    the pathogenic side -- but it must still reach the evidence, because
    dropping it loses a record entirely.
    """
    from dissent.models import Significance

    stating = ("pathogenic", "likely pathogenic", "benign", "likely benign", "uncertain")
    for p in FIXTURES:
        for a in assertions_in(p):
            raw = (a.get("significance_raw") or "").lower()
            if not raw.startswith(stating):
                continue
            sig = Significance(a.get("significance"))
            assert sig.is_actionable, (
                f"{p.name} {a.get('accession')}: {raw!r} stored as {sig.value!r}, "
                "which is filtered out as not actionable -- the record vanishes "
                "from the evidence"
            )
