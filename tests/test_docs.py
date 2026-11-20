"""The README has to agree with the code.

Twice now a figure in the documentation has drifted from the thing it
describes: the detector count after detectors were added, and the split
between computable and human-judged ACMG criteria, which the README stated two
different ways in two paragraphs. Both are the kind of error a reader cannot
catch, because the claim looks like a fact about the code and nothing checks
it. These tests make the code the source of the number.

They deliberately assert the *counts*, not the prose. Rewriting a sentence
should not fail a test; changing what the software does, without saying so,
should.
"""

import re
from pathlib import Path

import pytest

from dissent.acmg.criteria import evaluate_all
from dissent.conflicts import detectors as detectors_module
from dissent.panel import TIERS
from tests.conftest import available_bundles, load_bundle

ROOT = Path(__file__).resolve().parent.parent
README = (ROOT / "README.md").read_text()


# --------------------------------------------------------------------------
# Counts taken from the code
# --------------------------------------------------------------------------


def detector_count() -> int:
    """Distinct `kind=` values a detector can emit."""
    src = Path(detectors_module.__file__).read_text()
    return len(set(re.findall(r'kind="([a-z_]+)"', src)))


def criterion_counts() -> tuple[int, int]:
    """(total, needing a human) — the second measured across every bundle.

    A criterion counts as needing a human only if it is marked so for every
    recorded variant. One that is computable for some variants and not others
    belongs on the computable side; calling it manual would overstate what the
    engine cannot do.
    """
    always_manual: set[str] | None = None
    total = 0
    for name in available_bundles():
        results = evaluate_all(load_bundle(name))
        total = len(results)
        manual = {r.code for r in results if r.status == "manual"}
        always_manual = manual if always_manual is None else always_manual & manual
    return total, len(always_manual or set())


# --------------------------------------------------------------------------
# The claims
# --------------------------------------------------------------------------


def test_the_criteria_total_is_twenty_eight():
    """The ACMG/AMP 2015 framework has 28 criteria. If this changes, every
    count below it in the README is wrong too."""
    total, _ = criterion_counts()
    assert total == 28


def test_readme_states_the_right_split_of_computable_criteria():
    total, manual = criterion_counts()
    computable = total - manual
    assert f"{computable} are computed from public data" in README
    assert f"the other {manual} say what evidence would" in README


def test_readme_limits_section_agrees_with_itself():
    """The same figure was stated as 16 in one paragraph and 17 in another."""
    _, manual = criterion_counts()
    assert f"{manual} of 28 criteria need a human" in README


def test_readme_states_the_right_number_of_detectors():
    words = {7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten", 11: "Eleven", 12: "Twelve",
             13: "Thirteen", 14: "Fourteen", 15: "Fifteen"}
    n = detector_count()
    assert words.get(n, str(n)) + " detectors" in README, (
        f"the code has {n} detectors; the README says otherwise"
    )


def test_the_detector_table_lists_every_detector():
    """A detector missing from the table is a feature nobody is told about."""
    n = detector_count()
    # lstrip because the split leaves the newline that ended the header row,
    # so the first element of splitlines() is an empty string.
    after = README.split("| Detector | What it catches |", 1)[1].lstrip("\n")
    rows = []
    for line in after.splitlines():
        if not line.startswith("|"):
            # The table ends at the first line that is not part of it. Reading
            # to the end of the file instead would sweep in every later table
            # in the README and count those rows too.
            break
        if "---" in line:
            continue
        rows.append(line)
    assert len(rows) == n


def test_readme_names_every_panel_tier():
    for tier in TIERS:
        if tier.key == "unresolved":
            continue  # a failure state, not one of the four worklist tiers
        assert f"**{tier.name}**" in README, f"tier {tier.name!r} is not in the README"


def test_the_quoted_test_count_matches_the_suite(pytestconfig):
    """The README quotes a figure for the size of this suite.

    Read from the collector rather than hardcoded, so adding a test cannot
    leave the claim stale -- which is exactly what happened at 308.

    Only meaningful when the whole suite is being collected. Running a subset
    (`pytest tests/test_docs.py`) collects a handful, and failing on that would
    make the test about how it was invoked rather than about the README.
    """
    quoted = {int(m) for m in re.findall(r"Tests — (\d+) of them", README)}
    quoted |= {int(m) for m in re.findall(r"the (\d+) tests run offline", README)}
    if not quoted:
        pytest.skip("the README no longer quotes a test count")

    testpaths = set(pytestconfig.getini("testpaths") or [])
    if set(pytestconfig.args) - testpaths - {str(ROOT)}:
        pytest.skip("a subset of the suite was selected; the count is not comparable")

    reporter = pytestconfig.pluginmanager.get_plugin("terminalreporter")
    collected = getattr(reporter, "_numcollected", 0) if reporter else 0
    if not collected:
        pytest.skip("the collected count is not available in this run mode")

    assert quoted == {collected}, (
        f"README says {sorted(quoted)}, the suite collects {collected}"
    )
