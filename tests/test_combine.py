"""ACMG/AMP combining rules, Richards 2015 Table 5.

Counting is by applied strength, not each criterion's default, because PP3 can
reach Strong under the Pejaver calibration and PM2 is applied at Supporting.
"""

import pytest

from dissent.acmg.combine import (
    BENIGN,
    CONFLICTING,
    LIKELY_BENIGN,
    LIKELY_PATHOGENIC,
    PATHOGENIC,
    UNCERTAIN,
    classify,
)
from dissent.models import CriterionResult, Strength

VS, S, M, P, SA = (
    Strength.VERY_STRONG,
    Strength.STRONG,
    Strength.MODERATE,
    Strength.SUPPORTING,
    Strength.STAND_ALONE,
)


def crit(code, direction, strength, status="met"):
    return CriterionResult(
        code=code, name=code, direction=direction, status=status,
        strength=strength, rationale="",
    )


def path(code, strength):
    return crit(code, "pathogenic", strength)


def ben(code, strength):
    return crit(code, "benign", strength)


@pytest.mark.parametrize(
    "criteria,expected,label",
    [
        ([path("PVS1", VS), path("PS1", S)], PATHOGENIC, "1 very strong + 1 strong"),
        ([path("PVS1", VS), path("PM1", M), path("PM4", M)], PATHOGENIC, "1 VS + 2 moderate"),
        ([path("PVS1", VS), path("PM4", M), path("PP3", P)], PATHOGENIC, "1 VS + 1 mod + 1 supp"),
        ([path("PVS1", VS), path("PP2", P), path("PP3", P)], PATHOGENIC, "1 VS + 2 supporting"),
        ([path("PS1", S), path("PS3", S)], PATHOGENIC, "2 strong"),
        ([path("PS1", S)] + [path(f"PM{i}", M) for i in (1, 3, 4)], PATHOGENIC, "1 strong + 3 mod"),
        ([path("PVS1", VS), path("PM4", M)], LIKELY_PATHOGENIC, "1 VS + 1 moderate"),
        ([path("PS1", S), path("PM4", M)], LIKELY_PATHOGENIC, "1 strong + 1 moderate"),
        ([path("PS1", S), path("PP2", P), path("PP3", P)], LIKELY_PATHOGENIC, "1 strong + 2 supp"),
        ([path(f"PM{i}", M) for i in (1, 3, 4)], LIKELY_PATHOGENIC, "3 moderate"),
        ([ben("BA1", SA)], BENIGN, "BA1 stand-alone"),
        ([ben("BS1", S), ben("BS2", S)], BENIGN, "2 benign strong"),
        ([ben("BS1", S), ben("BP4", P)], LIKELY_BENIGN, "1 benign strong + 1 supp"),
        ([ben("BP4", P), ben("BP7", P)], LIKELY_BENIGN, "2 benign supporting"),
        ([path("PM2", P)], UNCERTAIN, "1 supporting only"),
        ([], UNCERTAIN, "no criteria at all"),
    ],
)
def test_combining_rules(criteria, expected, label):
    assert classify(criteria).call == expected, label


def test_contradictory_evidence_yields_uncertain():
    """Richards 2015: contradictory criteria give uncertain significance."""
    result = classify([ben("BA1", SA), path("PVS1", VS), path("PS1", S)])
    assert result.call == CONFLICTING
    assert "contradictory" in result.reasoning.lower()


def test_not_met_criteria_do_not_count():
    criteria = [
        crit("PVS1", "pathogenic", VS, status="not_met"),
        crit("PS1", "pathogenic", S, status="insufficient"),
        crit("PM1", "pathogenic", M, status="manual"),
    ]
    assert classify(criteria).call == UNCERTAIN


def test_reasoning_names_the_rule_that_fired():
    """The user must be able to check the arithmetic, not just read a label."""
    result = classify([path("PVS1", VS), path("PS1", S)])
    # The rule that fired, verbatim, not merely words that appear in it.
    assert "1 very strong + >=1 strong" in result.reasoning
    assert "1 very strong, 1 strong" in result.reasoning


@pytest.mark.parametrize(
    "criteria,label",
    [
        ([path("PS1", S), path("PP3", P)], "1 strong + 1 supporting"),
        ([path("PM1", M), path("PM4", M)], "2 moderate"),
        ([path("PM1", M)] + [path(f"PP{i}", P) for i in (1, 2, 3)], "1 mod + 3 supp"),
        ([path("PP1", P)], "1 supporting"),
        ([ben("BS1", S)], "1 benign strong"),
        ([ben("BP4", P)], "1 benign supporting"),
        ([path("PVS1", VS)], "1 very strong alone"),
    ],
)
def test_sub_threshold_combinations_stay_uncertain(criteria, label):
    """An over-permissive rule must fail here, not slip through unnoticed."""
    assert classify(criteria).call == UNCERTAIN, label


@pytest.mark.parametrize(
    "criteria,expected,label",
    [
        (
            [path("PS1", S), path("PM1", M), path("PM4", M), path("PP2", P), path("PP3", P)],
            PATHOGENIC,
            "1 strong + 2 moderate + 2 supporting",
        ),
        (
            [path("PS1", S), path("PM1", M)] + [path(f"PP{i}", P) for i in (1, 2, 3, 4)],
            PATHOGENIC,
            "1 strong + 1 moderate + 4 supporting",
        ),
        (
            [path("PM1", M)] + [path(f"PP{i}", P) for i in (1, 2, 3, 4)],
            LIKELY_PATHOGENIC,
            "1 moderate + 4 supporting",
        ),
        (
            [path("PM1", M), path("PM4", M), path("PP2", P), path("PP3", P)],
            LIKELY_PATHOGENIC,
            "2 moderate + 2 supporting",
        ),
    ],
)
def test_the_remaining_table_5_rows(criteria, expected, label):
    assert classify(criteria).call == expected, label


def test_benign_takes_precedence_over_a_pathogenic_rule():
    """BA1 is stand-alone; it decides even against met pathogenic criteria."""
    result = classify([ben("BA1", SA), path("PM2", P), path("PP3", M)])
    assert result.call == BENIGN


def test_the_framework_gap_is_reported_not_filled():
    """One very strong plus one supporting has no row in Table 5."""
    result = classify([path("PVS1", VS), path("PM2", P)])
    assert result.call == UNCERTAIN
    assert "no row in Richards 2015 Table 5" in result.reasoning
    assert "PM2" in result.reasoning


def test_the_gap_note_appears_only_for_that_combination():
    assert "Table 5" not in classify([path("PM2", P)]).reasoning
    assert "Table 5" not in classify([path("PVS1", VS), path("PM1", M)]).reasoning


def test_met_codes_are_reported():
    result = classify([path("PVS1", VS), path("PS1", S), crit("PM2", "pathogenic", P, "not_met")])
    assert set(result.met_codes) == {"PVS1", "PS1"}
