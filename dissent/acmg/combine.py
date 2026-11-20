"""Combine met criteria into a classification, per Richards 2015 Table 5.

Counting is by the strength actually *applied*, not each criterion's default.
That matters here: PP3 can reach Strong or Moderate under the Pejaver 2022
calibration, and PM2 is applied at Supporting on ClinGen SVI advice. A tool
that counted defaults would reach a different answer from the same evidence.

The returned reasoning names the specific rule that fired, so the user can
check the arithmetic rather than take the label on trust.
"""

from __future__ import annotations

from dataclasses import dataclass

from dissent.models import Classification, CriterionResult, Strength

PATHOGENIC = "Pathogenic"
LIKELY_PATHOGENIC = "Likely pathogenic"
UNCERTAIN = "Uncertain significance"
LIKELY_BENIGN = "Likely benign"
BENIGN = "Benign"
CONFLICTING = "Uncertain significance (conflicting evidence)"


@dataclass(frozen=True)
class Tally:
    very_strong: int = 0
    strong: int = 0
    moderate: int = 0
    supporting: int = 0
    stand_alone: int = 0

    @property
    def any(self) -> bool:
        return bool(
            self.very_strong or self.strong or self.moderate
            or self.supporting or self.stand_alone
        )


def _tally(criteria: list[CriterionResult], direction: str) -> Tally:
    counts = {s: 0 for s in Strength}
    for c in criteria:
        if c.status == "met" and c.direction == direction and c.strength is not None:
            counts[c.strength] += 1
    return Tally(
        very_strong=counts[Strength.VERY_STRONG],
        strong=counts[Strength.STRONG],
        moderate=counts[Strength.MODERATE],
        supporting=counts[Strength.SUPPORTING],
        stand_alone=counts[Strength.STAND_ALONE],
    )


def _pathogenic_rule(t: Tally) -> str | None:
    if t.very_strong >= 1:
        if t.strong >= 1:
            return "1 very strong + >=1 strong"
        if t.moderate >= 2:
            return "1 very strong + >=2 moderate"
        if t.moderate == 1 and t.supporting >= 1:
            return "1 very strong + 1 moderate + >=1 supporting"
        if t.supporting >= 2:
            return "1 very strong + >=2 supporting"
    if t.strong >= 2:
        return ">=2 strong"
    if t.strong == 1:
        if t.moderate >= 3:
            return "1 strong + >=3 moderate"
        if t.moderate == 2 and t.supporting >= 2:
            return "1 strong + 2 moderate + >=2 supporting"
        if t.moderate == 1 and t.supporting >= 4:
            return "1 strong + 1 moderate + >=4 supporting"
    return None


def _likely_pathogenic_rule(t: Tally) -> str | None:
    if t.very_strong >= 1 and t.moderate >= 1:
        return "1 very strong + 1 moderate"
    if t.strong == 1 and 1 <= t.moderate <= 2:
        return "1 strong + 1-2 moderate"
    if t.strong == 1 and t.supporting >= 2:
        return "1 strong + >=2 supporting"
    if t.moderate >= 3:
        return ">=3 moderate"
    if t.moderate == 2 and t.supporting >= 2:
        return "2 moderate + >=2 supporting"
    if t.moderate == 1 and t.supporting >= 4:
        return "1 moderate + >=4 supporting"
    return None


def _benign_rule(t: Tally) -> str | None:
    if t.stand_alone >= 1:
        return "1 stand-alone (BA1)"
    if t.strong >= 2:
        return ">=2 strong"
    return None


def _likely_benign_rule(t: Tally) -> str | None:
    if t.strong == 1 and t.supporting >= 1:
        return "1 strong + >=1 supporting"
    if t.supporting >= 2:
        return ">=2 supporting"
    return None


def _describe(t: Tally) -> str:
    bits = []
    for label, n in (
        ("very strong", t.very_strong),
        ("stand-alone", t.stand_alone),
        ("strong", t.strong),
        ("moderate", t.moderate),
        ("supporting", t.supporting),
    ):
        if n:
            bits.append(f"{n} {label}")
    return ", ".join(bits) if bits else "nothing"


def classify(criteria: list[CriterionResult]) -> Classification:
    path = _tally(criteria, "pathogenic")
    ben = _tally(criteria, "benign")

    met_codes = [c.code for c in criteria if c.status == "met"]

    path_rule = _pathogenic_rule(path)
    lp_rule = _likely_pathogenic_rule(path)
    ben_rule = _benign_rule(ben)
    lb_rule = _likely_benign_rule(ben)

    reached_pathogenic = path_rule or lp_rule
    reached_benign = ben_rule or lb_rule

    # Richards 2015: where pathogenic and benign criteria are both satisfied,
    # the variant is uncertain. Resolving the contradiction is a judgement call
    # about which evidence is stronger, and that is the user's to make.
    if reached_pathogenic and reached_benign:
        return Classification(
            call=CONFLICTING,
            reasoning=(
                f"Pathogenic criteria met ({_describe(path)}, satisfying "
                f"'{reached_pathogenic}') and benign criteria met ({_describe(ben)}, "
                f"satisfying '{reached_benign}'). Richards 2015 directs that "
                "contradictory evidence yields uncertain significance. The conflict "
                "panel below shows what is pulling in each direction."
            ),
            met_codes=met_codes,
        )

    if ben_rule:
        return Classification(
            call=BENIGN,
            reasoning=f"Benign: {_describe(ben)} benign criteria, satisfying '{ben_rule}'.",
            met_codes=met_codes,
        )
    if path_rule:
        return Classification(
            call=PATHOGENIC,
            reasoning=f"Pathogenic: {_describe(path)} pathogenic criteria, satisfying '{path_rule}'.",
            met_codes=met_codes,
        )
    if lp_rule:
        return Classification(
            call=LIKELY_PATHOGENIC,
            reasoning=f"Likely pathogenic: {_describe(path)}, satisfying '{lp_rule}'.",
            met_codes=met_codes,
        )
    if lb_rule:
        return Classification(
            call=LIKELY_BENIGN,
            reasoning=f"Likely benign: {_describe(ben)}, satisfying '{lb_rule}'.",
            met_codes=met_codes,
        )

    note = _framework_gap_note(path, met_codes)
    return Classification(
        call=UNCERTAIN,
        reasoning=(
            f"No combining rule is satisfied. Pathogenic evidence: {_describe(path)}. "
            f"Benign evidence: {_describe(ben)}. This is the expected outcome for most "
            "variants from public data alone; the criteria needing a human are listed "
            f"in the worksheet.{note}"
        ),
        met_codes=met_codes,
    )


def _framework_gap_note(path: Tally, met_codes: list[str]) -> str:
    """Flag a known gap rather than quietly inventing a rule to fill it.

    Richards 2015 Table 5 has no row for one very strong criterion plus a single
    supporting one. That combination did not arise in 2015 because PM2 was a
    Moderate criterion, and PVS1 + PM2 therefore reached Likely Pathogenic via
    'one very strong + one moderate'. ClinGen SVI later recommended applying PM2
    at Supporting strength, which drops this evidence into a gap the published
    table does not cover.

    Many laboratories do call PVS1 + PM2_Supporting likely pathogenic. This tool
    does not silently adopt that convention, because a worksheet whose
    arithmetic cannot be traced to a citable rule is worse than one that says
    where the rule runs out.
    """
    if path.very_strong >= 1 and path.moderate == 0 and path.supporting == 1:
        return (
            " Note: one very strong criterion plus a single supporting one has no row "
            "in Richards 2015 Table 5. The gap is an artefact of ClinGen SVI's later "
            "downgrade of PM2 from Moderate to Supporting -- under the original table "
            "this evidence would have reached Likely Pathogenic. Many laboratories do "
            "apply that combination; this worksheet does not, because the rule is not "
            "in the published framework. Worth a human decision."
        )
    return ""
