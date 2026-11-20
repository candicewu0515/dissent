"""Panel triage: interpret a list of variants and order them by what needs a human.

A single lookup saves a reader fifteen minutes of tab-switching. A panel
changes which variants get read at all. A gene panel returns dozens of
candidates and whoever reviews them has to choose where to start; that choice
is usually made on gut feeling, or on whatever order the sequencer's output
happened to arrive in. This module makes it explicit.

Two rules govern the ordering, and both exist to keep it auditable:

1. A variant's tier is computed from structured fields only -- ClinVar's own
   assertion labels, their review status and evaluation dates, and this tool's
   computed call. Never from the prose of a conflict message, which is written
   for a human and would make the ranking depend on wording.
2. There is no score. A number like "7.4" invites the reader to trust a
   precision that does not exist here. Four named tiers, each with the reasons
   that put the variant in it, can be disagreed with -- which is the point.
"""

from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from datetime import date, datetime, timezone

import httpx
from pydantic import BaseModel, Field

from dissent.acmg.combine import (
    BENIGN,
    CONFLICTING,
    LIKELY_BENIGN,
    LIKELY_PATHOGENIC,
    PATHOGENIC,
)
from dissent.acmg.thresholds import DEFAULT, Thresholds
from dissent.cache.store import ResponseCache
from dissent.models import UNKNOWN_REVIEW_STARS, Report, Significance
from dissent.report import interpret

# How many variants are interpreted at once. Every one of them fans out to six
# public APIs, none of which the project pays for, so this is a politeness
# limit as much as a performance one: four in flight is roughly 24 concurrent
# requests, which the sources absorb without rate-limiting. Raising it makes
# the whole panel slower, not faster, once a source starts shedding load.
MAX_CONCURRENCY = 4

# A cap per request. The web UI chunks a long panel into fives so the reader
# sees rows arrive, so this ceiling is reached only by other callers -- the
# CLI and anything using the API directly. Whatever it drops is counted and
# reported rather than silently discarded.
MAX_PANEL_SIZE = 50

# Three stars is ClinVar's "reviewed by expert panel"; four is "practice
# guideline". At or above it, a record carrying no assertion criteria is noise
# beside one a panel has signed off, and saying so would push well-curated
# variants into the thin-evidence tier. Below it there is nothing to outweigh
# the gap. Named here because it decides a tier and was previously a bare
# literal in the middle of the rule.
EXPERT_REVIEW_STARS = 3

# Conflict kinds the tier rule already accounts for through its own facts, so
# finding one is not new information. Everything else a detector can raise is
# a fact the tier facts do not encode -- a high-severity one of those must
# therefore be able to lift a row out of the bottom tier, or "No disagreement
# found" is a name the rule cannot keep.
#
# Matched on `kind`, which is a fixed slug on the Conflict model, never on the
# headline: the rule must not depend on wording written for a human.
TIER_ACCOUNTS_FOR = frozenset(
    {
        "record_discordance",
        "computed_versus_reported",
        "weak_evidence_base",
        "stale_assertion",
        "unrecognised_review_status",
    }
)

# Short, checkable names for the rest. A reason has to say what was found, and
# a slug in a worklist is not that.
CONFLICT_LABELS = {
    "annotator_disagreement": "two annotators disagree about the molecular consequence",
    "internal_contradiction": "pathogenic and benign criteria are both met",
    "frequency_versus_claim": "the population frequency contradicts the reported claim",
    "somatic_germline_misuse": "a somatic variant is being read against germline criteria",
    "predictor_split": "the prediction tools are split",
    "ambiguous_input": "the identifier covers more than one allele",
    "unreviewed_literature": "substantial literature sits above weak curation",
}


# --------------------------------------------------------------------------
# Tiers
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Tier:
    rank: int
    key: str
    name: str
    meaning: str


CONTRADICTED = Tier(
    1,
    "contradicted",
    "Contradicted",
    "Two parties land on opposite sides of the clinical line: either the "
    "submitted records disagree pathogenic-versus-benign, or this tool's own "
    "call opposes theirs. One of the two is wrong and a report may already "
    "have gone out.",
)
SPLIT = Tier(
    2,
    "split",
    "Records disagree",
    "ClinVar's records for this variant classify it differently, but not "
    "across the clinical line -- pathogenic against uncertain, typically. "
    "The aggregate label a report quotes picks one of them. Where the records "
    "cover different conditions, part of that spread may be different "
    "clinical questions rather than disagreement, and the row says which.",
)
THIN = Tier(
    3,
    "thin",
    "Thin evidence",
    "Nobody disagrees, but the agreement rests on little: no stated criteria, "
    "one star of four, or an evaluation older than the thresholds and "
    "predictors now in use.",
)
CONCORDANT = Tier(
    4,
    "concordant",
    "No disagreement found",
    "The submitted records agree with each other, state their criteria, and "
    "are current. Not the same as 'correct' -- where the computed call is "
    "more cautious than theirs, that is flagged on the row instead.",
)
UNRESOLVED = Tier(
    5,
    "unresolved",
    "Unresolved",
    "The input did not resolve to a variant, or a lookup failed. Read as "
    "'not yet triaged', never as 'nothing found'.",
)

TIERS = (CONTRADICTED, SPLIT, THIN, CONCORDANT, UNRESOLVED)
TIERS_BY_KEY = {t.key: t for t in TIERS}

# Why "computed uncertain, reported pathogenic" is not a tier.
#
# It was, in the first version of this module, and it swallowed 21 of 30 rows
# in a test panel -- which is the correct number and makes the tier useless.
# A clinical laboratory classifies from segregation in the family, functional
# assays and internal case series. None of that is in a public API, so a tool
# built only on public data reaching "uncertain" where a lab reached
# "pathogenic" is the expected outcome, not a finding. Ranking a worklist by
# it would hand the reader the entire panel and call it a priority list.
#
# It is still worth saying, so it is a flag on the row and a count in the
# summary, labelled as normal.
EVIDENCE_NOT_PUBLIC = (
    "Reported call is confident, computed call is uncertain: the submitter's "
    "evidence is not in public data. Expected, not a flag."
)


# The five labels that make a claim about pathogenicity. Everything else a
# submitter can send -- drug response, risk factor, low penetrance, an
# aggregate "conflicting" label that is not any one submitter's opinion -- is
# an answer to a different question and cannot disagree with these.
AXIS_LABELS = frozenset(
    {
        Significance.PATHOGENIC,
        Significance.LIKELY_PATHOGENIC,
        Significance.UNCERTAIN,
        Significance.LIKELY_BENIGN,
        Significance.BENIGN,
    }
)


def _call_side(call: str | None) -> str:
    """Which side of the clinical line a classification sits on."""
    if not call:
        return "none"
    if call in (PATHOGENIC, LIKELY_PATHOGENIC):
        return "pathogenic"
    if call in (BENIGN, LIKELY_BENIGN):
        return "benign"
    return "uncertain"


# --------------------------------------------------------------------------
# Row
# --------------------------------------------------------------------------


class PanelRow(BaseModel):
    """One variant's line in the worklist."""

    query: str
    resolved: bool = False
    display_name: str | None = None
    gene: str | None = None
    computed_call: str | None = None
    reported_call: str | None = None
    best_stars: int | None = None
    records: int = 0
    conflicts: int = 0
    high_severity: int = 0
    conflict_headlines: list[str] = Field(default_factory=list)
    met_codes: list[str] = Field(default_factory=list)
    tier: str = UNRESOLVED.key
    tier_rank: int = UNRESOLVED.rank
    # True when the reported call is confident and the computed one is not.
    # Deliberately not a tier; see EVIDENCE_NOT_PUBLIC above.
    evidence_not_public: bool = False
    # Why this row is in this tier. Short, plural, and the only justification
    # the ordering offers -- so it has to be complete.
    reasons: list[str] = Field(default_factory=list)
    error: str | None = None
    elapsed_ms: int | None = None


class PanelResult(BaseModel):
    rows: list[PanelRow] = Field(default_factory=list)
    tier_counts: dict[str, int] = Field(default_factory=dict)
    evidence_not_public: int = 0
    # Entries the parser refused: over the per-request cap, or too long to
    # be an identifier. Reported so a truncated worklist cannot be read as
    # a complete one.
    dropped: int = 0
    generated_at: str | None = None
    elapsed_ms: int | None = None


# --------------------------------------------------------------------------
# Triage
# --------------------------------------------------------------------------


@dataclass
class _Facts:
    """The structured facts the tier is computed from, kept separate so the
    rule can be read without the model plumbing around it."""

    computed_side: str = "none"
    reported_sides: set[str] = field(default_factory=set)
    reported_label: str | None = None
    best_stars: int | None = None
    records: int = 0
    axis_records: int = 0
    off_axis_records: int = 0
    no_criteria: int = 0
    rated_records: int = 0
    # A ClinVar record is one aggregate per condition, so disagreement
    # inside a condition and disagreement across conditions are different
    # findings and the row has to say which it is.
    axis_conditions: int = 0
    split_within_condition: bool = False
    stale_years: float | None = None
    opposed_records: bool = False
    # High-severity findings the tier facts do not otherwise encode, by
    # `kind`. Carried so the bottom tier cannot claim there is nothing to
    # read while a detector is saying there is.
    unaccounted_high: list[str] = field(default_factory=list)
    off_axis_labels: list[str] = field(default_factory=list)
    # The engine satisfied pathogenic and benign criteria at once. Carried
    # separately because combine.py returns it as a flavour of "Uncertain
    # significance", where _call_side cannot tell it from an ordinary VUS.
    computed_conflicting: bool = False
    # Whether any annotator recognised the locus, and whether the reader
    # asked by protein change -- the two things that decide what an empty
    # ClinVar result means.
    annotated: bool = False
    asked_by_protein: bool = False


def _assertion_side(a) -> str:
    """Which side of the clinical line one record sits on."""
    if a.significance.is_pathogenic_side:
        return "pathogenic"
    if a.significance.is_benign_side:
        return "benign"
    return "uncertain"


def _primary_label(a) -> str:
    """The classification to print for a record.

    ClinVar joins additional classifications to the primary one with a
    semicolon, so a raw label can read "Uncertain significance; association;
    other". The primary clause is what this row was sorted on, and it is what
    the reader needs beside the computed call; the whole string is still on
    the record for anyone who opens it.
    """
    raw = (a.significance_raw or "").strip()
    if raw:
        return raw.split(";")[0].strip() or raw
    return a.significance.value


def _gather_facts(report: Report, today: date) -> _Facts:
    f = _Facts()
    f.unaccounted_high = sorted(
        {
            c.kind
            for c in (report.conflicts or [])
            if c.severity == "high" and c.kind not in TIER_ACCOUNTS_FOR
        }
    )
    call = report.classification.call if report.classification else None
    f.computed_side = _call_side(call)
    f.computed_conflicting = call == CONFLICTING
    f.asked_by_protein = bool(report.query.gene and report.query.protein_change)

    bundle = report.bundle
    f.annotated = bool(bundle and (bundle.clinvar or bundle.consequence))

    cv = bundle.clinvar if bundle else None
    if cv is None:
        return f

    actionable = cv.actionable_assertions
    f.records = len(actionable)
    if not actionable:
        return f

    # Only records that make a claim on the pathogenic/benign axis count
    # towards a split. "Drug response", "risk factor" and "pathogenic, low
    # penetrance" are answers to different questions, and counting them made
    # MTHFR A222V read as "8 records disagree: benign vs uncertain" when the
    # records agreed about pathogenicity and some of them were about warfarin
    # dosing. Reporting that as disagreement is the error this tool exists to
    # find, committed by this tool.
    on_axis = [a for a in actionable if a.significance in AXIS_LABELS]
    for a in on_axis:
        f.reported_sides.add(_assertion_side(a))
    f.axis_records = len(on_axis)
    f.off_axis_records = len(actionable) - len(on_axis)
    f.opposed_records = "pathogenic" in f.reported_sides and "benign" in f.reported_sides

    by_condition: dict[str, set[str]] = {}
    for a in on_axis:
        condition = (a.condition or "").strip().lower() or "(unstated)"
        by_condition.setdefault(condition, set()).add(_assertion_side(a))
    f.axis_conditions = len(by_condition)
    f.split_within_condition = any(len(sides) > 1 for sides in by_condition.values())

    # Everything below asks "how strong is the pathogenicity claim, and when
    # was it last looked at", so it is measured over the records that make
    # that claim. Measuring it over all of them let an unrelated record answer
    # for one: a three-star expert-panel record about methotrexate response
    # donated its stars to a one-star pathogenic assertion, moving the row out
    # of "thin evidence" entirely, and a drug-response record re-evaluated
    # this year hid a pathogenicity claim last examined in 2005.
    #
    # Where nothing is on the axis there is no pathogenicity claim to rate, so
    # the off-axis records answer for themselves rather than for a claim that
    # was never made.
    rated = on_axis or actionable

    known = [a for a in rated if a.stars != UNKNOWN_REVIEW_STARS]
    if known:
        f.best_stars = max(a.stars for a in known)
        best = max(known, key=lambda a: a.stars)
        f.reported_label = _primary_label(best)
        f.no_criteria = sum(1 for a in known if a.stars == 0)
        # The denominator is the set these were counted from, not every record
        # in the payload: "1 of 3 records state no criteria" was being printed
        # where one of the three had an unrecognised status and another was
        # not a pathogenicity claim at all.
        f.rated_records = len(known)
    elif rated:
        # Every record on the axis carries a status this build cannot read.
        # best_stars stays None, which the rule reads as "cannot assess".
        f.reported_label = _primary_label(rated[0])

    dated = [a for a in rated if a.last_evaluated]
    if dated:
        newest = max(a.last_evaluated for a in dated)
        f.stale_years = (today - newest).days / 365.25

    # What the off-axis records actually say, so the row can name them rather
    # than assert a category. The hard-coded "(drug response, risk factor or
    # low penetrance)" was printed over two ClinVar records that both read
    # "Conflicting interpretations of pathogenicity".
    f.off_axis_labels = sorted(
        {
            (a.significance_raw or a.significance.value).strip().lower()
            for a in actionable
            if a.significance not in AXIS_LABELS
        }
    )

    return f


def _reported_text(f: _Facts) -> str:
    """What to call the reported classification in a sentence.

    Falls back to the sides when every record carries a review status this
    build cannot read, because `reported_label` is only taken from records
    whose status it could rate.
    """
    if f.reported_label:
        return f.reported_label
    if f.reported_sides:
        return " and ".join(sorted(f.reported_sides))
    return "no classification"


def _off_axis_note(f: _Facts) -> list[str]:
    """Name the records that were left out of the comparison.

    The wording used to hard-code "(drug response, risk factor or low
    penetrance)", which was printed over two MTHFR records that both read
    "Conflicting interpretations of pathogenicity" -- telling the reader those
    answer a different question when they say the opposite.
    """
    if not f.off_axis_records:
        return []
    labels = ", ".join(f.off_axis_labels[:3])
    if len(f.off_axis_labels) > 3:
        labels += ", ..."
    return [
        f"{f.off_axis_records} further record(s) are not a pathogenic/benign "
        f"claim and are not counted here: {labels}"
    ]


def _unaccounted_note(f: _Facts) -> list[str]:
    """Name a high-severity finding the tier facts did not already cover.

    Without this the rule read only ClinVar's labels, so a variant where
    snpEff and Ensembl VEP disagree about the molecular consequence -- which
    means the worksheet may be classifying the wrong change -- could sit in
    "No disagreement found" with nothing to say about it.
    """
    if not f.unaccounted_high:
        return []
    named = [CONFLICT_LABELS.get(kind, kind.replace("_", " ")) for kind in f.unaccounted_high]
    return [f"Also flagged: {'; '.join(named)}"]


def _triage_facts(f: _Facts, resolved: bool, t: Thresholds) -> tuple[Tier, list[str]]:
    """The rule itself, over nothing but the facts above."""
    if not resolved:
        return UNRESOLVED, ["Input did not resolve to a variant"]

    reasons: list[str] = []

    # --- Tier 1: opposite sides of the clinical line -----------------------
    if f.opposed_records:
        reasons.append(
            f"ClinVar records disagree pathogenic vs benign "
            f"({f.axis_records} records)"
        )
    crosses = f.computed_side in ("pathogenic", "benign") and (
        ("benign" in f.reported_sides and f.computed_side == "pathogenic")
        or ("pathogenic" in f.reported_sides and f.computed_side == "benign")
    )
    if crosses:
        reasons.append(f"Computed {f.computed_side}, reported {_reported_text(f)}")
    if f.computed_conflicting:
        # The engine satisfied pathogenic *and* benign criteria on the same
        # variant. That is the strongest internal signal it has, and it was
        # invisible here: combine.py returns it as a flavour of "Uncertain
        # significance", so _call_side bucketed it with every ordinary VUS and
        # the row could land in "No disagreement found" carrying two
        # high-severity conflicts.
        reasons.append(
            "Pathogenic and benign criteria are both met on this variant"
        )
    if reasons:
        return CONTRADICTED, reasons

    # --- Tier 2: the records disagree with each other ---------------------
    if len(f.reported_sides) > 1:
        reasons.append(
            f"{f.axis_records} records disagree: " + " vs ".join(sorted(f.reported_sides))
        )
        # A ClinVar record is an RCV: one aggregate per variant-and-condition,
        # not one submitter's opinion. Two records differing may be two
        # readings of the same question, or the same variant assessed against
        # two different diseases -- for BRAF V600E, fifteen somatic tumour
        # conditions beside one germline RASopathy. The conflict detectors
        # already draw that distinction and the tier rule was ignoring it,
        # reporting "submitters split" over records that were not in
        # disagreement at all.
        if f.split_within_condition:
            reasons.append(
                "At least one condition has records on both sides, so this is "
                "disagreement about the same question"
            )
        elif f.axis_conditions > 1:
            reasons.append(
                f"The records cover {f.axis_conditions} distinct conditions, so "
                "part of the spread may be different clinical questions rather "
                "than disagreement"
            )
        reasons.extend(_off_axis_note(f))
        return SPLIT, reasons

    # --- Tier 3: agreement resting on little ------------------------------
    if f.records == 0:
        if not f.annotated:
            # Coordinates parse into a chrom and a position whether or not a
            # variant exists there, so `chr9:999999999A>T` used to arrive here
            # as "Thin evidence" and sort above a fully curated variant. If no
            # annotator recognised the locus either, nothing resolved.
            return UNRESOLVED, [
                "Parsed as coordinates, but no annotation source recognised "
                "this position. Check the assembly and the position."
            ]
        reasons = ["No ClinVar record for this variant"]
        if f.asked_by_protein:
            # Only worth saying where it could be the cause. It was printed
            # for coordinate queries too, which carry no protein numbering to
            # be wrong about.
            reasons.append(
                "Resolved from the canonical transcript; if this variant is "
                "known clinically under a different protein number, search "
                "the c. notation or the rsID instead"
            )
        return THIN, reasons
    if f.best_stars == 0:
        reasons.append("No record states assertion criteria (0 of 4 stars)")
    elif f.best_stars == 1:
        # ClinVar's one star means either a single submitter or conflicting
        # classifications. Both are thin, so the star rating is enough on its
        # own; `number_submitters` comes off an RCV aggregate and counts
        # submissions to a condition, which is not the same question.
        reasons.append("1 of 4 review stars (single submitter, or conflicting)")
    if f.no_criteria and f.best_stars and f.best_stars < EXPERT_REVIEW_STARS:
        reasons.append(
            f"{f.no_criteria} of {f.rated_records} records state no criteria"
        )
    if f.stale_years is not None and f.stale_years >= t.stale_assertion_years:
        reasons.append(f"Last evaluated {f.stale_years:.0f} years ago")
    if f.best_stars is None:
        # Every record carries a review status this build does not recognise,
        # so there is no evidence either way about how well reviewed it is.
        # Reading that as four stars would be generous and as zero would be
        # slander; it belongs with the records nobody has vouched for.
        return THIN, [
            f"Review status not recognised on any of {f.records} record(s); "
            "cannot tell how well reviewed this is"
        ]
    reasons.extend(_unaccounted_note(f))
    if reasons:
        return THIN, reasons

    # --- Tier 4 -----------------------------------------------------------
    return CONCORDANT, [f"Agrees with {f.records} record(s), {f.best_stars} of 4 stars"]


def _is_resolved(report: Report) -> bool:
    if report.query.chrom and report.query.pos:
        return True
    return bool(report.bundle and (report.bundle.clinvar or report.bundle.consequence))


def triage(report: Report, t: Thresholds = DEFAULT, today: date | None = None) -> tuple[Tier, list[str]]:
    """Assign a tier and the reasons for it.

    Pure function over a finished report, so the rule is testable without a
    network and the UI, the CLI and the CSV export cannot drift apart.
    """
    today = today or date.today()
    return _triage_facts(_gather_facts(report, today), _is_resolved(report), t)


def evidence_not_public(f: _Facts) -> bool:
    """Reported confidently, computed uncertain: the lab has evidence we cannot see.

    Requires the records to agree. `reported_sides <= {"pathogenic","benign"}`
    was also satisfied by a set holding *both*, so a contradicted row printed
    "the submitter's evidence is not public -- expected, not a flag" directly
    underneath "these records disagree pathogenic vs benign". Two submitters
    on opposite sides is not one submitter holding private evidence.

    Also excludes the case where the engine itself met both sides: that is a
    conflict of its own, not deference to a laboratory.
    """
    if f.computed_conflicting or f.opposed_records:
        return False
    return f.computed_side == "uncertain" and len(f.reported_sides) == 1 and (
        f.reported_sides <= {"pathogenic", "benign"}
    )


def row_from_report(
    query: str, report: Report, t: Thresholds = DEFAULT, today: date | None = None
) -> PanelRow:
    f = _gather_facts(report, today or date.today())
    tier, reasons = _triage_facts(f, _is_resolved(report), t)
    conflicts = report.conflicts or []
    resolved = tier is not UNRESOLVED
    return PanelRow(
        query=query,
        resolved=resolved,
        display_name=report.display_name,
        gene=report.gene_symbol,
        # Suppressed when nothing resolved. The engine still produces
        # "Uncertain significance" over an empty bundle, and printing that next
        # to a failed lookup reads as a result rather than as the absence of
        # one.
        computed_call=(report.classification.call if report.classification else None)
        if resolved
        else None,
        reported_call=f.reported_label,
        best_stars=f.best_stars,
        records=f.records,
        conflicts=len(conflicts),
        high_severity=sum(1 for c in conflicts if c.severity == "high"),
        conflict_headlines=[c.headline for c in conflicts],
        met_codes=report.classification.met_codes if report.classification else [],
        tier=tier.key,
        tier_rank=tier.rank,
        evidence_not_public=evidence_not_public(f),
        reasons=reasons,
    )


def sort_rows(rows: list[PanelRow]) -> list[PanelRow]:
    """Worklist order: tier first, then how much there is to read.

    Within a tier, more high-severity conflicts come first, then more
    conflicts, then *fewer* review stars -- a 1-star record with a problem is a
    worse place to be than a 3-star one with the same problem, because nobody
    has checked it.
    """
    return sorted(
        rows,
        key=lambda r: (
            r.tier_rank,
            -r.high_severity,
            -r.conflicts,
            r.best_stars if r.best_stars is not None else 9,
            # Case-folded, because the page breaks this tie with
            # localeCompare, which ignores case, while a raw Python string
            # comparison puts every capital before every lower-case letter.
            # The two disagreed whenever a query was typed in lower case, so
            # the downloaded CSV was not in the order on screen -- in a file
            # presented as the audit trail for what the reader just looked at.
            (r.display_name or r.query).casefold(),
            r.query,
        ),
    )


def tier_counts(rows: list[PanelRow]) -> dict[str, int]:
    counts = {t.key: 0 for t in TIERS}
    for r in rows:
        counts[r.tier] = counts.get(r.tier, 0) + 1
    return counts


# --------------------------------------------------------------------------
# Running a panel
# --------------------------------------------------------------------------


# A query longer than this is not a variant identifier. The longest real one
# here is an HGVS string with a transcript accession and a large indel, which
# stays well inside it.
MAX_QUERY_LENGTH = 200


def parse_panel(text: str, limit: int = MAX_PANEL_SIZE) -> list[str]:
    """Split pasted text into queries.

    Accepts one per line, or comma-separated, and tolerates the shapes a
    spreadsheet column arrives in: surrounding quotes, a trailing comma, a
    leading '#' comment, blank lines. Duplicates are dropped, keeping the
    first, because a panel exported twice is a common accident and running a
    variant twice costs the public APIs a request for nothing.

    Over-long entries are dropped rather than sent on: the string reaches six
    public APIs, and nothing a variant identifier can legitimately be is 200
    characters.
    """
    return parse_panel_counted(text, limit)[0]


def parse_panel_counted(text: str, limit: int = MAX_PANEL_SIZE) -> tuple[list[str], int]:
    """As `parse_panel`, and also how many entries were dropped at the limit.

    Truncating silently meant a 200-line file produced a 50-row worklist with
    nothing saying so, in a tool whose whole purpose is to be read as a
    complete account of what is known.
    """
    out: list[str] = []
    seen: set[str] = set()
    dropped = 0
    for line in text.replace(",", "\n").splitlines():
        q = line.strip().strip('"').strip("'").strip()
        if not q or q.startswith("#"):
            continue
        if len(q) > MAX_QUERY_LENGTH:
            dropped += 1
            continue
        key = q.upper()
        if key in seen:
            continue
        seen.add(key)
        if len(out) >= limit:
            dropped += 1
            continue
        out.append(q)
    return out, dropped


async def interpret_panel(
    queries: list[str],
    assembly: str = "hg38",
    cache: ResponseCache | None = None,
    client: httpx.AsyncClient | None = None,
    thresholds: Thresholds = DEFAULT,
    concurrency: int = MAX_CONCURRENCY,
    per_variant_timeout: float = 70.0,
    today: date | None = None,
) -> PanelResult:
    """Interpret every query, then order the results.

    One variant failing never fails the panel. A panel that returns nothing
    because the fortieth lookup timed out would be useless in exactly the
    situation it exists for, so a failure becomes an Unresolved row carrying
    its own reason.
    """
    started = datetime.now(timezone.utc)
    sem = asyncio.Semaphore(max(1, concurrency))

    def _elapsed_ms(since: datetime) -> int:
        return int((datetime.now(timezone.utc) - since).total_seconds() * 1000)

    async def one(q: str) -> PanelRow:
        async with sem:
            # Timed from here, not from the call. Measuring across the wait for
            # a semaphore slot reported a variant that queued for three minutes
            # as a three-minute lookup, which says nothing about the sources.
            t0 = datetime.now(timezone.utc)
            try:
                report = await asyncio.wait_for(
                    interpret(q, assembly=assembly, cache=cache, client=client, thresholds=thresholds),
                    timeout=per_variant_timeout,
                )
            except asyncio.TimeoutError:
                return PanelRow(
                    query=q,
                    error="Sources did not respond in time",
                    reasons=["Lookup timed out"],
                    elapsed_ms=_elapsed_ms(t0),
                )
            except Exception as exc:  # noqa: BLE001 - one bad row must not sink the panel
                return PanelRow(
                    query=q,
                    error=f"{type(exc).__name__}: {exc}"[:200],
                    reasons=["Lookup failed"],
                    elapsed_ms=_elapsed_ms(t0),
                )

            # Inside the try as well. This is pure computation over a finished
            # report, but it was the one call in the coroutine that could still
            # raise -- and because gather() defaults to propagating, a single
            # unexpected exception here discarded every other row in the panel.
            try:
                row = row_from_report(q, report, t=thresholds, today=today)
            except Exception as exc:  # noqa: BLE001
                return PanelRow(
                    query=q,
                    display_name=report.display_name,
                    gene=report.gene_symbol,
                    error=f"Could not triage: {type(exc).__name__}: {exc}"[:200],
                    reasons=["Interpreted, but the triage rule failed on it"],
                    elapsed_ms=_elapsed_ms(t0),
                )
            row.elapsed_ms = _elapsed_ms(t0)
            return row

    # return_exceptions so that a failure the per-row handling did not
    # anticipate -- a cancellation, an error raised between the handlers --
    # costs that row rather than the whole panel.
    settled = await asyncio.gather(*(one(q) for q in queries), return_exceptions=True)
    rows: list[PanelRow] = []
    for q, outcome in zip(queries, settled):
        if isinstance(outcome, PanelRow):
            rows.append(outcome)
        else:
            rows.append(
                PanelRow(
                    query=q,
                    error=f"{type(outcome).__name__}: {outcome}"[:200],
                    reasons=["Lookup failed"],
                )
            )
    ordered = sort_rows(rows)
    return PanelResult(
        rows=ordered,
        tier_counts=tier_counts(ordered),
        evidence_not_public=sum(1 for r in ordered if r.evidence_not_public),
        generated_at=started.isoformat(timespec="seconds"),
        elapsed_ms=int((datetime.now(timezone.utc) - started).total_seconds() * 1000),
    )


# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------

CSV_COLUMNS = (
    "tier",
    "query",
    "variant",
    "gene",
    "computed_call",
    "reported_call",
    "review_stars",
    "clinvar_records",
    "conflicts",
    "high_severity_conflicts",
    "acmg_codes",
    "reasons",
    "evidence_not_public",
    "error",
)


# Spreadsheets treat a cell beginning with any of these as a formula, not as
# text, and open it without asking. The cells below carry a reader's own query
# string and text derived from public sources, so neither is trusted input.
_FORMULA_LEAD = ("=", "+", "-", "@", "\t", "\r")


def _csv_safe(value: object) -> object:
    """Keep a cell a cell.

    Excel, LibreOffice and Google Sheets all read a leading `=`, `+`, `-` or
    `@` as the start of a formula, so a variant name or an error message could
    execute on open -- in a file this tool hands over as the record of what a
    reviewer looked at. Prefixing an apostrophe is the documented way to force
    text; it is stripped by the spreadsheet on display.

    Numbers are returned untouched: they are generated here, not received, and
    quoting them would make every count a string in the sheet.
    """
    if not isinstance(value, str) or not value.startswith(_FORMULA_LEAD):
        return value
    return "'" + value


def rows_to_csv(rows: list[PanelRow]) -> str:
    """A worklist a reviewer can open in a spreadsheet and annotate.

    Written by hand rather than through pandas: the dependency would be there
    for one function, and the quoting rules for this data are not subtle.
    """
    import csv
    import io

    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(CSV_COLUMNS)
    for r in rows:
        tier = TIERS_BY_KEY.get(r.tier)
        w.writerow(
            _csv_safe(cell)
            for cell in [
                tier.name if tier else r.tier,
                r.query,
                r.display_name or "",
                r.gene or "",
                r.computed_call or "",
                r.reported_call or "",
                "" if r.best_stars is None else r.best_stars,
                r.records,
                r.conflicts,
                r.high_severity,
                " ".join(r.met_codes),
                "; ".join(r.reasons),
                "yes" if r.evidence_not_public else "",
                r.error or "",
            ]
        )
    return buf.getvalue()
