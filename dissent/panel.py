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

from dissent.acmg.combine import BENIGN, LIKELY_BENIGN, LIKELY_PATHOGENIC, PATHOGENIC
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

# A cap per request. The web UI sends a long panel in chunks so the reader sees
# rows arrive; this is the ceiling on one chunk.
MAX_PANEL_SIZE = 50


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
    "Submitters split",
    "The submitted records disagree with each other, but not across the "
    "clinical line -- pathogenic against uncertain, typically. The aggregate "
    "label a report quotes picks one of them.",
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
_AXIS = frozenset(
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
    stale_years: float | None = None
    opposed_records: bool = False


def _gather_facts(report: Report, t: Thresholds, today: date) -> _Facts:
    f = _Facts()
    f.computed_side = _call_side(report.classification.call if report.classification else None)

    cv = report.bundle.clinvar if report.bundle else None
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
    on_axis = [a for a in actionable if a.significance in _AXIS]
    for a in on_axis:
        if a.significance.is_pathogenic_side:
            f.reported_sides.add("pathogenic")
        elif a.significance.is_benign_side:
            f.reported_sides.add("benign")
        else:
            f.reported_sides.add("uncertain")
    f.axis_records = len(on_axis)
    f.off_axis_records = len(actionable) - len(on_axis)
    f.opposed_records = "pathogenic" in f.reported_sides and "benign" in f.reported_sides

    # The label shown to the reader: the assertion with the most review stars,
    # preferring one that is on the axis, because "drug response" as the
    # headline next to a computed benign call reads as a contradiction and is
    # not one. Unknown statuses are excluded so an unrecognised string cannot
    # win the comparison at -1.
    known = [a for a in actionable if a.stars != UNKNOWN_REVIEW_STARS]
    if known:
        f.best_stars = max(a.stars for a in known)
        preferred = [a for a in known if a.significance in _AXIS] or known
        best = max(preferred, key=lambda a: a.stars)
        f.reported_label = best.significance_raw or best.significance.value
        f.no_criteria = sum(1 for a in known if a.stars == 0)

    dated = [a for a in actionable if a.last_evaluated]
    if dated:
        newest = max(a.last_evaluated for a in dated)
        f.stale_years = (today - newest).days / 365.25

    return f


def _reported_text(f: _Facts) -> str:
    return f.reported_label or "/".join(sorted(f.reported_sides))


def _triage_facts(f: _Facts, resolved: bool, t: Thresholds) -> tuple[Tier, list[str]]:
    """The rule itself, over nothing but the facts above."""
    if not resolved:
        return UNRESOLVED, ["Input did not resolve to a variant"]

    reasons: list[str] = []

    # --- Tier 1: opposite sides of the clinical line -----------------------
    if f.opposed_records:
        reasons.append(
            f"Submitted records disagree pathogenic vs benign ({f.axis_records} records)"
        )
    crosses = f.computed_side in ("pathogenic", "benign") and (
        ("benign" in f.reported_sides and f.computed_side == "pathogenic")
        or ("pathogenic" in f.reported_sides and f.computed_side == "benign")
    )
    if crosses:
        reasons.append(f"Computed {f.computed_side}, reported {_reported_text(f)}")
    if reasons:
        return CONTRADICTED, reasons

    # --- Tier 2: the records disagree with each other ---------------------
    if len(f.reported_sides) > 1:
        reasons.append(
            f"{f.axis_records} records disagree: " + " vs ".join(sorted(f.reported_sides))
        )
        if f.off_axis_records:
            reasons.append(
                f"{f.off_axis_records} further record(s) answer a different question "
                "(drug response, risk factor or low penetrance) and are not counted here"
            )
        return SPLIT, reasons

    # --- Tier 3: agreement resting on little ------------------------------
    if f.records == 0:
        return THIN, [
            "No ClinVar record. Resolved from the canonical transcript; if this "
            "variant is known clinically under a different protein number, "
            "search the c. notation instead"
        ]
    if f.best_stars == 0:
        reasons.append("No record states assertion criteria (0 of 4 stars)")
    elif f.best_stars == 1:
        # ClinVar's one star means either a single submitter or conflicting
        # classifications. Both are thin, so the star rating is enough on its
        # own; `number_submitters` comes off an RCV aggregate and counts
        # submissions to a condition, which is not the same question.
        reasons.append("1 of 4 review stars (single submitter, or conflicting)")
    if f.no_criteria and f.best_stars and f.best_stars < 3:
        reasons.append(f"{f.no_criteria} of {f.records} records state no criteria")
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
    return _triage_facts(_gather_facts(report, t, today), _is_resolved(report), t)


def evidence_not_public(f: _Facts) -> bool:
    """Reported confidently, computed uncertain: the lab has evidence we cannot see."""
    return f.computed_side == "uncertain" and bool(
        f.reported_sides and f.reported_sides <= {"pathogenic", "benign"}
    )


def row_from_report(
    query: str, report: Report, t: Thresholds = DEFAULT, today: date | None = None
) -> PanelRow:
    f = _gather_facts(report, t, today or date.today())
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
            r.display_name or r.query,
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


def parse_panel(text: str, limit: int = MAX_PANEL_SIZE) -> list[str]:
    """Split pasted text into queries.

    Accepts one per line, or comma-separated, and tolerates the shapes a
    spreadsheet column arrives in: surrounding quotes, a trailing comma, a
    leading '#' comment, blank lines. Duplicates are dropped, keeping the
    first, because a panel exported twice is a common accident and running a
    variant twice costs the public APIs a request for nothing.
    """
    out: list[str] = []
    seen: set[str] = set()
    for line in text.replace(",", "\n").splitlines():
        q = line.strip().strip('"').strip("'").strip()
        if not q or q.startswith("#"):
            continue
        key = q.upper()
        if key in seen:
            continue
        seen.add(key)
        out.append(q)
        if len(out) >= limit:
            break
    return out


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

    async def one(q: str) -> PanelRow:
        t0 = datetime.now(timezone.utc)
        async with sem:
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
                )
            except Exception as exc:  # noqa: BLE001 - one bad row must not sink the panel
                return PanelRow(
                    query=q,
                    error=f"{type(exc).__name__}: {exc}"[:200],
                    reasons=["Lookup failed"],
                )
        row = row_from_report(q, report, t=thresholds, today=today)
        row.elapsed_ms = int((datetime.now(timezone.utc) - t0).total_seconds() * 1000)
        return row

    rows = await asyncio.gather(*(one(q) for q in queries))
    ordered = sort_rows(list(rows))
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
            [
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
