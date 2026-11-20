"""`interpret_panel` and the panel endpoints.

The module's headline claim is that one variant failing never fails the panel.
Nothing tested it: the whole runtime half -- concurrency, the semaphore, the
per-variant timeout, failure isolation, the aggregate counts -- had no tests at
all, and neither did any endpoint, which is how an unbounded export and a
silent truncation at fifty both survived.

Everything here stubs `interpret`, so it runs with the network severed.
"""

from __future__ import annotations

import asyncio
from datetime import date

import pytest
from fastapi.testclient import TestClient

from dissent.acmg.combine import PATHOGENIC, UNCERTAIN
from dissent.api.app import MAX_ROWS_PER_EXPORT, app
from dissent.models import (
    Classification,
    ClinVarAssertion,
    ClinVarEvidence,
    Consequence,
    EvidenceBundle,
    Report,
    Significance,
    VariantQuery,
)
from dissent.panel import (
    CONCORDANT,
    MAX_PANEL_SIZE,
    MAX_QUERY_LENGTH,
    UNRESOLVED,
    interpret_panel,
    parse_panel_counted,
)

TODAY = date(2026, 11, 20)


def a_report(query: str, call: str = PATHOGENIC) -> Report:
    q = VariantQuery(raw=query, kind="gene_protein", gene="TEST", protein_change="p.Ala1Val")
    bundle = EvidenceBundle(
        query=q,
        clinvar=ClinVarEvidence(
            variant_id=1,
            gene="TEST",
            assertions=[
                ClinVarAssertion(
                    accession="RCV1",
                    significance_raw="Pathogenic",
                    significance=Significance.PATHOGENIC,
                    review_status="criteria provided, multiple submitters, no conflicts",
                    stars=2,
                    last_evaluated=date(2025, 6, 1),
                )
            ],
        ),
        consequence=Consequence(most_severe="missense_variant"),
    )
    return Report(
        query=q,
        gene_symbol="TEST",
        display_name=query,
        classification=Classification(call=call, reasoning="t", met_codes=["PM2"]),
        bundle=bundle,
    )


@pytest.fixture
def stub_interpret(monkeypatch):
    """Replace the network half with a dictionary of behaviours.

    Keys map a query to "ok", "boom", "hang", or "untriageable".
    """
    behaviour: dict[str, str] = {}

    async def fake(raw, **_kwargs):
        what = behaviour.get(raw, "ok")
        if what == "boom":
            raise RuntimeError("source exploded")
        if what == "hang":
            await asyncio.sleep(30)
        if what == "untriageable":
            # A report the triage rule chokes on: no query object at all.
            return Report.model_construct(query=None, classification=None, bundle=None)
        return a_report(raw)

    monkeypatch.setattr("dissent.panel.interpret", fake)
    return behaviour


# --------------------------------------------------------------------------
# Failure isolation -- the module's headline claim
# --------------------------------------------------------------------------


async def test_a_raising_variant_does_not_sink_the_panel(stub_interpret):
    stub_interpret["BOOM"] = "boom"
    result = await interpret_panel(["ok1", "BOOM", "ok2"], today=TODAY)

    assert len(result.rows) == 3
    bad = next(r for r in result.rows if r.query == "BOOM")
    assert bad.tier == UNRESOLVED.key
    assert "source exploded" in bad.error
    assert sum(1 for r in result.rows if r.tier == CONCORDANT.key) == 2


async def test_a_hanging_variant_times_out_without_sinking_the_panel(stub_interpret):
    stub_interpret["SLOW"] = "hang"
    result = await interpret_panel(["ok1", "SLOW", "ok2"], per_variant_timeout=0.05, today=TODAY)

    slow = next(r for r in result.rows if r.query == "SLOW")
    assert slow.tier == UNRESOLVED.key
    assert "did not respond in time" in slow.error
    assert sum(1 for r in result.rows if r.tier == CONCORDANT.key) == 2


async def test_a_failure_inside_the_triage_rule_costs_one_row(stub_interpret):
    """`row_from_report` ran outside the try, and gather() propagates by
    default, so one unexpected exception there discarded every other row."""
    stub_interpret["WEIRD"] = "untriageable"
    result = await interpret_panel(["ok1", "WEIRD", "ok2"], today=TODAY)

    assert len(result.rows) == 3
    weird = next(r for r in result.rows if r.query == "WEIRD")
    assert weird.error
    assert sum(1 for r in result.rows if r.tier == CONCORDANT.key) == 2


async def test_every_query_gets_exactly_one_row(stub_interpret):
    stub_interpret["BOOM"] = "boom"
    queries = ["a", "BOOM", "b", "c"]
    result = await interpret_panel(queries, today=TODAY)
    assert sorted(r.query for r in result.rows) == sorted(queries)


# --------------------------------------------------------------------------
# Concurrency
# --------------------------------------------------------------------------


async def test_concurrency_is_bounded(monkeypatch):
    """Each variant fans out to six unpaid public APIs. A panel that ignored
    the limit would be rate-limited into being slower, not faster."""
    live = 0
    peak = 0

    async def fake(raw, **_kwargs):
        nonlocal live, peak
        live += 1
        peak = max(peak, live)
        try:
            await asyncio.sleep(0.01)
            return a_report(raw)
        finally:
            live -= 1

    monkeypatch.setattr("dissent.panel.interpret", fake)
    await interpret_panel([f"q{i}" for i in range(20)], concurrency=3, today=TODAY)
    assert peak <= 3


async def test_a_failing_variant_releases_its_slot(monkeypatch):
    """A leaked semaphore slot would stall the panel after a handful of
    failures rather than failing loudly."""
    live = 0
    peak = 0

    async def fake(raw, **_kwargs):
        nonlocal live, peak
        live += 1
        peak = max(peak, live)
        try:
            await asyncio.sleep(0.005)
            raise RuntimeError("always fails")
        finally:
            live -= 1

    monkeypatch.setattr("dissent.panel.interpret", fake)
    result = await asyncio.wait_for(
        interpret_panel([f"q{i}" for i in range(12)], concurrency=2, today=TODAY), timeout=10
    )
    assert len(result.rows) == 12
    assert peak <= 2


async def test_the_queue_wait_is_not_counted_against_the_lookup(monkeypatch):
    """Timing started before the semaphore, so a variant that queued for three
    minutes reported a three-minute lookup -- which says nothing about the
    sources."""

    async def fake(raw, **_kwargs):
        await asyncio.sleep(0.01)
        return a_report(raw)

    monkeypatch.setattr("dissent.panel.interpret", fake)
    result = await interpret_panel([f"q{i}" for i in range(8)], concurrency=1, today=TODAY)
    # Serialised, the last row waits ~70ms for its slot. Its own lookup is 10ms.
    assert all(r.elapsed_ms is not None and r.elapsed_ms < 60 for r in result.rows)


async def test_failed_rows_still_report_their_timing(stub_interpret):
    """Both failure paths returned before elapsed_ms was set, so it was absent
    on exactly the rows whose timing a reader would want."""
    stub_interpret["BOOM"] = "boom"
    stub_interpret["SLOW"] = "hang"
    result = await interpret_panel(["BOOM", "SLOW"], per_variant_timeout=0.05, today=TODAY)
    assert all(r.elapsed_ms is not None for r in result.rows)


# --------------------------------------------------------------------------
# Aggregates
# --------------------------------------------------------------------------


async def test_the_result_counts_agree_with_the_rows(stub_interpret):
    stub_interpret["BOOM"] = "boom"
    result = await interpret_panel(["a", "BOOM", "b"], today=TODAY)

    assert sum(result.tier_counts.values()) == len(result.rows)
    assert result.tier_counts[UNRESOLVED.key] == 1
    assert result.evidence_not_public == sum(1 for r in result.rows if r.evidence_not_public)
    assert result.generated_at
    assert result.elapsed_ms is not None


async def test_an_empty_panel_is_an_empty_result(stub_interpret):
    result = await interpret_panel([], today=TODAY)
    assert result.rows == []
    assert sum(result.tier_counts.values()) == 0


async def test_rows_come_back_in_worklist_order_not_request_order(monkeypatch):
    async def fake(raw, **_kwargs):
        return a_report(raw, call=UNCERTAIN if raw == "second" else PATHOGENIC)

    monkeypatch.setattr("dissent.panel.interpret", fake)
    result = await interpret_panel(["first", "second"], today=TODAY)
    assert [r.tier_rank for r in result.rows] == sorted(r.tier_rank for r in result.rows)


# --------------------------------------------------------------------------
# Parsing limits
# --------------------------------------------------------------------------


@pytest.mark.parametrize("n,kept,dropped", [(49, 49, 0), (50, 50, 0), (51, 50, 1), (200, 50, 150)])
def test_the_panel_cap_is_reported_not_silent(n, kept, dropped):
    queries, n_dropped = parse_panel_counted("\n".join(f"rs{i}" for i in range(n)))
    assert len(queries) == kept
    assert n_dropped == dropped


def test_an_over_long_entry_is_dropped_and_counted():
    """The string reaches six public APIs. Nothing a variant identifier can
    legitimately be is 200 characters."""
    queries, dropped = parse_panel_counted("A" * (MAX_QUERY_LENGTH + 1) + "\nrs334")
    assert queries == ["rs334"]
    assert dropped == 1


def test_an_entry_at_the_length_limit_is_kept():
    queries, dropped = parse_panel_counted("A" * MAX_QUERY_LENGTH)
    assert len(queries) == 1 and dropped == 0


# --------------------------------------------------------------------------
# The endpoints
# --------------------------------------------------------------------------


@pytest.fixture
def client(stub_interpret):
    with TestClient(app) as c:
        yield c


def test_panel_endpoint_returns_rows(client):
    r = client.post("/api/panel", json={"queries": ["BRAF V600E", "TP53 R175H"]})
    assert r.status_code == 200
    body = r.json()
    assert len(body["rows"]) == 2
    assert body["dropped"] == 0


def test_panel_endpoint_reports_what_it_dropped(client):
    """A 200-line file produced a 50-row worklist with nothing saying so."""
    r = client.post("/api/panel", json={"queries": [f"rs{i}" for i in range(80)]})
    assert r.status_code == 200
    body = r.json()
    assert len(body["rows"]) == MAX_PANEL_SIZE
    assert body["dropped"] == 30


@pytest.mark.parametrize(
    "payload",
    [
        {"queries": []},
        {"queries": 5},
        {"queries": [{"a": 1}]},
        {"queries": [None]},
        {"queries": ["rs334"], "assembly": "hg99"},
        {},
    ],
)
def test_panel_endpoint_rejects_bad_input(client, payload):
    assert client.post("/api/panel", json=payload).status_code == 422


def test_panel_endpoint_does_not_stringify_non_strings(client):
    """`str(q)` turned {"a": 1} into the query "{'a': 1}" and sent it on."""
    r = client.post("/api/panel", json={"queries": [{"a": 1}]})
    assert r.status_code == 422
    assert "string" in r.json()["detail"].lower()


def test_panel_endpoint_refuses_an_absurd_body(client):
    r = client.post("/api/panel", json={"queries": ["A" * 200] * 500})
    assert r.status_code == 413


def test_csv_endpoint_round_trips_a_panel(client):
    panel = client.post("/api/panel", json={"queries": ["BRAF V600E"]}).json()
    r = client.post("/api/panel/csv", json={"rows": panel["rows"]})
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/csv")
    assert "BRAF V600E" in r.text


def test_csv_endpoint_caps_the_number_of_rows(client):
    """Validated each row's shape, never how many: a request could ask for an
    eight-megabyte file to be built in memory."""
    rows = [{"query": f"rs{i}"} for i in range(MAX_ROWS_PER_EXPORT + 1)]
    assert client.post("/api/panel/csv", json={"rows": rows}).status_code == 413


def test_csv_endpoint_neutralises_a_formula(client):
    """A cell beginning `=`, `+`, `-` or `@` is a formula in Excel and Sheets,
    in a file this tool hands over as the record of what a reviewer saw."""
    rows = [{"query": "=HYPERLINK(\"http://x\")", "tier": "thin", "tier_rank": 3}]
    r = client.post("/api/panel/csv", json={"rows": rows})
    assert r.status_code == 200
    assert "'=HYPERLINK" in r.text
    assert "\n=HYPERLINK" not in r.text


@pytest.mark.parametrize("payload", [{"rows": []}, {"rows": [{"nope": 1}]}, {}])
def test_csv_endpoint_rejects_bad_input(client, payload):
    assert client.post("/api/panel/csv", json=payload).status_code == 422


def test_tiers_endpoint_describes_every_tier(client):
    body = client.get("/api/panel/tiers").json()
    assert {t["key"] for t in body["tiers"]} == {t.key for t in __import__(
        "dissent.panel", fromlist=["TIERS"]
    ).TIERS}
    assert body["max_per_request"] == MAX_PANEL_SIZE
    assert body["evidence_not_public"]


def test_examples_endpoint_offers_runnable_panels(client):
    panels = client.get("/api/panel/examples").json()
    assert panels
    for p in panels:
        assert p["name"] and p["note"]
        assert 0 < len(p["queries"]) <= MAX_PANEL_SIZE
        assert all(isinstance(q, str) and q.strip() for q in p["queries"])
