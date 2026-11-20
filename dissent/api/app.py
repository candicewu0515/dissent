"""FastAPI application: JSON API plus the static single-page UI."""

from __future__ import annotations

import asyncio
from contextlib import asynccontextmanager
import hashlib
import json
from pathlib import Path

import httpx
from fastapi import Body, FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

from dissent.baseline.rates import load_baseline
from dissent.cache.store import ResponseCache
from dissent.models import Report
from dissent.acmg.thresholds import DEFAULT as DEFAULT_THRESHOLDS
from dissent.panel import (
    AXIS_LABELS,
    EVIDENCE_NOT_PUBLIC,
    EXPERT_REVIEW_STARS,
    MAX_PANEL_SIZE,
    MAX_QUERY_LENGTH,
    TIERS,
    PanelResult,
    PanelRow,
    interpret_panel,
    parse_panel_counted,
    rows_to_csv,
    sort_rows,
)
from dissent.report import interpret
from dissent.sources.base import USER_AGENT

WEB_DIR = Path(__file__).resolve().parent.parent / "web"

# A whole-chunk ceiling, so a request cannot hold a connection for as long
# as its per-variant timeouts multiplied together allow.
PANEL_TIMEOUT_SECONDS = 180.0

# A panel is assembled from chunks, so an export can legitimately be larger
# than one request's worth -- but not unbounded.
MAX_ROWS_PER_EXPORT = 2000


def _asset_version() -> str:
    """A short digest of the static assets, appended to their URLs.

    Without it a browser keeps serving the stylesheet it already has: a CSS
    change looked like it had not applied at all, which is an easy hour to
    lose and would be worse after a deploy than in development.
    """
    digest = hashlib.sha256()
    for path in sorted(WEB_DIR.glob("*")):
        if path.suffix in (".css", ".js"):
            digest.update(path.read_bytes())
    return digest.hexdigest()[:10]


ASSET_VERSION = _asset_version()


def _page(name: str) -> HTMLResponse:
    """Serve a page with its asset links versioned."""
    html = (WEB_DIR / name).read_text()
    html = html.replace('.css"', f'.css?v={ASSET_VERSION}"')
    html = html.replace('.js"', f'.js?v={ASSET_VERSION}"')
    return HTMLResponse(html)
DATA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "clinvar_baseline.json"

# Curated demo variants. Each was chosen because it exercises a different
# failure mode, and the set doubles as the walkthrough in the README.
EXAMPLES = [
    {
        "query": "MTHFR A222V",
        "label": "MTHFR A222V",
        "why": "49% allele frequency in one population, yet most predictors call it damaging. BA1 and PP3 are met at once.",
    },
    {
        "query": "HFE C282Y",
        "label": "HFE C282Y",
        "why": "Asserted pathogenic in ClinVar at ~7% allele frequency, above the BA1 stand-alone benign threshold.",
    },
    {
        "query": "BRAF V600E",
        "label": "BRAF V600E",
        "why": "Somatic and germline assertions mixed in one record; most carry no assertion criteria.",
    },
    {
        "query": "CFTR G542X",
        "label": "CFTR G542X",
        "why": "Nonsense variant in a recessive gene, where gnomAD constraint is uninformative by construction.",
    },
    {
        "query": "NF1 R1947*",
        "label": "NF1 R1947*",
        "why": "PVS1 is met, but the evidence falls into a gap in the published combining rules.",
    },
    {
        "query": "KCNH2 A561V",
        "label": "KCNH2 A561V",
        "why": "Reaches Likely pathogenic from public data alone — a long QT variant where the frequency, the gene's constraint and a calibrated predictor all line up.",
    },
    {
        "query": "BRCA1 C61G",
        "label": "BRCA1 C61G",
        "why": "Well curated and unanimous. Shown so the conflict panel can be seen staying quiet.",
    },
]


# Two panels a reader would recognise, each a real gene set rather than a
# selection of variants chosen to look good. Every variant below was checked
# against the live sources; the tiers they land in are whatever the data says,
# including the ones that do not resolve, which are left in on purpose.
PANEL_EXAMPLES = [
    {
        "name": "Cardiomyopathy and arrhythmia",
        "note": (
            "Fourteen variants from the genes on a standard inherited cardiac "
            "panel. Three sit on opposite sides of the clinical line from each "
            "other; most of the rest are reported confidently on one star."
        ),
        "queries": [
            "MYH7 R403Q",
            "MYBPC3 R502W",
            "TNNT2 R92Q",
            "TNNI3 R145W",
            "MYL2 R58Q",
            "KCNQ1 A341V",
            "KCNH2 A561V",
            "SCN5A E1784K",
            "RYR2 R176Q",
            "CASQ2 D307H",
            "LMNA R482W",
            "DSP R1267Q",
            "PKP2 R735X",
            "TTN R8604X",
        ],
    },
    {
        "name": "Carrier screening",
        "note": (
            "Recessive-disease variants of the kind a reproductive carrier "
            "screen reports. Includes HBB E6V, which does not resolve: the "
            "clinical name uses legacy protein numbering and HGVS calls the "
            "same variant p.Glu7Val. Left in because that failure is the point."
        ),
        "queries": [
            "CFTR G542X",
            "HFE C282Y",
            "GJB2 V37I",
            "G6PD V68M",
            "ATP7B H1069Q",
            "PAH R408W",
            "GALT Q188R",
            "BTD D444H",
            "MTHFR A222V",
            "HBB E6V",
        ],
    },
]


@asynccontextmanager
async def lifespan(app: FastAPI):
    # One shared connection pool and cache for the process. Opening a client
    # per request would discard connection reuse, which matters when six
    # sources are queried concurrently.
    app.state.cache = ResponseCache()
    app.state.client = httpx.AsyncClient(
        timeout=25.0, headers={"User-Agent": USER_AGENT}, follow_redirects=True
    )
    try:
        yield
    finally:
        await app.state.client.aclose()


app = FastAPI(
    title="Dissent",
    description="A variant interpretation engine that shows you the argument, not just the answer.",
    version="0.1.0",
    lifespan=lifespan,
)


@app.get("/api/examples")
async def examples() -> JSONResponse:
    return JSONResponse(EXAMPLES)


@app.get("/api/interpret", response_model=Report)
async def api_interpret(
    q: str = Query(..., min_length=1, max_length=200, description="rsID, HGVS, coordinates, or GENE p.Change"),
    assembly: str = Query("hg38", pattern="^(hg19|hg38)$"),
) -> Report:
    try:
        report = await asyncio.wait_for(
            interpret(q, assembly=assembly, cache=app.state.cache, client=app.state.client),
            # A hard ceiling so the UI never hangs. Individual sources have
            # their own budgets; this guards against all of them being slow.
            timeout=70.0,
        )
    except asyncio.TimeoutError:
        raise HTTPException(
            status_code=504,
            detail="Sources did not respond in time. Public APIs are rate-limited; try again.",
        ) from None

    if report.bundle is not None and not report.bundle.clinvar and not report.bundle.consequence:
        raise HTTPException(
            status_code=404,
            detail=(
                f"Could not resolve '{q}'. Try an rsID (rs113488022), genomic "
                "coordinates (chr7:140753336A>T), or GENE p.Change (BRAF V600E). "
                "Indels in protein notation such as CFTR F508del are not yet parsed. "
                "If the protein change came from clinical literature, many genes are "
                "still written in legacy numbering that omits the signal peptide while "
                "HGVS counts from the initiator methionine — Factor V Leiden is R506Q "
                "in the literature and p.Arg534Gln in HGVS. The two do not resolve to "
                "each other, so search the rsID or the c. notation."
            ),
        )
    return report


@app.get("/api/baseline")
async def api_baseline() -> JSONResponse:
    """Database-wide conflict rates, for the scale page.

    Serves the summary and the leaderboard, not the 15,326-row gene table:
    the page looks genes up one at a time.
    """
    baseline = load_baseline()
    if baseline is None:
        raise HTTPException(status_code=503, detail="Baseline data is not installed.")

    histogram = {}
    try:
        raw = json.loads((DATA_PATH).read_text())
        histogram = raw.get("submitters_per_variant", {})
    except (OSError, ValueError):
        histogram = {}

    return JSONResponse(
        {
            "generated": baseline.generated,
            "variants": baseline.variants,
            "submissions": baseline.submissions,
            "rates_pct": baseline.rates_pct,
            "counts": baseline.counts,
            "genes_with_a_rate": len(baseline.gene_table),
            "gene_min_variants": baseline.gene_min_variants,
            "leaderboard_min_variants": baseline.leaderboard_min_variants,
            "most_discordant_genes": baseline.most_discordant_genes,
            "submitters_per_variant": histogram,
        }
    )


@app.get("/api/baseline/opposed")
async def api_opposed(
    limit: int = Query(40, ge=1, le=200),
    gene: str | None = Query(None, max_length=40),
    min_stars: int = Query(0, ge=0, le=4),
) -> JSONResponse:
    """Variants one laboratory calls disease-causing and another calls harmless.

    The sharpest form the measurement takes: not a rate, but a list of
    identifiers anyone can open and check.
    """
    try:
        raw = json.loads(DATA_PATH.read_text())
    except (OSError, ValueError):
        raise HTTPException(status_code=503, detail="Baseline data is not installed.") from None

    rows = raw.get("opposed_variants", [])
    if gene:
        wanted = gene.strip().upper()
        rows = [r for r in rows if (r.get("gene") or "").upper() == wanted]
    if min_stars:
        rows = [r for r in rows if (r.get("best_stars") or 0) >= min_stars]

    return JSONResponse({"total": len(rows), "variants": rows[:limit]})


@app.get("/api/baseline/gene")
async def api_baseline_gene(symbol: str = Query(..., min_length=1, max_length=40)) -> JSONResponse:
    baseline = load_baseline()
    if baseline is None:
        raise HTTPException(status_code=503, detail="Baseline data is not installed.")

    rate = baseline.gene(symbol)
    if rate is None:
        raise HTTPException(
            status_code=404,
            detail=(
                f"No rate for '{symbol.upper()}'. A gene needs at least "
                f"{baseline.gene_min_variants} classified variants for the figure to "
                "mean anything."
            ),
        )
    rank = baseline.gene_rank(symbol)
    return JSONResponse(
        {
            "gene": rate.gene,
            "variants": rate.variants,
            "discordant": rate.discordant,
            "rate_pct": round(rate.rate_pct, 2),
            "rank": rank[0] if rank else None,
            "of_genes": rank[1] if rank else None,
            "overall_rate_pct": baseline.rates_pct.get(
                "any_label_disagreement_among_multi_submitter"
            ),
        }
    )


# --------------------------------------------------------------------------
# Panel triage
# --------------------------------------------------------------------------


@app.get("/api/panel/tiers")
async def api_panel_tiers() -> JSONResponse:
    """The triage rule, served so the UI states it rather than restating it.

    A ranking the reader cannot audit is worse than no ranking, so the
    definitions live in one place and both the page and the docs read them
    from here.
    """
    return JSONResponse(
        {
            "max_per_request": MAX_PANEL_SIZE,
            "max_query_length": MAX_QUERY_LENGTH,
            "evidence_not_public": EVIDENCE_NOT_PUBLIC,
            "tiers": [
                {"rank": t.rank, "key": t.key, "name": t.name, "meaning": t.meaning}
                for t in TIERS
            ],
            # The page says "the rule is below so you can disagree with it",
            # and served only the prose. These are the numbers the rule
            # actually turns on; without them a reader can disagree with the
            # description but not with the thing.
            "cutoffs": {
                "stale_after_years": DEFAULT_THRESHOLDS.stale_assertion_years,
                "thin_at_or_below_stars": 1,
                "no_criteria_ignored_at_or_above_stars": EXPERT_REVIEW_STARS,
                "pathogenicity_labels": sorted(s.value for s in AXIS_LABELS),
            },
        }
    )


@app.get("/api/panel/examples")
async def api_panel_examples() -> JSONResponse:
    return JSONResponse(PANEL_EXAMPLES)


@app.post("/api/panel", response_model=PanelResult)
async def api_panel(
    payload: dict = Body(...),
) -> PanelResult:
    """Interpret a chunk of a panel.

    The web UI sends a long list in chunks and renders rows as they land, so
    this returns a chunk's worth rather than streaming: a reader watching 40
    variants resolve wants to see the first ones, and a failed chunk should
    cost one chunk rather than the whole panel.
    """
    raw = payload.get("queries")
    if isinstance(raw, str):
        text = raw
    elif isinstance(raw, list):
        # Strings only. Stringifying whatever arrived turned {"a": 1} into the
        # query "{'a': 1}" and sent it to six public APIs.
        if not all(isinstance(q, str) for q in raw):
            raise HTTPException(
                status_code=422, detail="Every entry in 'queries' must be a string."
            )
        text = "\n".join(raw)
    else:
        raise HTTPException(
            status_code=422, detail="Send 'queries' as a list or a newline-separated string."
        )

    # Bound the input before parsing it. Without this a single request could
    # carry megabytes to be split and de-duplicated.
    if len(text) > MAX_QUERY_LENGTH * (MAX_PANEL_SIZE + 10):
        raise HTTPException(status_code=413, detail="That is more text than a panel can be.")

    queries, dropped = parse_panel_counted(text)
    if not queries:
        raise HTTPException(status_code=422, detail="No variants found in the input.")

    assembly = payload.get("assembly", "hg38")
    if assembly not in ("hg19", "hg38"):
        raise HTTPException(status_code=422, detail="assembly must be hg19 or hg38.")

    try:
        result = await asyncio.wait_for(
            interpret_panel(
                queries,
                assembly=assembly,
                cache=app.state.cache,
                client=app.state.client,
            ),
            # A ceiling for the whole chunk. Per-variant timeouts alone let a
            # full request hold a connection for ceil(n/concurrency) x 70
            # seconds, and the page's fetch has no timeout of its own, so the
            # progress bar would sit still rather than fail.
            timeout=PANEL_TIMEOUT_SECONDS,
        )
    except asyncio.TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=(
                "The panel did not finish in time. Public APIs are rate-limited; "
                "send fewer variants per request."
            ),
        ) from None

    result.dropped = dropped
    return result


@app.post("/api/panel/csv")
async def api_panel_csv(payload: dict = Body(...)) -> PlainTextResponse:
    """Turn finished rows into a spreadsheet.

    Takes rows the client already has rather than re-running the panel: a
    download should not cost the public APIs a second pass over the same
    variants, and the file then matches exactly what the reader is looking at.
    """
    raw_rows = payload.get("rows")
    if not isinstance(raw_rows, list) or not raw_rows:
        raise HTTPException(status_code=422, detail="Send 'rows' from a panel result.")
    # The endpoint validated each row's shape but never the number of them, so
    # a request could ask for an eight-megabyte file to be built in memory.
    if len(raw_rows) > MAX_ROWS_PER_EXPORT:
        raise HTTPException(
            status_code=413,
            detail=f"Send at most {MAX_ROWS_PER_EXPORT} rows per export.",
        )
    try:
        rows = [PanelRow.model_validate(r) for r in raw_rows]
    except Exception:  # noqa: BLE001
        raise HTTPException(status_code=422, detail="Rows did not match the panel schema.") from None

    body = rows_to_csv(sort_rows(rows))
    return PlainTextResponse(
        body,
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="dissent-panel.csv"'},
    )


@app.get("/panel")
async def panel_page() -> HTMLResponse:
    return _page("panel.html")


@app.get("/scale")
async def scale() -> HTMLResponse:
    return _page("scale.html")


@app.get("/api/health")
async def health() -> dict:
    return {
        "status": "ok",
        "assets": ASSET_VERSION,
        "cache": app.state.cache.stats(),
    }


@app.get("/")
async def index() -> HTMLResponse:
    return _page("index.html")


app.mount("/static", StaticFiles(directory=WEB_DIR), name="static")
