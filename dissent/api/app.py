"""FastAPI application: JSON API plus the static single-page UI."""

from __future__ import annotations

import asyncio
from contextlib import asynccontextmanager
import hashlib
import json
from pathlib import Path

import httpx
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from dissent.baseline.rates import load_baseline
from dissent.cache.store import ResponseCache
from dissent.models import Report
from dissent.report import interpret
from dissent.sources.base import USER_AGENT

WEB_DIR = Path(__file__).resolve().parent.parent / "web"


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
        "query": "BRCA1 C61G",
        "label": "BRCA1 C61G",
        "why": "Well-curated and unanimous. Shown so the conflict panel can be seen staying quiet.",
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
                "Indels in protein notation such as CFTR F508del are not yet parsed."
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
