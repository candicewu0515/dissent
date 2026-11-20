"""FastAPI application: JSON API plus the static single-page UI."""

from __future__ import annotations

import asyncio
from contextlib import asynccontextmanager
from pathlib import Path

import httpx
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from dissent.cache.store import ResponseCache
from dissent.models import Report
from dissent.report import interpret
from dissent.sources.base import USER_AGENT

WEB_DIR = Path(__file__).resolve().parent.parent / "web"

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


@app.get("/api/health")
async def health() -> dict:
    return {"status": "ok", "cache": app.state.cache.stats()}


@app.get("/")
async def index() -> FileResponse:
    return FileResponse(WEB_DIR / "index.html")


app.mount("/static", StaticFiles(directory=WEB_DIR), name="static")
