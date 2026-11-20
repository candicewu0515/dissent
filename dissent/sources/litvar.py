"""NCBI LitVar2 client.

Returns the PubMed IDs that mention a given variant. This is deliberately NOT
used to assert PS3/BS3 (well-established functional studies): a PMID count says
nothing about whether an assay was performed or what it showed. It is used to
tell the user where to read, and to flag the case where a variant has
substantial literature that no ClinVar submitter has cited.
"""

from __future__ import annotations

from urllib.parse import quote

from dissent.models import SourceStatus
from dissent.sources.base import SourceClient

API = "https://www.ncbi.nlm.nih.gov/research/litvar2-api"
PUBMED_URL = "https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
LITVAR_URL = "https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?text={rsid}"


class LitVarClient(SourceClient):
    name = "litvar"

    async def fetch_pmids(self, rsid: str | None) -> tuple[list[int], SourceStatus]:
        if not rsid:
            return [], SourceStatus(source=self.name, ok=True, detail="no rsID to query")

        # LitVar2 identifies variants as `litvar@rs113488022##`. The id must be
        # percent-encoded before it goes in the path: an unescaped `#` is read
        # as the start of a URL fragment and the server returns HTTP 400.
        litvar_id = quote(f"litvar@{rsid}##", safe="")
        data, status = await self.fetch_json(
            f"{API}/variant/get/{litvar_id}/publications",
            cache_key_parts=["pubs", rsid],
        )
        if not isinstance(data, dict):
            return [], status

        pmids: list[int] = []
        for value in data.get("pmids") or []:
            try:
                pmids.append(int(value))
            except (TypeError, ValueError):
                continue
        return pmids, status
