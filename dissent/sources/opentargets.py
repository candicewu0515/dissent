"""Open Targets Platform client.

Used for gene-level context only: what the gene does and which diseases it is
associated with. This feeds the report's orientation panel, not the ACMG
criteria — association evidence is gene-level and cannot support a
variant-level criterion.
"""

from __future__ import annotations

from dissent.models import SourceStatus
from dissent.sources.base import SourceClient, first

API = "https://api.platform.opentargets.org/api/v4/graphql"
TARGET_URL = "https://platform.opentargets.org/target/{ensembl_id}"

_SEARCH = """
query Search($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id name }
  }
}
"""

_TARGET = """
query Target($id: String!) {
  target(ensemblId: $id) {
    id
    approvedSymbol
    approvedName
    biotype
    functionDescriptions
    associatedDiseases(page: {index: 0, size: 8}) {
      rows { score disease { id name } }
    }
  }
}
"""


def _graphql_errors(data: object) -> str | None:
    """Return a joined error string if a GraphQL payload carries one.

    GraphQL reports schema and validation failures with HTTP 200 and an
    `errors` array. Without this check a renamed field -- exactly the kind of
    change gnomAD made moving to v4 -- reads as "this variant is not in the
    database": every common variant would silently satisfy PM2, the frequency
    criteria would go unresolved, and the sources panel would still say `ok`.
    """
    if not isinstance(data, dict):
        return None
    errors = data.get("errors")
    if not errors:
        return None
    messages = []
    for err in errors if isinstance(errors, list) else [errors]:
        if isinstance(err, dict) and err.get("message"):
            messages.append(str(err["message"]))
        elif err:
            messages.append(str(err))
    return "; ".join(messages)[:300] if messages else "GraphQL error"


class OpenTargetsClient(SourceClient):
    name = "opentargets"

    async def fetch_gene(self, symbol: str) -> tuple[dict | None, SourceStatus]:
        if not symbol:
            return None, SourceStatus(source=self.name, ok=True, detail="no gene symbol")

        found, status = await self.fetch_json(
            API,
            method="POST",
            json_body={"query": _SEARCH, "variables": {"q": symbol}},
            cache_key_parts=["search", symbol],
        )
        if not isinstance(found, dict):
            return None, status

        error = _graphql_errors(found)
        if error:
            return None, SourceStatus(source=self.name, ok=False, detail=error)

        hits = (((found.get("data") or {}).get("search") or {}).get("hits")) or []
        ensembl_id = first([h.get("id") for h in hits if isinstance(h, dict)])
        if not ensembl_id:
            return None, SourceStatus(
                source=self.name, ok=True, detail="gene not found", cached=status.cached
            )

        data, status2 = await self.fetch_json(
            API,
            method="POST",
            json_body={"query": _TARGET, "variables": {"id": ensembl_id}},
            cache_key_parts=["target", ensembl_id],
        )
        if not isinstance(data, dict):
            return None, status2

        error = _graphql_errors(data)
        if error:
            return None, SourceStatus(source=self.name, ok=False, detail=error)

        target = (data.get("data") or {}).get("target")
        if not isinstance(target, dict):
            return None, SourceStatus(
                source=self.name, ok=True, detail="no target record", cached=status2.cached
            )

        diseases = [
            {"name": row["disease"]["name"], "score": row.get("score")}
            for row in ((target.get("associatedDiseases") or {}).get("rows") or [])
            if isinstance(row, dict) and isinstance(row.get("disease"), dict)
        ]

        return (
            {
                "ensembl_id": target.get("id"),
                "symbol": target.get("approvedSymbol"),
                "name": target.get("approvedName"),
                "biotype": target.get("biotype"),
                "function": first(target.get("functionDescriptions")),
                "diseases": diseases,
                "url": TARGET_URL.format(ensembl_id=target.get("id")),
            },
            status2,
        )
