"""Base rates for each conflict type, measured across the whole of ClinVar.

A report that says "submitters disagree about this variant" leaves the reader
without the one thing they need to act on it: is that unusual? The rates here
come from ClinVar's own `submission_summary.txt.gz` release file, one row per
submitted record, so submitter-level disagreement is directly observable rather
than inferred.

Measured on the 2026-11-15 release: 6,895,437 submissions across 4,470,690
variants. The headline figure is that 34.5% of variants with more than one
submitter have submitters who disagree, and that the rate is *highest* in the
genes where the answer matters most -- BRCA2 33.4%, TP53 33.1%, BRCA1 29.6%.

Regenerate with `python analysis/clinvar_baseline.py <release.gz> <out.json>`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "clinvar_baseline.json"


@dataclass(frozen=True)
class GeneRate:
    gene: str
    variants: int
    discordant: int

    @property
    def rate_pct(self) -> float:
        return 100.0 * self.discordant / self.variants if self.variants else 0.0


@dataclass(frozen=True)
class Baseline:
    generated: str
    variants: int
    submissions: int
    rates_pct: dict[str, float]
    counts: dict[str, int]
    gene_table: dict[str, list[int]]
    most_discordant_genes: list[dict]
    gene_min_variants: int
    leaderboard_min_variants: int

    def gene(self, symbol: str | None) -> GeneRate | None:
        """Disagreement rate for one gene, if enough variants exist to quote one."""
        if not symbol:
            return None
        row = self.gene_table.get(symbol.upper())
        if not row:
            return None
        return GeneRate(gene=symbol.upper(), variants=row[0], discordant=row[1])

    def gene_rank(self, symbol: str | None) -> tuple[int, int] | None:
        """(rank, total) of a gene by disagreement rate, 1 = most discordant.

        Ranked against every gene with a quotable rate, which is a larger and
        noisier population than the leaderboard's. The leaderboard applies a
        higher variant floor so its ordering is stable, so its top entry is
        not necessarily rank 1 here -- both populations are named wherever a
        figure from them is shown.

        Computed on demand rather than stored: the ranking is only ever needed
        for the one gene being looked at.
        """
        target = self.gene(symbol)
        if target is None:
            return None
        rates = sorted(
            (100.0 * d / v for v, d in self.gene_table.values() if v),
            reverse=True,
        )
        mine = target.rate_pct
        rank = sum(1 for r in rates if r > mine) + 1
        return (rank, len(rates))


@lru_cache(maxsize=1)
def load_baseline(path: Path | None = None) -> Baseline | None:
    """Load the measured base rates, or None if the data file is absent.

    Absence is not an error: the tool works without it, and simply omits the
    comparative context. Shipping a 471 KB data file should not be a condition
    of the thing running.
    """
    target = path or DATA_PATH
    try:
        raw = json.loads(target.read_text())
    except (OSError, json.JSONDecodeError):
        return None

    return Baseline(
        generated=raw.get("generated", ""),
        variants=raw.get("variants", 0),
        submissions=raw.get("submission_rows", 0),
        rates_pct=raw.get("rates_pct", {}),
        counts=raw.get("counts", {}),
        gene_table=raw.get("gene_table", {}),
        most_discordant_genes=raw.get("most_discordant_genes", []),
        gene_min_variants=raw.get("gene_min_variants_for_rate", 30),
        leaderboard_min_variants=raw.get("leaderboard_min_variants", 300),
    )
