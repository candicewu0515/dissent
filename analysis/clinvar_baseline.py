"""Compute how common each conflict type is across the whole of ClinVar.

Turns "this variant has disagreement" into "this variant's disagreement is in
the worst N% of ClinVar". A claim about one variant is an anecdote; a base rate
makes it a measurement.

Input is ClinVar's own submission_summary.txt.gz, which has one row per
submitted record (SCV) rather than per variant, so submitter-level disagreement
is directly observable. Streamed from the gzip; never fully decompressed.
"""

from __future__ import annotations

import gzip
import json
import sys
from collections import Counter, defaultdict
from datetime import date, datetime

# Column names in the header line (which begins with '#').
WANT = [
    "VariationID",
    "ClinicalSignificance",
    "DateLastEvaluated",
    "ReviewStatus",
    "Submitter",
    "SubmittedGeneSymbol",
]

PATHOGENIC = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
BENIGN = {"benign", "likely benign", "benign/likely benign"}
UNCERTAIN = {"uncertain significance", "uncertain risk allele"}
PLACEHOLDER = {"not provided", "no classification provided", "-", ""}

NO_CRITERIA = {
    "no assertion criteria provided",
    "no assertion provided",
    "no classification provided",
    "no classifications from unflagged records",
    # A flagged submission has been withdrawn from the aggregate by ClinVar;
    # it carries no review weight.
    "flagged submission",
}

# Staleness is measured against the day the analysis runs, not a fixed date,
# so a re-run against a newer release reports the rate correctly.
TODAY = date.today()
STALE_YEARS = 5


def parse_date(value: str) -> date | None:
    value = (value or "").strip()
    if not value or value == "-":
        return None
    for fmt in ("%b %d, %Y", "%Y-%m-%d", "%d %b %Y"):
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            continue
    return None


def main(path: str, out_path: str) -> None:
    # Accumulate per variant, streaming. Holding every row would need several GB.
    sigs: dict[str, list[str]] = defaultdict(list)
    stars: dict[str, list[int]] = defaultdict(list)
    newest: dict[str, date] = {}
    submitters: dict[str, set[str]] = defaultdict(set)
    genes: dict[str, str] = {}

    header: list[str] | None = None
    idx: dict[str, int] = {}
    rows = 0

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                # The last commented line is the header.
                candidate = line.lstrip("#").rstrip("\n").split("\t")
                if "VariationID" in candidate:
                    header = candidate
                    idx = {name: header.index(name) for name in WANT if name in header}
                continue
            if header is None:
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            rows += 1

            vid = parts[idx["VariationID"]]
            sig = parts[idx["ClinicalSignificance"]].strip().lower()
            review = parts[idx["ReviewStatus"]].strip().lower()

            sigs[vid].append(sig)
            stars[vid].append(0 if review in NO_CRITERIA else 1)

            d = parse_date(parts[idx["DateLastEvaluated"]]) if "DateLastEvaluated" in idx else None
            if d and (vid not in newest or d > newest[vid]):
                newest[vid] = d

            if "Submitter" in idx:
                submitters[vid].add(parts[idx["Submitter"]])
            if "SubmittedGeneSymbol" in idx and vid not in genes:
                gene = parts[idx["SubmittedGeneSymbol"]].strip()
                if gene and gene != "-":
                    genes[vid] = gene

            if rows % 1_000_000 == 0:
                print(f"  {rows:,} rows, {len(sigs):,} variants", file=sys.stderr, flush=True)

    print(f"read {rows:,} submission rows covering {len(sigs):,} variants", file=sys.stderr)

    total = 0
    multi_submitter = 0
    discordant_pb = 0        # pathogenic-side and benign-side both present
    discordant_any = 0       # any disagreement among actionable labels
    vus_vs_path = 0
    all_zero_star = 0
    any_zero_star = 0
    stale = 0
    dated = 0
    per_gene_discordant: Counter[str] = Counter()
    per_gene_total: Counter[str] = Counter()
    submitter_counts: Counter[int] = Counter()

    for vid, labels in sigs.items():
        total += 1
        actionable = [s for s in labels if s not in PLACEHOLDER]
        n_sub = len(submitters.get(vid, ())) or len(labels)
        submitter_counts[min(n_sub, 10)] += 1
        if n_sub > 1:
            multi_submitter += 1

        gene = genes.get(vid)
        if gene:
            per_gene_total[gene] += 1

        if len(set(actionable)) > 1:
            discordant_any += 1
            if gene:
                per_gene_discordant[gene] += 1

        has_p = any(s in PATHOGENIC for s in actionable)
        has_b = any(s in BENIGN for s in actionable)
        has_u = any(s in UNCERTAIN for s in actionable)
        if has_p and has_b:
            discordant_pb += 1
        if has_p and has_u:
            vus_vs_path += 1

        st = stars[vid]
        if st and max(st) == 0:
            all_zero_star += 1
        if st and min(st) == 0:
            any_zero_star += 1

        d = newest.get(vid)
        if d:
            dated += 1
            if (TODAY - d).days / 365.25 >= STALE_YEARS:
                stale += 1

    def pct(n: int, d: int) -> float:
        return round(100.0 * n / d, 2) if d else 0.0

    # Genes where disagreement is most concentrated, among those with enough
    # variants for the rate to mean anything.
    gene_rates = [
        {
            "gene": g,
            "variants": per_gene_total[g],
            "discordant": per_gene_discordant[g],
            "rate_pct": pct(per_gene_discordant[g], per_gene_total[g]),
        }
        for g in per_gene_total
        if per_gene_total[g] >= 300
    ]
    gene_rates.sort(key=lambda r: r["rate_pct"], reverse=True)

    # A compact lookup for every gene with enough variants to quote a rate.
    # The report needs the rate for whichever gene the user asked about, not
    # only for the worst offenders. Stored as [variants, discordant] pairs to
    # keep the file small enough to ship in the repository.
    MIN_FOR_RATE = 30
    gene_table = {
        g: [per_gene_total[g], per_gene_discordant[g]]
        for g in per_gene_total
        if per_gene_total[g] >= MIN_FOR_RATE
    }

    result = {
        "generated": TODAY.isoformat(),
        "source": "ClinVar submission_summary.txt.gz",
        "submission_rows": rows,
        "variants": total,
        "variants_with_multiple_submitters": multi_submitter,
        "rates_pct": {
            "multiple_submitters": pct(multi_submitter, total),
            "any_label_disagreement": pct(discordant_any, total),
            "any_label_disagreement_among_multi_submitter": pct(discordant_any, multi_submitter),
            "pathogenic_vs_benign": pct(discordant_pb, total),
            "pathogenic_vs_uncertain": pct(vus_vs_path, total),
            "every_record_lacks_criteria": pct(all_zero_star, total),
            "some_record_lacks_criteria": pct(any_zero_star, total),
            "stale_5y": pct(stale, dated),
        },
        "counts": {
            "any_label_disagreement": discordant_any,
            "pathogenic_vs_benign": discordant_pb,
            "pathogenic_vs_uncertain": vus_vs_path,
            "every_record_lacks_criteria": all_zero_star,
            "some_record_lacks_criteria": any_zero_star,
            "stale_5y": stale,
            "variants_with_a_date": dated,
        },
        "submitters_per_variant": {str(k): v for k, v in sorted(submitter_counts.items())},
        "most_discordant_genes": gene_rates[:40],
        "gene_min_variants_for_rate": MIN_FOR_RATE,
        "gene_table": gene_table,
    }

    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=1)
    print(json.dumps(result["rates_pct"], indent=1))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
