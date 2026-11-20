"""MyVariant.info client.

The richest single source: it carries ClinVar RCV records, the dbNSFP predictor
panel, gnomAD allele frequencies, CADD and SpliceAI, in one request.

Assembly handling is explicit. MyVariant's default index is hg19; hg38 requires
`assembly=hg38`. Guessing here silently returns the wrong variant, so the
assembly travels with the query and is always sent.
"""

from __future__ import annotations

import re
from datetime import date, datetime
from typing import Any, NamedTuple

from dissent.models import (
    ClinVarAssertion,
    ClinVarEvidence,
    PopulationFrequency,
    PredictorCall,
    PredictorResult,
    SourceStatus,
    VariantQuery,
    normalize_significance,
    review_stars,
)
from dissent.sources.base import SourceClient, as_list, dig, first, to_float, to_float_all, to_int

API = "https://myvariant.info/v1"

# Reverse of the three-letter map, for rebuilding ClinVar-style protein names.
# ClinVar's `preferred_name` uses three-letter codes, e.g.
#   NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)
AA1_TO_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "U": "Sec",
}


def three_letter_protein(protein_change: str | None) -> str | None:
    """`p.C61G` -> `p.Cys61Gly`. Returns None if it cannot be rewritten.

    Needed because protein-position matching against dbNSFP's `hgvsp` field is
    unsafe: that field lists every transcript's numbering, so a query for
    residue 61 can match a variant that is residue 61 only on a non-canonical
    transcript. ClinVar's `preferred_name` is MANE-based and unambiguous.
    """
    if not protein_change:
        return None
    body = protein_change.removeprefix("p.")
    m = re.match(r"^([A-Z*])(\d+)([A-Z*])$", body)
    if not m:
        return None
    ref, pos, alt = m.groups()
    ref3, alt3 = AA1_TO_3.get(ref), AA1_TO_3.get(alt)
    if not ref3 or not alt3:
        return None
    return f"p.{ref3}{pos}{alt3}"

CLINVAR_URL = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/"
RCV_URL = "https://www.ncbi.nlm.nih.gov/clinvar/{acc}/"


class PredictorSpec(NamedTuple):
    """How to read one predictor out of dbNSFP.

    `codes` maps that predictor's own categorical letters to a call. These
    mappings are NOT interchangeable between predictors: `P` means "possibly
    damaging" in PolyPhen-2, "Pathogenic" in AlphaMissense, and "polymorphism"
    (i.e. benign) in MutationTaster. Each predictor therefore carries its own map.

    `threshold` is used only when the source gives no categorical call:
    (cutoff, direction) where direction says which side is damaging.
    """

    name: str
    path: tuple[str, ...]
    score_key: str | None = "score"
    pred_key: str | None = "pred"
    codes: dict[str, PredictorCall] = {}
    threshold: tuple[float, str] | None = None  # (cutoff, "above"|"below")


_D_T = {"D": "damaging", "T": "benign"}

PREDICTORS: tuple[PredictorSpec, ...] = (
    PredictorSpec("SIFT", ("dbnsfp", "sift"), codes=_D_T),
    PredictorSpec("SIFT4G", ("dbnsfp", "sift4g"), codes=_D_T),
    PredictorSpec(
        "PolyPhen-2 HVAR",
        ("dbnsfp", "polyphen2", "hvar"),
        codes={"D": "damaging", "P": "uncertain", "B": "benign"},
    ),
    PredictorSpec(
        "AlphaMissense",
        ("dbnsfp", "alphamissense"),
        codes={"P": "damaging", "A": "uncertain", "B": "benign"},
    ),
    PredictorSpec("REVEL", ("dbnsfp", "revel"), pred_key=None, threshold=(0.5, "above")),
    PredictorSpec("MetaSVM", ("dbnsfp", "metasvm"), codes=_D_T),
    PredictorSpec("MetaLR", ("dbnsfp", "metalr"), codes=_D_T),
    PredictorSpec("MetaRNN", ("dbnsfp", "metarnn"), codes=_D_T),
    PredictorSpec("ClinPred", ("dbnsfp", "clinpred"), codes=_D_T),
    PredictorSpec("M-CAP", ("dbnsfp", "m-cap"), codes=_D_T),
    PredictorSpec("PrimateAI", ("dbnsfp", "primateai"), codes=_D_T),
    PredictorSpec("DEOGEN2", ("dbnsfp", "deogen2"), codes=_D_T),
    PredictorSpec("LIST-S2", ("dbnsfp", "list-s2"), codes=_D_T),
    PredictorSpec("FATHMM", ("dbnsfp", "fathmm"), codes=_D_T),
    PredictorSpec("PROVEAN", ("dbnsfp", "provean"), codes={"D": "damaging", "N": "benign"}),
    PredictorSpec(
        "MutationTaster",
        ("dbnsfp", "mutationtaster"),
        codes={"A": "damaging", "D": "damaging", "N": "benign", "P": "benign"},
    ),
    PredictorSpec(
        "MutationAssessor",
        ("dbnsfp", "mutationassessor"),
        codes={"H": "damaging", "M": "damaging", "L": "benign", "N": "benign"},
    ),
    PredictorSpec("VEST4", ("dbnsfp", "vest4"), pred_key=None, threshold=(0.5, "above")),
    PredictorSpec(
        "EVE",
        ("dbnsfp", "eve"),
        pred_key="class",
        codes={"P": "damaging", "B": "benign", "U": "uncertain"},
    ),
    PredictorSpec("ESM1b", ("dbnsfp", "esm1b"), pred_key=None, threshold=(-7.5, "below")),
    PredictorSpec(
        "BayesDel", ("dbnsfp", "bayesdel", "add_af"), pred_key=None, threshold=(0.0692, "above")
    ),
)


def _parse_date(value: Any) -> date | None:
    raw = first(value)
    if not raw:
        return None
    for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%d %b %Y", "%Y"):
        try:
            return datetime.strptime(str(raw).strip(), fmt).date()
        except ValueError:
            continue
    return None


def _extreme(scores: list[float], direction: str) -> float | None:
    """Pick the most damaging score across transcripts.

    Reporting the mean would hide a transcript on which the variant is severe;
    reporting the first is arbitrary. The most damaging value is the one that
    matters for a worksheet a human will review.
    """
    if not scores:
        return None
    return max(scores) if direction == "above" else min(scores)


def extract_predictors(doc: dict[str, Any]) -> list[PredictorResult]:
    """Read the dbNSFP predictor panel out of a MyVariant document."""
    results: list[PredictorResult] = []

    for spec in PREDICTORS:
        node = doc
        for key in spec.path:
            node = node.get(key) if isinstance(node, dict) else None
            if node is None:
                break
        if node is None:
            continue
        if not isinstance(node, dict):
            # Some entries are a bare score rather than an object.
            node = {"score": node}

        direction = spec.threshold[1] if spec.threshold else "above"
        scores = to_float_all(node.get(spec.score_key)) if spec.score_key else []
        score = _extreme(scores, direction)
        rankscore = to_float(node.get("rankscore") or node.get("converted_rankscore"))

        call: PredictorCall | None = None
        from_threshold = False

        if spec.pred_key:
            codes = [str(c).strip().upper()[:1] for c in as_list(node.get(spec.pred_key)) if c]
            calls = [spec.codes.get(c) for c in codes]
            calls = [c for c in calls if c]
            if calls:
                # Any transcript called damaging makes the predictor damaging:
                # a variant severe on one transcript is a real finding.
                if "damaging" in calls:
                    call = "damaging"
                elif "uncertain" in calls:
                    call = "uncertain"
                else:
                    call = "benign"

        if call is None and spec.threshold and score is not None:
            cutoff, side = spec.threshold
            call = (
                "damaging"
                if (score >= cutoff if side == "above" else score <= cutoff)
                else "benign"
            )
            from_threshold = True

        if call is None and score is None and rankscore is None:
            continue

        results.append(
            PredictorResult(
                name=spec.name,
                score=score,
                rankscore=rankscore,
                call=call,
                call_from_threshold=from_threshold,
            )
        )

    return results


_ID_LOCUS = re.compile(r"^(chr[0-9XYMT]+):g\.(\d+)", re.I)


def _locus(variant_id: Any) -> tuple[str, str] | None:
    """`chr7:g.140753336A>T` -> `('CHR7', '140753336')`."""
    match = _ID_LOCUS.match(str(variant_id or ""))
    if not match:
        return None
    return (match.group(1).upper(), match.group(2))


def _same_position_siblings(best: dict[str, Any], hits: list[dict[str, Any]]) -> list[str]:
    """Other ALT alleles reported at the same locus as the chosen hit."""
    anchor = _locus(best.get("_id"))
    if anchor is None:
        return []
    out: list[str] = []
    for hit in hits:
        hit_id = hit.get("_id")
        if not hit_id or hit_id == best.get("_id"):
            continue
        if _locus(hit_id) == anchor:
            out.append(str(hit_id))
    return out


def _hit_genes(hit: dict[str, Any]) -> set[str]:
    """Every gene symbol a MyVariant document mentions, upper-cased."""
    genes: set[str] = set()
    cv = hit.get("clinvar")
    if isinstance(cv, dict):
        for g in as_list(cv.get("gene")):
            if isinstance(g, dict) and g.get("symbol"):
                genes.add(str(g["symbol"]).upper())
            elif isinstance(g, str):
                genes.add(g.upper())
    nsfp = hit.get("dbnsfp")
    if isinstance(nsfp, dict):
        for g in as_list(nsfp.get("genename")):
            if g:
                genes.add(str(g).upper())
    snpeff = dig(hit, "snpeff", "ann", "genename")
    for g in as_list(snpeff):
        if g:
            genes.add(str(g).upper())
    return genes


def _hit_preferred_names(hit: dict[str, Any]) -> list[str]:
    cv = hit.get("clinvar")
    if not isinstance(cv, dict):
        return []
    out: list[str] = []
    for rcv in as_list(cv.get("rcv")):
        if isinstance(rcv, dict):
            for name in as_list(rcv.get("preferred_name")):
                if name:
                    out.append(str(name))
    return out


def _preferred_name(cv: dict[str, Any]) -> str | None:
    """The MANE-based display name, e.g. `NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)`.

    Shown in the UI so the user can confirm at a glance that the tool resolved
    the variant they meant.
    """
    for rcv in as_list(cv.get("rcv")):
        if isinstance(rcv, dict):
            name = first(rcv.get("preferred_name"))
            if name:
                return str(name)
    return None


def extract_clinvar(doc: dict[str, Any]) -> ClinVarEvidence | None:
    cv = doc.get("clinvar")
    if not isinstance(cv, dict):
        return None

    assertions: list[ClinVarAssertion] = []
    for rcv in as_list(cv.get("rcv")):
        if not isinstance(rcv, dict):
            continue
        raw_sig = first(rcv.get("clinical_significance"))
        status = first(rcv.get("review_status"))
        condition = dig(rcv, "conditions", "name") or first(
            [c.get("name") for c in as_list(rcv.get("conditions")) if isinstance(c, dict)]
        )
        assertions.append(
            ClinVarAssertion(
                accession=first(rcv.get("accession")),
                significance_raw=raw_sig,
                significance=normalize_significance(raw_sig),
                review_status=status,
                stars=review_stars(status),
                number_submitters=to_int(rcv.get("number_submitters")) or 0,
                last_evaluated=_parse_date(rcv.get("last_evaluated")),
                condition=condition,
                origin=first(rcv.get("origin")),
            )
        )

    gene = dig(cv, "gene", "symbol") or first(
        [g.get("symbol") for g in as_list(cv.get("gene")) if isinstance(g, dict)]
    )

    return ClinVarEvidence(
        variant_id=to_int(cv.get("variant_id")),
        allele_id=to_int(cv.get("allele_id")),
        rsid=(f"rs{cv['rsid']}" if str(cv.get("rsid", "")).isdigit() else first(cv.get("rsid"))),
        gene=gene,
        preferred_name=_preferred_name(cv),
        assertions=assertions,
    )


_POP_LABELS = {
    "af_afr": "African/African-American",
    "af_amr": "Latino/Admixed American",
    "af_asj": "Ashkenazi Jewish",
    "af_eas": "East Asian",
    "af_fin": "Finnish",
    "af_nfe": "Non-Finnish European",
    "af_sas": "South Asian",
    "af_oth": "Other",
    "af_mid": "Middle Eastern",
}


def extract_frequencies(doc: dict[str, Any]) -> list[PopulationFrequency]:
    """Frequencies embedded in the MyVariant document.

    This is the fallback for when the live gnomAD query fails, and it is a
    weaker source than it looks: the embedded copy is gnomAD v2, and its
    per-population fields carry a frequency with no denominator. Without an
    allele number there is no way to tell 20% of five alleles from 20% of
    fifty thousand.

    So no popmax is derived here. Taking an unguarded maximum over those
    fields let a handful of alleles in `af_oth` drive a BA1 decision -- a
    false-benign path, and the worst direction to be wrong in. The global
    frequency, which does carry an allele number, is used instead, and the
    per-population figures are kept for display only.
    """
    out: list[PopulationFrequency] = []
    for key, label in (
        ("gnomad_exome", "gnomAD v2 exomes (via MyVariant)"),
        ("gnomad_genome", "gnomAD v2 genomes (via MyVariant)"),
    ):
        node = doc.get(key)
        if not isinstance(node, dict):
            continue
        af_node = node.get("af") if isinstance(node.get("af"), dict) else {}
        per_pop: dict[str, float] = {}
        for pop_key, pop_label in _POP_LABELS.items():
            value = to_float(af_node.get(pop_key))
            if value is not None:
                per_pop[pop_label] = value

        out.append(
            PopulationFrequency(
                source=label,
                af=to_float(af_node.get("af")),
                ac=to_int(dig(node, "ac", "ac")),
                an=to_int(dig(node, "an", "an")),
                homozygote_count=to_int(dig(node, "hom", "hom")),
                # Deliberately unset: see the docstring.
                popmax_af=None,
                popmax_population=None,
                per_population=per_pop,
            )
        )
    return out


class MyVariantClient(SourceClient):
    name = "myvariant"

    async def fetch(self, query: VariantQuery) -> tuple[dict[str, Any] | None, SourceStatus]:
        """Resolve a query to a single MyVariant document."""
        fields = "all"
        assembly = query.assembly

        if query.hgvs_g or (query.chrom and query.pos):
            hgvs = query.hgvs_g or f"chr{query.chrom}:g.{query.pos}{query.ref}>{query.alt}"
            params: dict[str, Any] = {"fields": fields}
            if assembly == "hg38":
                params["assembly"] = "hg38"
            doc, status = await self.fetch_json(
                f"{API}/variant/{hgvs}",
                params=params,
                cache_key_parts=["variant", hgvs, assembly],
            )
            if isinstance(doc, dict):
                return doc, status
            if not status.ok:
                return None, status

        terms = self._search_terms(query)
        if not terms:
            return None, SourceStatus(source=self.name, ok=True, detail="nothing to query")

        best: dict[str, Any] | None = None
        hits: list[dict[str, Any]] = []
        status = SourceStatus(source=self.name, ok=True, detail="no matching variant")

        # Candidate queries are tried in order of precision. The first is
        # constrained to the gene; only if that finds nothing do we fall back to
        # a free-text phrase, which ranks by relevance across the whole index.
        for term in terms:
            params = {"q": term, "fields": fields, "size": 25}
            if assembly == "hg38":
                params["assembly"] = "hg38"
            data, status = await self.fetch_json(
                f"{API}/query", params=params, cache_key_parts=["query", term, assembly]
            )
            if not isinstance(data, dict):
                continue
            hits = [h for h in data.get("hits", []) if isinstance(h, dict)]
            if not hits:
                continue
            best = self._best_hit(hits, query)
            if best is not None:
                break

        if best is None:
            return None, SourceStatus(
                source=self.name,
                ok=True,
                detail=(
                    f"no variant in {query.gene} matching {query.protein_change}"
                    if query.gene and query.protein_change
                    else "no matching variant"
                ),
                cached=status.cached,
            )
        # Record siblings so the report can warn that the input was ambiguous.
        # Only alleles at the *same genomic position* count: a phrase search for
        # `p.Val600Glu` also matches that residue change in unrelated genes, and
        # those are different variants, not alternative alleles of this one.
        best = dict(best)
        best["_siblings"] = _same_position_siblings(best, hits)
        return best, status

    @staticmethod
    def _search_terms(query: VariantQuery) -> list[str]:
        """Candidate queries, most precise first.

        For a gene-plus-protein query the gene-constrained form must come
        first. A bare phrase search for `"p.Arg273His"` returns 88 hits ranked
        by relevance across every gene that has that residue change, and TP53
        -- the gene actually asked for -- is not among the top matches. The
        most common TP53 hotspot in human cancer failed to resolve this way.

        The phrase search is kept as a fallback because it reaches ClinVar's
        `preferred_name`, which is MANE-numbered, whereas dbNSFP's `hgvsp`
        lists every transcript's numbering. Hits from either are verified
        against the requested gene and protein change by `_best_hit`.
        """
        if query.rsid:
            return [f"dbsnp.rsid:{query.rsid}"]

        if query.gene and query.protein_change:
            terms: list[str] = []
            three = three_letter_protein(query.protein_change)
            aa = query.protein_change.removeprefix("p.")
            if three:
                terms.append(f'dbnsfp.genename:{query.gene} AND dbnsfp.hgvsp:"{three}"')
            terms.append(f'dbnsfp.genename:{query.gene} AND dbnsfp.hgvsp:"p.{aa}"')
            if three:
                terms.append(f'"{three}"')
            return terms

        if query.transcript and query.hgvs_c:
            return [f'clinvar.hgvs.coding:"{query.transcript}:{query.hgvs_c}"']

        return [query.raw] if query.raw else []

    @staticmethod
    def _best_hit(hits: list[dict[str, Any]], query: VariantQuery) -> dict[str, Any] | None:
        """Prefer a hit that actually carries interpretable evidence.

        MyVariant's relevance order puts textual matches first, which can rank a
        bare record above one with ClinVar and predictor data for the same locus.
        """

        expected = three_letter_protein(query.protein_change)
        want_gene = (query.gene or "").upper()

        def score(hit: dict[str, Any]) -> tuple[int, int, int, int, float]:
            """Rank candidate hits.

            Ordering matters: a phrase search for `p.Cys61Gly` matches that
            residue change in every gene that has one. Gene agreement is
            therefore the first key, and an exact preferred_name match the
            second. Relevance score is only the tie-breaker.
            """
            gene_match = 1 if want_gene and want_gene in _hit_genes(hit) else 0

            name_match = 0
            if expected:
                for name in _hit_preferred_names(hit):
                    if expected.lower() in name.lower():
                        name_match = 1
                        break

            # Number of ClinVar assertions, not merely their presence: for a
            # multi-allelic rsID this is what distinguishes the allele people
            # have actually interpreted from its untouched neighbours.
            n_assertions = 0
            cv = hit.get("clinvar")
            if isinstance(cv, dict):
                n_assertions = len(as_list(cv.get("rcv")))
            has_dbnsfp = 1 if isinstance(hit.get("dbnsfp"), dict) else 0
            return (gene_match, name_match, n_assertions, has_dbnsfp, float(hit.get("_score") or 0))

        best = max(hits, key=score)
        # Refuse a confidently-wrong answer: if the caller named a gene and no
        # hit belongs to it, report nothing rather than a variant in some other
        # gene that shares the residue number.
        if want_gene and want_gene not in _hit_genes(best):
            return None
        return best
