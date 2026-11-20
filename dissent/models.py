"""Typed models for the evidence bundle, ACMG results, and conflicts.

Design rule: every value that reaches the UI is reachable from one of these
models, and anything asserted carries a Provenance. Nothing is rendered that
cannot be traced back to a source record.
"""

from __future__ import annotations

from datetime import date
from enum import Enum
from typing import Literal

from pydantic import BaseModel, Field

Assembly = Literal["hg19", "hg38"]
QueryKind = Literal["rsid", "hgvs", "gcoord", "gene_protein", "unknown"]


# --------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------


class Provenance(BaseModel):
    """Where a single claim came from. Required on every assertion."""

    source: str
    label: str
    url: str | None = None
    retrieved_at: str | None = None


# --------------------------------------------------------------------------
# Query
# --------------------------------------------------------------------------


class VariantQuery(BaseModel):
    raw: str
    kind: QueryKind = "unknown"
    assembly: Assembly = "hg38"
    rsid: str | None = None
    hgvs_g: str | None = None
    hgvs_c: str | None = None
    transcript: str | None = None
    chrom: str | None = None
    pos: int | None = None
    ref: str | None = None
    alt: str | None = None
    gene: str | None = None
    protein_change: str | None = None

    @property
    def gnomad_variant_id(self) -> str | None:
        """gnomAD wants `7-140753336-A-T`."""
        if self.chrom and self.pos and self.ref and self.alt:
            c = self.chrom.removeprefix("chr")
            return f"{c}-{self.pos}-{self.ref}-{self.alt}"
        return None


# --------------------------------------------------------------------------
# ClinVar
# --------------------------------------------------------------------------


class Significance(str, Enum):
    PATHOGENIC = "pathogenic"
    LIKELY_PATHOGENIC = "likely_pathogenic"
    UNCERTAIN = "uncertain_significance"
    LIKELY_BENIGN = "likely_benign"
    BENIGN = "benign"
    CONFLICTING = "conflicting"
    NOT_PROVIDED = "not_provided"
    OTHER = "other"

    @property
    def is_pathogenic_side(self) -> bool:
        return self in (Significance.PATHOGENIC, Significance.LIKELY_PATHOGENIC)

    @property
    def is_benign_side(self) -> bool:
        return self in (Significance.BENIGN, Significance.LIKELY_BENIGN)

    @property
    def is_actionable(self) -> bool:
        """Carries a real interpretation, as opposed to a placeholder."""
        return self not in (Significance.NOT_PROVIDED, Significance.OTHER)


_SIGNIFICANCE_MAP: dict[str, Significance] = {
    "pathogenic": Significance.PATHOGENIC,
    "pathogenic/likely pathogenic": Significance.PATHOGENIC,
    "likely pathogenic": Significance.LIKELY_PATHOGENIC,
    "uncertain significance": Significance.UNCERTAIN,
    "uncertain risk allele": Significance.UNCERTAIN,
    "likely benign": Significance.LIKELY_BENIGN,
    "benign": Significance.BENIGN,
    "benign/likely benign": Significance.BENIGN,
    "conflicting interpretations of pathogenicity": Significance.CONFLICTING,
    "conflicting data from submitters": Significance.CONFLICTING,
    "not provided": Significance.NOT_PROVIDED,
    "no classification provided": Significance.NOT_PROVIDED,
}


def normalize_significance(raw: str | None) -> Significance:
    if not raw:
        return Significance.NOT_PROVIDED
    key = raw.strip().lower()
    if key in _SIGNIFICANCE_MAP:
        return _SIGNIFICANCE_MAP[key]
    # Compound labels such as "Pathogenic, low penetrance".
    head = key.split(",")[0].strip()
    return _SIGNIFICANCE_MAP.get(head, Significance.OTHER)


# ClinVar review status -> star rating. Mirrors ClinVar's own documented scale.
_REVIEW_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
}


def review_stars(review_status: str | None) -> int:
    if not review_status:
        return 0
    return _REVIEW_STARS.get(review_status.strip().lower(), 0)


class ClinVarAssertion(BaseModel):
    accession: str | None = None
    significance_raw: str | None = None
    significance: Significance = Significance.NOT_PROVIDED
    review_status: str | None = None
    stars: int = 0
    number_submitters: int = 0
    last_evaluated: date | None = None
    condition: str | None = None
    origin: str | None = None

    @property
    def has_assertion_criteria(self) -> bool:
        return self.stars >= 1


class ClinVarEvidence(BaseModel):
    variant_id: int | None = None
    allele_id: int | None = None
    rsid: str | None = None
    gene: str | None = None
    preferred_name: str | None = None
    assertions: list[ClinVarAssertion] = Field(default_factory=list)

    @property
    def actionable_assertions(self) -> list[ClinVarAssertion]:
        return [a for a in self.assertions if a.significance.is_actionable]

    @property
    def distinct_significances(self) -> set[Significance]:
        return {a.significance for a in self.actionable_assertions}

    @property
    def origins(self) -> set[str]:
        return {a.origin.lower() for a in self.assertions if a.origin}


# --------------------------------------------------------------------------
# Predictors
# --------------------------------------------------------------------------

PredictorCall = Literal["damaging", "benign", "uncertain"]


class PredictorResult(BaseModel):
    name: str
    score: float | None = None
    rankscore: float | None = None
    call: PredictorCall | None = None
    # True when `call` was derived from a numeric threshold rather than read
    # from the source's own categorical prediction.
    call_from_threshold: bool = False


# --------------------------------------------------------------------------
# Population frequency / constraint
# --------------------------------------------------------------------------


class PopulationFrequency(BaseModel):
    source: str
    af: float | None = None
    ac: int | None = None
    an: int | None = None
    homozygote_count: int | None = None
    # Highest AF across continental populations, and which one.
    popmax_af: float | None = None
    popmax_population: str | None = None
    per_population: dict[str, float] = Field(default_factory=dict)


class GeneConstraint(BaseModel):
    gene_symbol: str | None = None
    gene_id: str | None = None
    pli: float | None = None
    oe_lof: float | None = None
    oe_lof_upper: float | None = None  # LOEUF
    oe_mis: float | None = None
    mis_z: float | None = None
    lof_z: float | None = None


# --------------------------------------------------------------------------
# Consequence
# --------------------------------------------------------------------------

# Consequence terms that truncate or abolish the protein product.
LOF_CONSEQUENCES = frozenset(
    {
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "start_lost",
    }
)


class Consequence(BaseModel):
    most_severe: str | None = None
    transcript_id: str | None = None
    gene_symbol: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    protein_position: int | None = None
    amino_acids: str | None = None
    canonical: bool = False
    biotype: str | None = None
    terms: list[str] = Field(default_factory=list)

    @property
    def is_lof(self) -> bool:
        return bool(LOF_CONSEQUENCES & set(self.terms or ([self.most_severe] if self.most_severe else [])))

    @property
    def is_missense(self) -> bool:
        return "missense_variant" in (self.terms or []) or self.most_severe == "missense_variant"

    @property
    def is_synonymous(self) -> bool:
        return "synonymous_variant" in (self.terms or []) or self.most_severe == "synonymous_variant"

    @property
    def is_inframe_indel(self) -> bool:
        t = set(self.terms or [])
        return bool(t & {"inframe_insertion", "inframe_deletion"})


# --------------------------------------------------------------------------
# Source status — surfaced so the user knows what was NOT consulted
# --------------------------------------------------------------------------


class SourceStatus(BaseModel):
    source: str
    ok: bool
    detail: str | None = None
    cached: bool = False
    elapsed_ms: int | None = None


# --------------------------------------------------------------------------
# Evidence bundle
# --------------------------------------------------------------------------


class EvidenceBundle(BaseModel):
    query: VariantQuery
    gene_symbol: str | None = None
    clinvar: ClinVarEvidence | None = None
    predictors: list[PredictorResult] = Field(default_factory=list)
    frequencies: list[PopulationFrequency] = Field(default_factory=list)
    constraint: GeneConstraint | None = None
    # Primary consequence, from snpEff inside the MyVariant document.
    consequence: Consequence | None = None
    # Independent VEP call, kept separate so the two can be compared rather
    # than one silently overwriting the other.
    consequence_vep: Consequence | None = None
    # Capped for payload size; `literature_total` keeps the true count, which
    # the unreviewed-literature detector needs. BRAF V600E has over 32,000
    # linked publications and storing them all bloated a fixture to 400 KB.
    literature_pmids: list[int] = Field(default_factory=list)
    literature_total: int = 0
    gene_description: str | None = None
    splice_ai_max: float | None = None
    cadd_phred: float | None = None
    # Other ALT alleles sharing this rsID. An rsID is a *position* identifier,
    # not an allele identifier, so a multi-allelic rsID is ambiguous input.
    alternative_alleles: list[str] = Field(default_factory=list)
    sources: list[SourceStatus] = Field(default_factory=list)

    def best_frequency(self) -> PopulationFrequency | None:
        """The frequency record with the largest allele number.

        Largest AN means the best-powered estimate, which is what the
        BA1/BS1/PM2 thresholds assume.
        """
        with_an = [f for f in self.frequencies if f.an]
        if not with_an:
            return self.frequencies[0] if self.frequencies else None
        return max(with_an, key=lambda f: f.an or 0)

    def max_af(self) -> float | None:
        """Highest credible AF across sources, preferring popmax.

        BA1/BS1 are about the *most common* population, not the global mean.
        """
        candidates: list[float] = []
        for f in self.frequencies:
            if f.popmax_af is not None:
                candidates.append(f.popmax_af)
            elif f.af is not None:
                candidates.append(f.af)
        return max(candidates) if candidates else None


# --------------------------------------------------------------------------
# ACMG
# --------------------------------------------------------------------------


class Strength(str, Enum):
    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"
    STAND_ALONE = "stand_alone"


CriterionStatus = Literal["met", "not_met", "insufficient", "manual"]
Direction = Literal["pathogenic", "benign"]


class CriterionResult(BaseModel):
    code: str
    name: str
    direction: Direction
    status: CriterionStatus
    strength: Strength | None = None
    rationale: str
    sources: list[Provenance] = Field(default_factory=list)
    what_would_change: str | None = None


class Classification(BaseModel):
    call: str
    reasoning: str
    met_codes: list[str] = Field(default_factory=list)


# --------------------------------------------------------------------------
# Conflicts
# --------------------------------------------------------------------------

Severity = Literal["high", "medium", "low"]


class Conflict(BaseModel):
    kind: str
    severity: Severity
    headline: str
    detail: str
    sources: list[Provenance] = Field(default_factory=list)


class Report(BaseModel):
    query: VariantQuery
    gene_symbol: str | None = None
    display_name: str | None = None
    criteria: list[CriterionResult] = Field(default_factory=list)
    classification: Classification | None = None
    conflicts: list[Conflict] = Field(default_factory=list)
    bundle: EvidenceBundle | None = None
    generated_at: str | None = None
