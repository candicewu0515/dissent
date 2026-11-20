"""One pure function per ACMG/AMP criterion.

Every function takes the evidence bundle and the thresholds, and returns a
CriterionResult. No function performs I/O or holds state, which is what makes
the engine testable offline and auditable by a reviewer.

Four statuses, and the distinction between the last two matters:

    met          the evidence supports applying this criterion
    not_met      the evidence was available and argues against it
    insufficient the evidence needed is missing from the sources consulted
    manual       this criterion cannot be decided from public data at all

A tool that silently omits what it cannot evaluate invites the reader to treat
absence as absence of evidence. All 28 criteria are therefore always returned,
and the ones requiring a human say what evidence would settle them.
"""

from __future__ import annotations

from datetime import date
from typing import Callable

from dissent.acmg.thresholds import DEFAULT, Thresholds
from dissent.models import (
    CriterionResult,
    EvidenceBundle,
    Provenance,
    Strength,
)

CriterionFn = Callable[[EvidenceBundle, Thresholds], CriterionResult]


def _fmt_af(af: float | None) -> str:
    if af is None:
        return "unknown"
    if af == 0:
        return "0"
    if af < 1e-4:
        return f"{af:.2e}"
    return f"{af * 100:.3g}%"


def _gnomad_provenance(bundle: EvidenceBundle) -> list[Provenance]:
    vid = bundle.query.gnomad_variant_id
    if not vid:
        return []
    return [
        Provenance(
            source="gnomAD",
            label=f"gnomAD v4 variant {vid}",
            url=f"https://gnomad.broadinstitute.org/variant/{vid}?dataset=gnomad_r4",
        )
    ]


def _clinvar_provenance(bundle: EvidenceBundle) -> list[Provenance]:
    cv = bundle.clinvar
    if cv is None or cv.variant_id is None:
        return []
    return [
        Provenance(
            source="ClinVar",
            label=f"ClinVar variation {cv.variant_id}",
            url=f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv.variant_id}/",
        )
    ]


def _gene_provenance(bundle: EvidenceBundle) -> list[Provenance]:
    if not bundle.gene_symbol:
        return []
    return [
        Provenance(
            source="gnomAD",
            label=f"gnomAD gene constraint for {bundle.gene_symbol}",
            url=f"https://gnomad.broadinstitute.org/gene/{bundle.gene_symbol}?dataset=gnomad_r4",
        )
    ]


def _manual(
    code: str,
    name: str,
    direction: str,
    strength: Strength,
    needs: str,
) -> CriterionResult:
    """A criterion that public data cannot decide."""
    return CriterionResult(
        code=code,
        name=name,
        direction=direction,  # type: ignore[arg-type]
        status="manual",
        strength=strength,
        rationale="Not determinable from the public sources consulted.",
        what_would_change=needs,
    )


# ==========================================================================
# Pathogenic — Very Strong
# ==========================================================================


def pvs1(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Null variant in a gene where loss of function causes disease"
    cons = bundle.consequence
    constraint = bundle.constraint

    if cons is None:
        return CriterionResult(
            code="PVS1", name=name, direction="pathogenic", status="insufficient",
            strength=Strength.VERY_STRONG,
            rationale="No molecular consequence was resolved for this variant.",
            what_would_change="A transcript-level consequence call.",
        )

    if not cons.is_lof:
        return CriterionResult(
            code="PVS1", name=name, direction="pathogenic", status="not_met",
            strength=Strength.VERY_STRONG,
            rationale=(
                f"Consequence is {cons.most_severe or 'unclear'}, which is not a "
                "null variant (nonsense, frameshift, canonical splice +/-1,2, or "
                "start loss)."
            ),
        )

    # The criterion requires that loss of function be a *known mechanism of
    # disease* for the gene. Public constraint metrics do not say that; they say
    # the gene is depleted of loss-of-function variation in the population,
    # which is suggestive but not the same claim.
    lof_intolerant = False
    detail = "no gene constraint data available"
    if constraint:
        pli, loeuf = constraint.pli, constraint.oe_lof_upper
        parts = []
        if pli is not None:
            parts.append(f"pLI {pli:.3g}")
        if loeuf is not None:
            parts.append(f"LOEUF {loeuf:.3g}")
        detail = ", ".join(parts) if parts else detail
        lof_intolerant = (pli is not None and pli >= t.pvs1_pli) or (
            loeuf is not None and loeuf <= t.pvs1_loeuf
        )

    if lof_intolerant:
        return CriterionResult(
            code="PVS1", name=name, direction="pathogenic", status="met",
            strength=Strength.VERY_STRONG,
            rationale=(
                f"Null variant ({cons.most_severe}) in {bundle.gene_symbol}, which is "
                f"strongly depleted of loss-of-function variation ({detail}). That "
                "depletion is a proxy for haploinsufficiency, not a statement of "
                "mechanism: the criterion asks whether loss of function is an "
                "established cause of this gene's disease, which needs curated "
                "gene-disease evidence. High constraint supports the inference; it "
                "does not complete it."
            ),
            sources=_gene_provenance(bundle),
            what_would_change=(
                "A ClinGen gene-disease validity record confirming loss of function "
                "as the mechanism, and confirmation the variant is not in the last "
                "exon or otherwise escaping nonsense-mediated decay."
            ),
        )

    # A low pLI is NOT evidence against PVS1, and treating it that way would
    # misclassify most recessive disease genes. pLI and LOEUF measure selection
    # against *heterozygous* loss of function. In a recessive condition carriers
    # are healthy, so null alleles persist in the population and the gene scores
    # as unconstrained. CFTR is the clearest example: pLI 3.8e-65, LOEUF 1.19,
    # while loss of function is unquestionably the mechanism of cystic fibrosis.
    # The same holds for late-onset tumour suppressors such as BRCA1.
    return CriterionResult(
        code="PVS1", name=name, direction="pathogenic", status="insufficient",
        strength=Strength.VERY_STRONG,
        rationale=(
            f"This is a null variant ({cons.most_severe}) in {bundle.gene_symbol}, which "
            f"satisfies the first half of PVS1. The gene is not loss-of-function "
            f"constrained ({detail}), but that is not evidence against the criterion: "
            "pLI and LOEUF measure selection against heterozygous loss of function, so "
            "recessive disease genes score as unconstrained because carriers are "
            "healthy. CFTR scores pLI 3.8e-65 while loss of function is certainly the "
            "mechanism of cystic fibrosis. PVS1 is therefore left open rather than "
            "declined, and needs the gene's mechanism confirmed."
        ),
        sources=_gene_provenance(bundle),
        what_would_change=(
            "The gene's established disease mechanism and inheritance pattern -- a "
            "ClinGen gene-disease validity record, or an OMIM entry stating whether "
            "loss of function is causal. Also confirm the variant is not in the last "
            "exon, where nonsense-mediated decay may be escaped."
        ),
    )


# ==========================================================================
# Pathogenic — Strong
# ==========================================================================


def ps1(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PS1",
        "Same amino acid change as an established pathogenic variant",
        "pathogenic",
        Strength.STRONG,
        "A ClinVar search for other nucleotide changes producing the same amino "
        "acid substitution, with at least one classified pathogenic on review "
        "status of two stars or better.",
    )


def ps2(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PS2", "De novo, with maternity and paternity confirmed", "pathogenic",
        Strength.STRONG,
        "Trio sequencing confirming the variant is absent in both biological parents.",
    )


def ps3(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    result = _manual(
        "PS3", "Well-established functional study shows a damaging effect",
        "pathogenic", Strength.STRONG,
        "A functional assay validated against known pathogenic and benign controls. "
        "A count of publications mentioning the variant is not evidence of one.",
    )
    total = bundle.literature_total or len(bundle.literature_pmids)
    if total:
        result.rationale = (
            f"{total} publications mention this variant, but "
            "whether any reports a validated functional assay cannot be determined "
            "without reading them."
        )
        result.sources = [
            Provenance(
                source="LitVar",
                label=f"{total} linked publications",
                url=f"https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?text={bundle.query.rsid}",
            )
        ] if bundle.query.rsid else []
    return result


def ps4(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PS4", "Prevalence in affected individuals significantly increased", "pathogenic",
        Strength.STRONG,
        "A case-control study, or enough independent probands for ClinGen's "
        "proband-counting alternative.",
    )


# ==========================================================================
# Pathogenic — Moderate
# ==========================================================================


def pm1(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PM1", "In a mutational hot spot or critical functional domain", "pathogenic",
        Strength.MODERATE,
        "Domain annotation for the residue plus the distribution of known "
        "pathogenic and benign variants across that domain.",
    )


def pm2(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Absent, or at extremely low frequency, in population databases"
    freq = bundle.best_frequency()
    af = bundle.max_af()

    # ClinGen SVI recommends PM2 at Supporting strength rather than Moderate.
    strength = Strength.SUPPORTING

    if freq is None or freq.an is None:
        return CriterionResult(
            code="PM2", name=name, direction="pathogenic", status="insufficient",
            strength=strength,
            rationale="No population frequency data was returned for this variant.",
            sources=_gnomad_provenance(bundle),
            what_would_change="A gnomAD record, or a note that the site is not callable there.",
        )

    if freq.an < t.min_allele_number:
        return CriterionResult(
            code="PM2", name=name, direction="pathogenic", status="insufficient",
            strength=strength,
            rationale=(
                f"Only {freq.an:,} alleles were genotyped in {freq.source}, below the "
                f"{t.min_allele_number:,} needed for a frequency statement."
            ),
            sources=_gnomad_provenance(bundle),
            what_would_change="Frequency from a better-powered cohort.",
        )

    if af is None or af <= t.pm2_af:
        return CriterionResult(
            code="PM2", name=name, direction="pathogenic", status="met",
            strength=strength,
            rationale=(
                f"Highest population frequency is {_fmt_af(af)} in {freq.source} "
                f"(AN {freq.an:,}), at or below the {_fmt_af(t.pm2_af)} threshold. "
                "Applied at Supporting per ClinGen SVI rather than Moderate."
            ),
            sources=_gnomad_provenance(bundle),
        )

    return CriterionResult(
        code="PM2", name=name, direction="pathogenic", status="not_met",
        strength=strength,
        rationale=(
            f"Population frequency {_fmt_af(af)} exceeds the {_fmt_af(t.pm2_af)} "
            f"threshold ({freq.popmax_population or freq.source})."
        ),
        sources=_gnomad_provenance(bundle),
    )


def pm3(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PM3", "For recessive disorders, detected in trans with a pathogenic variant",
        "pathogenic", Strength.MODERATE,
        "Phased genotypes for this patient showing a pathogenic variant on the other allele.",
    )


def pm4(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Protein length change from an in-frame indel or stop loss"
    cons = bundle.consequence
    if cons is None:
        return CriterionResult(
            code="PM4", name=name, direction="pathogenic", status="insufficient",
            strength=Strength.MODERATE,
            rationale="No molecular consequence was resolved.",
            what_would_change="A transcript-level consequence call.",
        )

    terms = set(cons.terms or [])
    length_changing = bool(terms & {"inframe_insertion", "inframe_deletion", "stop_lost"})
    if length_changing:
        return CriterionResult(
            code="PM4", name=name, direction="pathogenic", status="met",
            strength=Strength.MODERATE,
            rationale=f"Consequence {cons.most_severe} alters protein length.",
        )
    return CriterionResult(
        code="PM4", name=name, direction="pathogenic", status="not_met",
        strength=Strength.MODERATE,
        rationale=f"Consequence {cons.most_severe} does not change protein length.",
    )


def pm5(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PM5", "Novel missense at a residue where a different missense is pathogenic",
        "pathogenic", Strength.MODERATE,
        "A ClinVar search for other missense variants at this same residue with a "
        "pathogenic classification.",
    )


def pm6(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PM6", "Assumed de novo, without confirmation of parentage", "pathogenic",
        Strength.MODERATE,
        "Family history indicating the variant arose de novo.",
    )


# ==========================================================================
# Pathogenic — Supporting
# ==========================================================================


def pp1(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PP1", "Cosegregation with disease in multiple affected family members",
        "pathogenic", Strength.SUPPORTING,
        "Segregation data across affected and unaffected relatives; three or more "
        "informative meioses allow this criterion to be upgraded.",
    )


def pp2(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Missense in a gene with a low rate of benign missense variation"
    cons, constraint = bundle.consequence, bundle.constraint

    if cons is None or not cons.is_missense:
        return CriterionResult(
            code="PP2", name=name, direction="pathogenic", status="not_met",
            strength=Strength.SUPPORTING,
            rationale=(
                f"Consequence is {cons.most_severe if cons else 'unresolved'}; "
                "PP2 applies only to missense variants."
            ),
        )
    if constraint is None or constraint.mis_z is None:
        return CriterionResult(
            code="PP2", name=name, direction="pathogenic", status="insufficient",
            strength=Strength.SUPPORTING,
            rationale="No missense constraint (mis_z) available for this gene.",
            sources=_gene_provenance(bundle),
            what_would_change="Gene-level missense constraint, or a curated statement "
            "that missense is a common mechanism in this gene.",
        )

    if constraint.mis_z >= t.pp2_mis_z:
        return CriterionResult(
            code="PP2", name=name, direction="pathogenic", status="met",
            strength=Strength.SUPPORTING,
            rationale=(
                f"{bundle.gene_symbol} is missense-constrained (mis_z "
                f"{constraint.mis_z:.2f} >= {t.pp2_mis_z}). Constraint is a proxy for "
                "the criterion's requirement that missense be a common disease mechanism."
            ),
            sources=_gene_provenance(bundle),
        )
    return CriterionResult(
        code="PP2", name=name, direction="pathogenic", status="not_met",
        strength=Strength.SUPPORTING,
        rationale=(
            f"{bundle.gene_symbol} is not missense-constrained (mis_z "
            f"{constraint.mis_z:.2f} < {t.pp2_mis_z})."
        ),
        sources=_gene_provenance(bundle),
    )


def _revel(bundle: EvidenceBundle) -> float | None:
    for p in bundle.predictors:
        if p.name == "REVEL":
            return p.score
    return None


def _alphamissense(bundle: EvidenceBundle) -> float | None:
    for p in bundle.predictors:
        if p.name == "AlphaMissense":
            return p.score
    return None


def pp3(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Multiple lines of computational evidence support a deleterious effect"
    revel = _revel(bundle)
    am = _alphamissense(bundle)
    cons = bundle.consequence

    source = [
        Provenance(
            source="dbNSFP",
            label="REVEL, via MyVariant.info",
            url=f"https://myvariant.info/v1/variant/{bundle.query.hgvs_g}" if bundle.query.hgvs_g else None,
        )
    ]

    if cons is not None and not cons.is_missense:
        return CriterionResult(
            code="PP3", name=name, direction="pathogenic", status="not_met",
            strength=Strength.SUPPORTING,
            rationale=(
                f"The calibrated predictor used here (REVEL) applies to missense "
                f"variants; this variant is {cons.most_severe}."
            ),
            what_would_change="A predictor calibrated for this consequence class, "
            "such as SpliceAI for splice-region variants.",
        )

    if revel is None:
        return CriterionResult(
            code="PP3", name=name, direction="pathogenic", status="insufficient",
            strength=Strength.SUPPORTING,
            rationale="No REVEL score available for this variant.",
            what_would_change="A REVEL score, or another predictor with published "
            "strength calibration.",
        )

    # ClinGen SVI: use one calibrated predictor at its calibrated strength,
    # rather than counting how many uncalibrated tools agree.
    strength: Strength | None = None
    for cutoff, level in t.revel_strength_ladder:
        if revel >= cutoff:
            strength = level
            break

    if strength is None:
        return CriterionResult(
            code="PP3", name=name, direction="pathogenic", status="not_met",
            strength=Strength.SUPPORTING,
            rationale=(
                f"REVEL {revel:.3f} is below the supporting threshold "
                f"({t.revel_pp3_supporting}) from the Pejaver 2022 calibration."
            ),
            sources=source,
        )

    corroboration = ""
    if am is not None:
        agrees = am >= t.alphamissense_pathogenic
        corroboration = (
            f" AlphaMissense {am:.3f} {'agrees' if agrees else 'disagrees'}."
        )

    return CriterionResult(
        code="PP3", name=name, direction="pathogenic", status="met",
        strength=strength,
        rationale=(
            f"REVEL {revel:.3f} reaches the {strength.value} threshold from the "
            f"Pejaver 2022 calibration.{corroboration} Strength follows the single "
            "calibrated predictor, not a tally of agreeing tools."
        ),
        sources=source,
    )


def pp4(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "PP4", "Patient phenotype highly specific for a single-gene disease",
        "pathogenic", Strength.SUPPORTING,
        "The patient's phenotype, assessed against the gene's established "
        "phenotypic spectrum.",
    )


def pp5(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    cv = bundle.clinvar
    detail = "No ClinVar record."
    if cv and cv.actionable_assertions:
        path = [a for a in cv.actionable_assertions if a.significance.is_pathogenic_side]
        detail = (
            f"{len(path)} of {len(cv.actionable_assertions)} ClinVar assertions are "
            "pathogenic or likely pathogenic."
        )
    return CriterionResult(
        code="PP5", name="Reputable source reports pathogenic (retired)",
        direction="pathogenic", status="not_met", strength=Strength.SUPPORTING,
        rationale=(
            f"Retired by ClinGen SVI and deliberately not applied: using another "
            f"lab's conclusion as evidence double-counts the evidence they used. "
            f"{detail} Shown for context only."
        ),
        sources=_clinvar_provenance(bundle),
    )


# ==========================================================================
# Benign
# ==========================================================================


def ba1(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Allele frequency above 5% in a population database"
    freq, af = bundle.best_frequency(), bundle.max_af()

    if freq is None or af is None:
        return CriterionResult(
            code="BA1", name=name, direction="benign", status="insufficient",
            strength=Strength.STAND_ALONE,
            rationale="No population frequency data was returned.",
            sources=_gnomad_provenance(bundle),
            what_would_change="A gnomAD frequency for this variant.",
        )
    if freq.an is not None and freq.an < t.min_allele_number:
        return CriterionResult(
            code="BA1", name=name, direction="benign", status="insufficient",
            strength=Strength.STAND_ALONE,
            rationale=f"Allele number {freq.an:,} is too small for a frequency statement.",
            sources=_gnomad_provenance(bundle),
            what_would_change="Frequency from a better-powered cohort.",
        )

    if af >= t.ba1_af:
        return CriterionResult(
            code="BA1", name=name, direction="benign", status="met",
            strength=Strength.STAND_ALONE,
            rationale=(
                f"Frequency {_fmt_af(af)} in {freq.popmax_population or freq.source} "
                f"is at or above {_fmt_af(t.ba1_af)}. BA1 is stand-alone: on its own "
                "it classifies the variant benign."
            ),
            sources=_gnomad_provenance(bundle),
        )
    return CriterionResult(
        code="BA1", name=name, direction="benign", status="not_met",
        strength=Strength.STAND_ALONE,
        rationale=f"Highest population frequency {_fmt_af(af)} is below {_fmt_af(t.ba1_af)}.",
        sources=_gnomad_provenance(bundle),
    )


def bs1(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Allele frequency greater than expected for the disorder"
    freq, af = bundle.best_frequency(), bundle.max_af()

    if af is None:
        return CriterionResult(
            code="BS1", name=name, direction="benign", status="insufficient",
            strength=Strength.STRONG,
            rationale="No population frequency data was returned.",
            sources=_gnomad_provenance(bundle),
            what_would_change="A gnomAD frequency, plus the disorder's prevalence "
            "and penetrance to set a gene-specific threshold.",
        )

    if af >= t.bs1_af:
        return CriterionResult(
            code="BS1", name=name, direction="benign", status="met",
            strength=Strength.STRONG,
            rationale=(
                f"Frequency {_fmt_af(af)} exceeds the {_fmt_af(t.bs1_af)} working "
                f"threshold in {freq.popmax_population if freq else 'the population'}. "
                "This threshold is a default, not a gene-specific one: the criterion "
                "asks what is expected for the particular disorder, which depends on "
                "its prevalence, penetrance and inheritance."
            ),
            sources=_gnomad_provenance(bundle),
            what_would_change="A disorder-specific threshold from ClinGen's "
            "variant-curation expert panel for this gene.",
        )
    return CriterionResult(
        code="BS1", name=name, direction="benign", status="not_met",
        strength=Strength.STRONG,
        rationale=f"Frequency {_fmt_af(af)} is below the {_fmt_af(t.bs1_af)} working threshold.",
        sources=_gnomad_provenance(bundle),
    )


def bs2(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Observed in a healthy individual in the expected genotype state"
    freqs = [f for f in bundle.frequencies if f.homozygote_count is not None]
    if not freqs:
        return CriterionResult(
            code="BS2", name=name, direction="benign", status="insufficient",
            strength=Strength.STRONG,
            rationale="No homozygote counts were returned.",
            sources=_gnomad_provenance(bundle),
            what_would_change="Homozygote counts, and the disorder's inheritance mode.",
        )

    best = max(freqs, key=lambda f: f.homozygote_count or 0)
    hom = best.homozygote_count or 0

    if hom >= t.bs2_min_homozygotes:
        return CriterionResult(
            code="BS2", name=name, direction="benign", status="met",
            strength=Strength.STRONG,
            rationale=(
                f"{hom:,} homozygotes in {best.source}. For a recessive disorder with "
                "full penetrance this argues benign. Whether it applies depends on the "
                "inheritance mode, which is not encoded here, and gnomAD is not a "
                "fully healthy cohort."
            ),
            sources=_gnomad_provenance(bundle),
            what_would_change="The disorder's inheritance mode and age of onset, to "
            "confirm homozygotes in gnomAD would be expected to be affected.",
        )
    return CriterionResult(
        code="BS2", name=name, direction="benign", status="not_met",
        strength=Strength.STRONG,
        rationale=f"Only {hom} homozygotes in {best.source}.",
        sources=_gnomad_provenance(bundle),
    )


def bs3(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "BS3", "Well-established functional study shows no damaging effect", "benign",
        Strength.STRONG,
        "A validated functional assay reporting normal function.",
    )


def bs4(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "BS4", "Lack of segregation in affected family members", "benign",
        Strength.STRONG,
        "Segregation data showing affected relatives who do not carry the variant.",
    )


def bp1(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "BP1", "Missense in a gene where only truncating variants cause disease",
        "benign", Strength.SUPPORTING,
        "A curated statement of the gene's disease mechanism. Constraint metrics "
        "cannot establish that only truncating variants are causal.",
    )


def bp2(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "BP2", "Observed in trans with a pathogenic variant in a dominant gene",
        "benign", Strength.SUPPORTING,
        "Phased genotypes for this patient.",
    )


def bp3(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "BP3", "In-frame indel in a repetitive region of no known function", "benign",
        Strength.SUPPORTING,
        "Repeat annotation for the affected region and confirmation it has no known "
        "functional role.",
    )


def bp4(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Computational evidence suggests no impact"
    revel = _revel(bundle)
    cons = bundle.consequence

    if cons is not None and not cons.is_missense:
        return CriterionResult(
            code="BP4", name=name, direction="benign", status="not_met",
            strength=Strength.SUPPORTING,
            rationale=f"REVEL applies to missense variants; this is {cons.most_severe}.",
        )
    if revel is None:
        return CriterionResult(
            code="BP4", name=name, direction="benign", status="insufficient",
            strength=Strength.SUPPORTING,
            rationale="No REVEL score available.",
            what_would_change="A REVEL score, or another calibrated predictor.",
        )

    strength: Strength | None = None
    if revel <= t.revel_bp4_strong:
        strength = Strength.STRONG
    elif revel <= t.revel_bp4_moderate:
        strength = Strength.MODERATE
    elif revel <= t.revel_bp4_supporting:
        strength = Strength.SUPPORTING

    if strength is None:
        return CriterionResult(
            code="BP4", name=name, direction="benign", status="not_met",
            strength=Strength.SUPPORTING,
            rationale=(
                f"REVEL {revel:.3f} is above the benign supporting threshold "
                f"({t.revel_bp4_supporting})."
            ),
        )
    return CriterionResult(
        code="BP4", name=name, direction="benign", status="met", strength=strength,
        rationale=(
            f"REVEL {revel:.3f} reaches the {strength.value} benign threshold from the "
            "Pejaver 2022 calibration."
        ),
    )


def bp5(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    return _manual(
        "BP5", "Found in a case with an alternate molecular basis for disease", "benign",
        Strength.SUPPORTING,
        "The patient's other findings, specifically a causal variant elsewhere.",
    )


def bp6(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    cv = bundle.clinvar
    detail = "No ClinVar record."
    if cv and cv.actionable_assertions:
        benign = [a for a in cv.actionable_assertions if a.significance.is_benign_side]
        detail = (
            f"{len(benign)} of {len(cv.actionable_assertions)} ClinVar assertions are "
            "benign or likely benign."
        )
    return CriterionResult(
        code="BP6", name="Reputable source reports benign (retired)", direction="benign",
        status="not_met", strength=Strength.SUPPORTING,
        rationale=(
            f"Retired by ClinGen SVI and deliberately not applied. {detail} "
            "Shown for context only."
        ),
        sources=_clinvar_provenance(bundle),
    )


def bp7(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> CriterionResult:
    name = "Synonymous variant with no predicted splice impact"
    cons = bundle.consequence
    if cons is None:
        return CriterionResult(
            code="BP7", name=name, direction="benign", status="insufficient",
            strength=Strength.SUPPORTING,
            rationale="No molecular consequence was resolved.",
            what_would_change="A transcript-level consequence call.",
        )
    if not cons.is_synonymous:
        return CriterionResult(
            code="BP7", name=name, direction="benign", status="not_met",
            strength=Strength.SUPPORTING,
            rationale=f"BP7 applies to synonymous variants; this is {cons.most_severe}.",
        )

    splice = bundle.splice_ai_max
    if splice is None:
        return CriterionResult(
            code="BP7", name=name, direction="benign", status="insufficient",
            strength=Strength.SUPPORTING,
            rationale="Synonymous, but no SpliceAI score was available to rule out a splice effect.",
            what_would_change="A SpliceAI delta score, and nucleotide conservation at the position.",
        )
    if splice < t.spliceai_no_effect:
        return CriterionResult(
            code="BP7", name=name, direction="benign", status="met",
            strength=Strength.SUPPORTING,
            rationale=(
                f"Synonymous with SpliceAI max delta {splice:.3f}, below "
                f"{t.spliceai_no_effect}. The criterion also asks that the nucleotide "
                "not be conserved, which is not checked here."
            ),
        )
    return CriterionResult(
        code="BP7", name=name, direction="benign", status="not_met",
        strength=Strength.SUPPORTING,
        rationale=f"Synonymous, but SpliceAI max delta {splice:.3f} suggests a possible splice effect.",
    )


# ==========================================================================
# Registry — order follows the ACMG/AMP tables
# ==========================================================================

ALL_CRITERIA: tuple[CriterionFn, ...] = (
    pvs1,
    ps1, ps2, ps3, ps4,
    pm1, pm2, pm3, pm4, pm5, pm6,
    pp1, pp2, pp3, pp4, pp5,
    ba1,
    bs1, bs2, bs3, bs4,
    bp1, bp2, bp3, bp4, bp5, bp6, bp7,
)


def evaluate_all(bundle: EvidenceBundle, t: Thresholds = DEFAULT) -> list[CriterionResult]:
    """Evaluate every criterion. Always returns all 28."""
    return [fn(bundle, t) for fn in ALL_CRITERIA]
