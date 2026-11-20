"""Numeric thresholds used by the criteria, in one place with their sources.

Every number a criterion depends on lives here rather than inline, for two
reasons: a reviewer can audit the assumptions in one screen, and a lab with
gene-specific thresholds can override them without touching the rule logic.

References
----------
Richards S, et al. Standards and guidelines for the interpretation of sequence
variants. Genet Med. 2015;17(5):405-424.  -- the original ACMG/AMP framework.

Pejaver V, et al. Calibration of computational tools for missense variant
pathogenicity classification. Am J Hum Genet. 2022;109(12):2163-2177.
-- the calibrated REVEL thresholds used for PP3/BP4 below.

ClinGen Sequence Variant Interpretation working group recommendations:
  * PM2 should be applied at Supporting strength, not Moderate.
  * PP5 and BP6 are retired and should not contribute to a classification.
  * PP3/BP4 should rest on one calibrated predictor at a calibrated strength,
    rather than on a count of tools that happen to agree.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from dissent.models import Strength


@dataclass(frozen=True)
class Thresholds:
    # ---- allele frequency -------------------------------------------------
    # BA1: Richards 2015 gives 5% as the stand-alone benign cutoff.
    ba1_af: float = 0.05
    # BS1: "greater than expected for the disorder" is disease-specific. 1% is
    # a conventional default for a dominant condition with moderate penetrance
    # and is flagged as an assumption in the rationale.
    bs1_af: float = 0.01
    # PM2: ClinGen SVI treats absence or near-absence as Supporting only.
    pm2_af: float = 0.0001
    # A frequency estimate from a small cohort is not informative either way.
    min_allele_number: int = 2000

    # ---- gene constraint --------------------------------------------------
    # pLI / LOEUF stand in for "loss of function is a known disease mechanism".
    # This is a proxy, not the statement PVS1 actually requires, and the
    # criterion says so in its rationale.
    pvs1_pli: float = 0.9
    pvs1_loeuf: float = 0.35
    # PP2 requires a low rate of benign missense variation in the gene.
    pp2_mis_z: float = 3.09

    # ---- REVEL, calibrated (Pejaver 2022) --------------------------------
    revel_pp3_strong: float = 0.932
    revel_pp3_moderate: float = 0.773
    revel_pp3_supporting: float = 0.644
    revel_bp4_supporting: float = 0.290
    revel_bp4_moderate: float = 0.183
    revel_bp4_strong: float = 0.016

    # ---- AlphaMissense, used only as corroboration -----------------------
    alphamissense_pathogenic: float = 0.564
    alphamissense_benign: float = 0.340

    # ---- splicing ---------------------------------------------------------
    # SpliceAI delta below this is treated as no predicted splice effect.
    spliceai_no_effect: float = 0.10
    spliceai_effect: float = 0.50

    # ---- BS2 --------------------------------------------------------------
    # One homozygote can be a genotyping artefact; two is harder to dismiss.
    bs2_min_homozygotes: int = 2

    # ---- predictor consensus (reported, never scored) --------------------
    # Used for the conflict panel's "predictors split" detector, not for PP3.
    predictor_split_min_tools: int = 6
    predictor_split_min_minority_fraction: float = 0.25

    # ---- staleness --------------------------------------------------------
    stale_assertion_years: int = 5

    revel_strength_ladder: tuple[tuple[float, Strength], ...] = field(
        default=(
            (0.932, Strength.STRONG),
            (0.773, Strength.MODERATE),
            (0.644, Strength.SUPPORTING),
        )
    )


DEFAULT = Thresholds()
