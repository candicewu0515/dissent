"""Transcript selection from snpEff annotations.

The regression guarded here is real: snpEff emits structural annotations whose
feature_id embeds a transcript accession late in the string, and an unanchored
match selected a PDB interaction record as the consequence for TP53 p.Arg175His.
"""

from dissent.models import Consequence
from dissent.normalize.consequence import (
    _bare_accession,
    disagreement,
    from_snpeff,
    transcript_from_preferred_name,
)


def test_preferred_name_yields_the_clinvar_transcript():
    name = "NM_007294.4(BRCA1):c.181T>G (p.Cys61Gly)"
    assert transcript_from_preferred_name(name) == "NM_007294"


def test_accession_match_is_anchored_at_the_start():
    """A PDB-style feature_id must not be read as its embedded transcript."""
    assert _bare_accession("1GZH:A_175-A_237:NM_000546.5") is None
    assert _bare_accession("NM_000546.5") == "NM_000546"


def test_structural_annotations_are_not_chosen_as_the_consequence():
    doc = {
        "snpeff": {
            "ann": [
                {
                    "effect": "structural_interaction_variant",
                    "feature_id": "1GZH:A_175-A_237:NM_000546.5",
                    "feature_type": "interaction",
                    "putative_impact": "HIGH",
                    "genename": "TP53",
                },
                {
                    "effect": "missense_variant",
                    "feature_id": "NM_000546.5",
                    "feature_type": "transcript",
                    "transcript_biotype": "protein_coding",
                    "putative_impact": "MODERATE",
                    "genename": "TP53",
                    "hgvs_c": "c.524G>A",
                    "hgvs_p": "p.Arg175His",
                    "protein": {"position": "175"},
                },
            ]
        }
    }
    cons = from_snpeff(doc, preferred_transcript="NM_000546")
    assert cons is not None
    assert cons.most_severe == "missense_variant"
    assert cons.hgvs_p == "p.Arg175His"


def test_clinvar_transcript_wins_over_a_more_severe_one():
    """Residue numbering must follow the transcript ClinVar used."""
    doc = {
        "snpeff": {
            "ann": [
                {
                    "effect": "missense_variant",
                    "feature_id": "NM_007297.3",
                    "feature_type": "transcript",
                    "transcript_biotype": "protein_coding",
                    "putative_impact": "MODERATE",
                    "hgvs_p": "p.Cys14Gly",
                    "protein": {"position": "14", "length": "1816"},
                },
                {
                    "effect": "missense_variant",
                    "feature_id": "NM_007294.3",
                    "feature_type": "transcript",
                    "transcript_biotype": "protein_coding",
                    "putative_impact": "MODERATE",
                    "hgvs_p": "p.Cys61Gly",
                    "protein": {"position": "61", "length": "1863"},
                },
            ]
        }
    }
    cons = from_snpeff(doc, preferred_transcript="NM_007294")
    assert cons is not None and cons.hgvs_p == "p.Cys61Gly"
    assert cons.canonical is True


def test_compound_effects_are_split():
    doc = {
        "snpeff": {
            "ann": [
                {
                    "effect": "splice_region_variant&intron_variant",
                    "feature_id": "NM_1.1",
                    "feature_type": "transcript",
                    "transcript_biotype": "protein_coding",
                    "putative_impact": "LOW",
                }
            ]
        }
    }
    cons = from_snpeff(doc)
    assert cons is not None
    assert cons.terms == ["splice_region_variant", "intron_variant"]


def test_no_usable_annotation_returns_none():
    assert from_snpeff({}) is None
    assert from_snpeff({"snpeff": {"ann": []}}) is None
    assert (
        from_snpeff(
            {"snpeff": {"ann": [{"effect": "sequence_feature", "feature_type": "other"}]}}
        )
        is None
    )


def test_disagreement_reports_only_material_differences():
    missense = Consequence(most_severe="missense_variant", terms=["missense_variant"])
    other_tx = Consequence(
        most_severe="missense_variant", terms=["missense_variant"], transcript_id="NM_9.9"
    )
    stop = Consequence(most_severe="stop_gained", terms=["stop_gained"])

    # Same effect class on a different transcript is not a disagreement.
    assert disagreement(missense, other_tx) is None
    # A different effect class changes which criteria apply, so it is.
    assert disagreement(missense, stop) == "missense vs loss of function"
    assert disagreement(missense, None) is None
