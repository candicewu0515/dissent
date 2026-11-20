"""The parser must accept every notation a user might actually type."""

import pytest

from dissent.normalize.query import parse_variant
from dissent.sources.myvariant import three_letter_protein


@pytest.mark.parametrize(
    "raw,kind",
    [
        ("rs113488022", "rsid"),
        ("RS113488022", "rsid"),
        ("rs1801133", "rsid"),
        ("chr7:140753336A>T", "gcoord"),
        ("chr7:g.140753336A>T", "gcoord"),
        ("7:140753336A>T", "gcoord"),
        ("7-140753336-A-T", "gcoord"),
        ("7:140753336:A:T", "gcoord"),
        ("BRAF V600E", "gene_protein"),
        ("BRAF p.Val600Glu", "gene_protein"),
        ("braf p.V600E", "gene_protein"),
        ("BRAF:p.Val600Glu", "gene_protein"),
        ("NM_004333.6:c.1799T>A", "hgvs"),
        ("ENST00000288602.11:c.1799T>A", "hgvs"),
    ],
)
def test_recognised_notations(raw, kind):
    assert parse_variant(raw).kind == kind


@pytest.mark.parametrize("raw", ["", "   ", "total garbage", "not a variant at all", "???"])
def test_unrecognised_input_is_not_guessed_at(raw):
    """Better to report `unknown` than to invent a parse."""
    assert parse_variant(raw).kind == "unknown"


def test_genomic_coordinates_normalise_identically():
    forms = [
        "chr7:140753336A>T",
        "chr7:g.140753336A>T",
        "7-140753336-A-T",
        "7:140753336:A:T",
    ]
    parsed = [parse_variant(f) for f in forms]
    assert {(q.chrom, q.pos, q.ref, q.alt) for q in parsed} == {("7", 140753336, "A", "T")}
    assert {q.gnomad_variant_id for q in parsed} == {"7-140753336-A-T"}


def test_mitochondrial_contig_normalises_to_mt():
    assert parse_variant("chrM:500T>C").chrom == "MT"
    assert parse_variant("chrMT:500T>C").chrom == "MT"


def test_three_letter_and_one_letter_protein_are_equivalent():
    a = parse_variant("BRAF V600E")
    b = parse_variant("BRAF p.Val600Glu")
    assert a.protein_change == b.protein_change == "p.V600E"


def test_hgvs_transcript_and_cdot_are_kept():
    q = parse_variant("NM_004333.6:c.1799T>A")
    assert q.transcript == "NM_004333.6"
    assert q.hgvs_c == "c.1799T>A"


@pytest.mark.parametrize(
    "one,three",
    [("p.C61G", "p.Cys61Gly"), ("p.V600E", "p.Val600Glu"), ("p.Q100*", "p.Gln100Ter")],
)
def test_protein_change_rewrites_to_clinvar_notation(one, three):
    assert three_letter_protein(one) == three


def test_protein_rewrite_refuses_unparseable_input():
    assert three_letter_protein("garbage") is None
    assert three_letter_protein(None) is None
