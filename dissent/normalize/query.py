"""Parse free-form variant input into a structured VariantQuery.

People write variants in whatever notation their field uses. All of these refer
to the same BRAF variant and all must work:

    rs113488022
    BRAF V600E
    BRAF p.Val600Glu
    chr7:140753336A>T
    7-140753336-A-T
    NM_004333.6:c.1799T>A

Parsing is deliberately strict about what it claims to recognise: anything it
cannot confidently classify becomes `kind="unknown"` and is handed to the
source layer as a free-text search rather than being guessed at.
"""

from __future__ import annotations

import re

from dissent.models import Assembly, VariantQuery

# Three-letter to one-letter amino acid codes, plus Ter/* for stops.
AA3_TO_1 = {
    "ala": "A", "arg": "R", "asn": "N", "asp": "D", "cys": "C",
    "gln": "Q", "glu": "E", "gly": "G", "his": "H", "ile": "I",
    "leu": "L", "lys": "K", "met": "M", "phe": "F", "pro": "P",
    "ser": "S", "thr": "T", "trp": "W", "tyr": "Y", "val": "V",
    "ter": "*", "sec": "U", "asx": "B", "glx": "Z", "xaa": "X",
}

_RSID = re.compile(r"^rs\d+$", re.I)

# chr7:140753336A>T | chr7:g.140753336A>T | 7:140753336A>T
_GCOORD_HGVS = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|[XYMT]{1,2})"
    r"\s*[:\s]\s*(?:g\.)?(?P<pos>\d+)\s*"
    r"(?P<ref>[ACGT]+)\s*>\s*(?P<alt>[ACGT]+)$",
    re.I,
)

# 7-140753336-A-T  |  7:140753336:A:T
_GCOORD_DELIM = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|[XYMT]{1,2})"
    r"[-:](?P<pos>\d+)[-:](?P<ref>[ACGT]+)[-:](?P<alt>[ACGT]+)$",
    re.I,
)

# NM_004333.6:c.1799T>A  | ENST00000288602.11:c.1799T>A
_HGVS_C = re.compile(
    r"^(?P<tx>(?:NM_|NR_|ENST)[0-9.]+)\s*:\s*(?P<cdot>c\.\S+)$", re.I
)

# BRAF V600E | BRAF p.V600E | BRAF p.Val600Glu | BRAF:p.Val600Glu
_GENE_PROTEIN = re.compile(
    r"^(?P<gene>[A-Z][A-Z0-9\-]{1,19})"
    r"[\s:]+(?:p\.)?"
    r"(?P<ref>[A-Z][a-z]{2}|[A-Z])"
    r"(?P<pos>\d+)"
    r"(?P<alt>[A-Z][a-z]{2}|[A-Z*]|Ter)$",
    re.I,
)


def _norm_chrom(c: str) -> str:
    c = c.upper().removeprefix("CHR")
    return "MT" if c in ("M", "MT") else c


def _aa_to_one(aa: str) -> str:
    if len(aa) == 1:
        return aa.upper()
    return AA3_TO_1.get(aa.lower(), aa.upper())


def parse_variant(raw: str, assembly: Assembly = "hg38") -> VariantQuery:
    """Parse user input. Never raises; falls back to kind='unknown'."""
    text = (raw or "").strip()
    if not text:
        return VariantQuery(raw=raw or "", kind="unknown", assembly=assembly)

    collapsed = re.sub(r"\s+", " ", text)

    if _RSID.match(collapsed):
        return VariantQuery(
            raw=text, kind="rsid", assembly=assembly, rsid=collapsed.lower()
        )

    for pattern in (_GCOORD_DELIM, _GCOORD_HGVS):
        m = pattern.match(collapsed)
        if m:
            chrom = _norm_chrom(m.group("chrom"))
            pos = int(m.group("pos"))
            ref = m.group("ref").upper()
            alt = m.group("alt").upper()
            return VariantQuery(
                raw=text,
                kind="gcoord",
                assembly=assembly,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                hgvs_g=f"chr{chrom}:g.{pos}{ref}>{alt}",
            )

    m = _HGVS_C.match(collapsed)
    if m:
        return VariantQuery(
            raw=text,
            kind="hgvs",
            assembly=assembly,
            transcript=m.group("tx"),
            hgvs_c=m.group("cdot"),
        )

    m = _GENE_PROTEIN.match(collapsed)
    if m:
        ref1 = _aa_to_one(m.group("ref"))
        alt1 = _aa_to_one(m.group("alt"))
        return VariantQuery(
            raw=text,
            kind="gene_protein",
            assembly=assembly,
            gene=m.group("gene").upper(),
            protein_change=f"p.{ref1}{m.group('pos')}{alt1}",
        )

    return VariantQuery(raw=text, kind="unknown", assembly=assembly)
