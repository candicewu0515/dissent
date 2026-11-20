"""Fixture loading. No test in this suite touches the network."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dissent.models import EvidenceBundle

FIXTURE_DIR = Path(__file__).parent / "fixtures"


def load_bundle(name: str) -> EvidenceBundle:
    path = FIXTURE_DIR / f"{name}.json"
    return EvidenceBundle.model_validate_json(path.read_text())


def available_bundles() -> list[str]:
    return sorted(p.stem for p in FIXTURE_DIR.glob("*.json"))


@pytest.fixture
def braf() -> EvidenceBundle:
    """BRAF V600E: somatic/germline mix, discordant submitters, split predictors."""
    return load_bundle("braf_v600e")


@pytest.fixture
def mthfr() -> EvidenceBundle:
    """MTHFR A222V: common polymorphism that predictors call damaging."""
    return load_bundle("mthfr_a222v")


@pytest.fixture
def brca1() -> EvidenceBundle:
    """BRCA1 C61G: well-curated pathogenic missense."""
    return load_bundle("brca1_c61g")


@pytest.fixture
def tp53() -> EvidenceBundle:
    """TP53 R175H: the variant whose snpEff record includes a PDB annotation."""
    return load_bundle("tp53_r175h")
