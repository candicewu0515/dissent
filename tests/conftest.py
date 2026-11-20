"""Fixture loading, and the switch that proves the suite is offline.

Every external response in this repository is a recorded fixture. `--no-network`
replaces the socket constructor so any attempt to open one fails loudly: it
turns "these tests do not need the network" from a claim into something CI
enforces.
"""

from __future__ import annotations

import json
import socket
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


def pytest_addoption(parser):
    parser.addoption(
        "--no-network",
        action="store_true",
        default=False,
        help="Fail any test that tries to open a socket.",
    )


@pytest.fixture(autouse=True)
def _forbid_network(request, monkeypatch):
    if not request.config.getoption("--no-network"):
        return

    real_socket = socket.socket

    def refuse(*args, **kwargs):
        raise RuntimeError(
            "This suite runs offline. Something tried to open a network "
            "connection; it should be reading a fixture instead."
        )

    def guarded_socket(family=socket.AF_INET, *args, **kwargs):
        """Block the internet, not every socket.

        Replacing the constructor outright also blocked AF_UNIX socketpairs,
        which is how asyncio builds the self-pipe an event loop needs -- so no
        async code could be tested offline at all, and the concurrency and
        failure-isolation paths went untested as a result.

        Blocking by address family is the stronger guarantee anyway: it names
        what is forbidden (an IP socket) rather than what happens to be the
        usual way of getting one.
        """
        if family in (socket.AF_INET, socket.AF_INET6):
            refuse()
        return real_socket(family, *args, **kwargs)

    monkeypatch.setattr(socket, "socket", guarded_socket)
    monkeypatch.setattr(socket, "create_connection", refuse)
