"""The suite's offline claim, asserted rather than assumed.

The README and the CI workflow both say this suite cannot reach the internet.
That is only true if the guard in conftest.py works, and the guard was loosened
to let asyncio build the self-pipe an event loop needs -- a change that could
quietly have opened the door it exists to close. These tests hold it shut.

They run only under --no-network, which is the mode the claim is about.
"""

import socket

import httpx
import pytest


def _offline(request) -> bool:
    return bool(request.config.getoption("--no-network"))


@pytest.fixture(autouse=True)
def _only_offline(request):
    if not _offline(request):
        pytest.skip("this file asserts what --no-network guarantees")


@pytest.mark.parametrize("family", [socket.AF_INET, socket.AF_INET6])
def test_an_ip_socket_cannot_be_opened(family):
    with pytest.raises(RuntimeError, match="offline"):
        socket.socket(family, socket.SOCK_STREAM)


def test_create_connection_is_refused():
    with pytest.raises(RuntimeError, match="offline"):
        socket.create_connection(("example.com", 80), timeout=1)


def test_an_actual_http_request_fails():
    """The end the guard exists for: no test can reach a public API, whatever
    route it takes to try."""
    with pytest.raises(Exception):
        httpx.get("https://myvariant.info/v1/query?q=rs334", timeout=5)


def test_a_local_socketpair_is_still_allowed():
    """asyncio needs one for its self-pipe. Blocking it meant no async code
    could be tested offline at all, which is how the panel's concurrency and
    failure-isolation paths went untested."""
    a, b = socket.socketpair()
    try:
        assert a.fileno() >= 0
    finally:
        a.close()
        b.close()


async def test_an_event_loop_runs_under_the_guard():
    import asyncio

    await asyncio.sleep(0)
