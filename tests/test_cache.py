"""The response cache.

A cache is an optimisation, so the property that matters most is not that it
stores things -- it is that it can never fail a request. Every test below that
breaks the storage asserts the caller still gets an answer.
"""

import sqlite3

import pytest

from dissent.cache.store import ResponseCache


@pytest.fixture
def cache(tmp_path) -> ResponseCache:
    return ResponseCache(path=tmp_path / "c.sqlite")


# ---- the ordinary path ----------------------------------------------------


def test_round_trip(cache):
    cache.set("k", "clinvar", {"a": 1})
    assert cache.get("k") == {"a": 1}


def test_a_miss_is_none_not_an_error(cache):
    assert cache.get("never-stored") is None


def test_keys_differ_by_source_and_request():
    a = ResponseCache.make_key("clinvar", {"q": "BRAF"})
    b = ResponseCache.make_key("gnomad", {"q": "BRAF"})
    c = ResponseCache.make_key("clinvar", {"q": "TP53"})
    assert a != b and a != c


def test_keys_are_stable_across_dict_ordering():
    """Two requests that differ only in key order are the same request. If the
    digest disagreed, every lookup would miss and the rate limit would be hit
    on a tool that looks cached."""
    a = ResponseCache.make_key("s", {"x": 1, "y": 2})
    b = ResponseCache.make_key("s", {"y": 2, "x": 1})
    assert a == b


def test_expired_entries_are_not_returned(tmp_path):
    c = ResponseCache(path=tmp_path / "c.sqlite", ttl=0)
    c.set("k", "s", {"a": 1})
    assert c.get("k") is None


def test_a_fresh_entry_within_the_ttl_is_returned(tmp_path):
    c = ResponseCache(path=tmp_path / "c.sqlite", ttl=3600)
    c.set("k", "s", {"a": 1})
    assert c.get("k") == {"a": 1}


def test_clear_by_source_leaves_the_others(cache):
    cache.set("a", "clinvar", {"n": 1})
    cache.set("b", "gnomad", {"n": 2})
    assert cache.clear("clinvar") == 1
    assert cache.get("a") is None
    assert cache.get("b") == {"n": 2}


def test_stats_counts_per_source(cache):
    cache.set("a", "clinvar", {})
    cache.set("b", "clinvar", {})
    cache.set("c", "gnomad", {})
    assert cache.stats() == {"clinvar": 2, "gnomad": 1}


def test_corrupt_payload_reads_as_a_miss(cache):
    """A row that will not parse is a miss, never an exception thrown at the
    caller mid-interpretation."""
    cache.set("k", "s", {"a": 1})
    with sqlite3.connect(cache.path) as conn:
        conn.execute("UPDATE responses SET payload = ? WHERE key = ?", ("{not json", "k"))
    assert cache.get("k") is None


# ---- the failure path -----------------------------------------------------


def test_an_unwritable_location_disables_rather_than_raises():
    """A read-only filesystem is the normal case in a container. Raising here
    would take the whole server down at startup."""
    c = ResponseCache(path="/proc/no/such/dir/c.sqlite")
    assert c.enabled is False
    assert c.disabled_reason


def test_a_disabled_cache_still_answers_every_call():
    c = ResponseCache(path="/proc/no/such/dir/c.sqlite")
    assert c.get("k") is None
    c.set("k", "s", {"a": 1})  # must not raise
    assert c.get("k") is None
    assert c.clear() == 0
    assert c.stats()["enabled"] is False


def test_stats_says_why_it_is_off():
    c = ResponseCache(path="/proc/no/such/dir/c.sqlite")
    stats = c.stats()
    assert stats["enabled"] is False
    assert stats["reason"]


def test_a_sqlite_error_mid_run_turns_the_cache_off_without_raising(cache, monkeypatch):
    """Disk full, database locked, file deleted underneath us. The caller gets
    a miss and the server keeps serving."""
    cache.set("k", "s", {"a": 1})

    def boom(*_args, **_kwargs):
        raise sqlite3.OperationalError("database is locked")

    monkeypatch.setattr(sqlite3, "connect", boom)
    assert cache.get("k") is None
    assert cache.enabled is False


def test_an_unserialisable_payload_does_not_disable_the_cache(cache):
    """One bad response is a fault in that response, not in the cache. The
    cache must stay on for everything after it."""
    cache.set("bad", "s", {"f": object()})
    assert cache.enabled is True
    cache.set("good", "s", {"a": 1})
    assert cache.get("good") == {"a": 1}


def test_the_path_comes_from_the_environment(tmp_path, monkeypatch):
    """How a deployment points the cache at a writable volume."""
    target = tmp_path / "vol" / "c.sqlite"
    monkeypatch.setenv("DISSENT_CACHE_PATH", str(target))
    import importlib

    from dissent.cache import store

    importlib.reload(store)
    try:
        c = store.ResponseCache()
        assert c.path == target
        assert c.enabled is True
    finally:
        monkeypatch.delenv("DISSENT_CACHE_PATH", raising=False)
        importlib.reload(store)
