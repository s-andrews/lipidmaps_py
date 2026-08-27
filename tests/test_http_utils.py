"""Tests for the rate-limit-aware HTTP retry helper."""

import requests

from lipidmaps import http_utils


class FakeResponse:
    def __init__(self, status_code, headers=None):
        self.status_code = status_code
        self.headers = headers or {}


def test_injects_default_user_agent(monkeypatch):
    seen = {}

    def fake_post(url, headers=None, **kwargs):
        seen["headers"] = headers
        return FakeResponse(200)

    monkeypatch.setattr(requests, "post", fake_post)
    http_utils.post_with_retry("http://example/api", json={"a": 1})
    assert seen["headers"]["User-Agent"] == "lipidmaps_py"


def test_caller_headers_merge_over_defaults(monkeypatch):
    seen = {}

    def fake_get(url, headers=None, **kwargs):
        seen["headers"] = headers
        return FakeResponse(200)

    monkeypatch.setattr(requests, "get", fake_get)
    http_utils.get_with_retry("http://example/api", headers={"X-Test": "1"})
    assert seen["headers"]["User-Agent"] == "lipidmaps_py"
    assert seen["headers"]["X-Test"] == "1"


def test_retries_on_429_then_succeeds(monkeypatch):
    calls = []
    sleeps = []
    monkeypatch.setattr(http_utils.time, "sleep", lambda s: sleeps.append(s))

    def fake_post(url, headers=None, **kwargs):
        calls.append(1)
        if len(calls) < 3:
            return FakeResponse(429, {"Retry-After": "2"})
        return FakeResponse(200)

    monkeypatch.setattr(requests, "post", fake_post)
    resp = http_utils.post_with_retry("http://example/api", json={})
    assert resp.status_code == 200
    assert len(calls) == 3
    # Retry-After honored (2s), not the exponential default.
    assert sleeps == [2.0, 2.0]


def test_gives_up_after_max_retries(monkeypatch):
    monkeypatch.setattr(http_utils.time, "sleep", lambda s: None)

    def fake_post(url, headers=None, **kwargs):
        return FakeResponse(429, {"Retry-After": "1"})

    monkeypatch.setattr(requests, "post", fake_post)
    resp = http_utils.post_with_retry("http://example/api", json={}, max_retries=2)
    assert resp.status_code == 429


def test_non_retryable_status_returned_immediately(monkeypatch):
    calls = []

    def fake_post(url, headers=None, **kwargs):
        calls.append(1)
        return FakeResponse(400)

    monkeypatch.setattr(requests, "post", fake_post)
    resp = http_utils.post_with_retry("http://example/api", json={})
    assert resp.status_code == 400
    assert len(calls) == 1
