"""Shared HTTP helpers: identify this client and tolerate dev-API rate limits.

The dev.lipidmaps.org API applies application-level rate limiting (Flask-Limiter
style): once a per-window budget is exhausted it returns HTTP 429 with a
``Retry-After`` header, fronted by — but not enforced by — Cloudflare. A single
pipeline run stays under the budget, but batched or repeated runs can trip it.

These helpers send the package's User-Agent and retry retryable responses (429
plus transient 5xx), honoring ``Retry-After`` so callers don't fail on transient
throttling. They call ``requests.post``/``requests.get`` by attribute so existing
``monkeypatch.setattr(requests, "post", ...)`` test doubles still intercept.
"""

from __future__ import annotations

import logging
import time

import requests

from .config import DEFAULT_HEADERS

logger = logging.getLogger(__name__)

# Rate limit + transient upstream errors are worth retrying; 4xx (other than
# 429) are the caller's problem and returned as-is.
_RETRY_STATUS = {429, 502, 503, 504}
_MAX_RETRIES = 4
# Cap so a hostile/misconfigured Retry-After can't stall a run indefinitely.
_MAX_BACKOFF_SECONDS = 60.0
_DEFAULT_BACKOFF_SECONDS = 5.0


def _retry_after_seconds(response: requests.Response, attempt: int) -> float:
    """Seconds to wait before the next attempt, honoring ``Retry-After``."""
    raw = response.headers.get("Retry-After")
    if raw:
        try:
            return min(float(raw), _MAX_BACKOFF_SECONDS)
        except ValueError:
            pass
    # No/invalid header: exponential backoff, capped.
    return min(_DEFAULT_BACKOFF_SECONDS * (2 ** attempt), _MAX_BACKOFF_SECONDS)


def request_with_retry(
    method: str,
    url: str,
    *,
    max_retries: int = _MAX_RETRIES,
    headers: dict | None = None,
    **kwargs,
) -> requests.Response:
    """Issue an HTTP request, retrying rate-limit / transient-error responses.

    ``headers`` are merged over the package defaults (so the User-Agent is always
    sent). All other keyword arguments are forwarded to ``requests``.
    """
    merged_headers = {**DEFAULT_HEADERS, **(headers or {})}
    caller = getattr(requests, method.lower())
    for attempt in range(max_retries + 1):
        response = caller(url, headers=merged_headers, **kwargs)
        if response.status_code not in _RETRY_STATUS or attempt == max_retries:
            return response
        wait = _retry_after_seconds(response, attempt)
        logger.warning(
            "HTTP %s from %s; retry %d/%d after %.1fs",
            response.status_code,
            url,
            attempt + 1,
            max_retries,
            wait,
        )
        time.sleep(wait)
    return response  # pragma: no cover - loop always returns above


def post_with_retry(url: str, **kwargs) -> requests.Response:
    """POST with rate-limit-aware retry. See :func:`request_with_retry`."""
    return request_with_retry("POST", url, **kwargs)


def get_with_retry(url: str, **kwargs) -> requests.Response:
    """GET with rate-limit-aware retry. See :func:`request_with_retry`."""
    return request_with_retry("GET", url, **kwargs)
