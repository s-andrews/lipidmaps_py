"""Locate (and optionally download/cache) the molecule database for annotation.

imzML references only metadata ontologies (psi-ms.obo / unit.obo), not a compound
database, so auto-annotation needs a molecule list (CoreMetabolome-style
``id,name,formula`` table). This resolver implements a "use it if present, else
provision and cache" policy so the CLI works without bundling the 3.6 MB file.

Resolution order:
  1. explicit ``arg`` path
  2. ``$LIPIDMAPS_METABOLOME_DB``
  3. user cache ``~/.cache/lipidmaps/core_metabolome_v3.csv``
  4. in-repo ``tests/data/core_metabolome_v3.csv`` (when running from source)
  5. download from ``download_url`` / ``$LIPIDMAPS_METABOLOME_DB_URL`` into the cache

No hardcoded default URL: there is no stable public CoreMetabolome download, so the
user must supply one (or a path) — otherwise a clear error explains how to obtain it.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

DEFAULT_DB_NAME = "core_metabolome_v3.csv"
ENV_DB_PATH = "LIPIDMAPS_METABOLOME_DB"
ENV_DB_URL = "LIPIDMAPS_METABOLOME_DB_URL"


def _cache_path() -> Path:
    root = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(root) / "lipidmaps" / DEFAULT_DB_NAME


def _repo_fallback() -> Optional[Path]:
    """tests/data/<db> relative to the source tree, if present (editable installs)."""
    # src/lipidmaps/data/annotation/db_provision.py -> parents[4] == repo root
    try:
        repo_root = Path(__file__).resolve().parents[4]
    except IndexError:  # pragma: no cover
        return None
    candidate = repo_root / "tests" / "data" / DEFAULT_DB_NAME
    return candidate if candidate.exists() else None


def _download(url: str, dest: Path) -> Path:
    from ...http_utils import get_with_retry

    logger.info("Downloading metabolome database from %s", url)
    resp = get_with_retry(url, timeout=120)
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(resp.content)
    logger.info("Cached metabolome database at %s (%d bytes)", dest, len(resp.content))
    return dest


def resolve_metabolome_db(
    arg: Optional[str] = None, download_url: Optional[str] = None
) -> Path:
    """Return a path to the molecule database, provisioning into the cache if needed.

    Raises ``FileNotFoundError`` with guidance when nothing resolves and no URL is
    available to download from.
    """
    if arg:
        p = Path(arg).expanduser()
        if p.exists():
            return p
        raise FileNotFoundError(f"--database path does not exist: {p}")

    env_path = os.environ.get(ENV_DB_PATH)
    if env_path and Path(env_path).expanduser().exists():
        return Path(env_path).expanduser()

    cache = _cache_path()
    if cache.exists():
        return cache

    repo = _repo_fallback()
    if repo is not None:
        return repo

    url = download_url or os.environ.get(ENV_DB_URL)
    if url:
        return _download(url, cache)

    raise FileNotFoundError(
        "No molecule database found for auto-annotation. Provide one of:\n"
        f"  --database PATH, or ${ENV_DB_PATH}=PATH (a CoreMetabolome-style "
        "id,name,formula table),\n"
        f"  or --database-url URL / ${ENV_DB_URL}=URL to download and cache it at "
        f"{cache}.\n"
        "CoreMetabolome is available from METASPACE (https://metaspace2020.eu)."
    )
