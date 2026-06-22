"""Tests for the warm BioPAN service orchestration.

These exercise the service-specific logic (request validation, the warm session
registry, background jobs and status files, Stage C waiting on Stage B) without
network or heavy compute by stubbing ``run_session`` — the export work itself is
already covered by the CLI / exporter tests and is byte-identical here since the
service calls the same function.
"""

import importlib.util
import json
import os
import sys
import time

import pytest

from lipidmaps.biopan_cli import RunParams, RunResult

# The service is a standalone BioPAN script that lives with the web tool
# (laravel/scripts/biopan/biopan_service.py), not inside the lipidmaps package.
# Load it by path so these tests can run from the lipidmaps_py venv (which has the
# lipidmaps deps the script imports). Overridable via BIOPAN_SERVICE_PATH.
_SERVICE_PATH = os.environ.get(
    "BIOPAN_SERVICE_PATH",
    "/lipidmaps/lipidmaps/www/laravel/scripts/biopan/biopan_service.py",
)
if not os.path.exists(_SERVICE_PATH):
    pytest.skip(f"biopan_service.py not found at {_SERVICE_PATH}", allow_module_level=True)

_spec = importlib.util.spec_from_file_location("biopan_service", _SERVICE_PATH)
biopan_service = importlib.util.module_from_spec(_spec)
# Register before exec so the @dataclass annotation resolution (which looks the
# module up in sys.modules by __module__) works.
sys.modules[_spec.name] = biopan_service
_spec.loader.exec_module(biopan_service)

BioPANRequestError = biopan_service.BioPANRequestError
_params_from_payload = biopan_service._params_from_payload
_resolve_session_dir = biopan_service._resolve_session_dir
_resolve_status_file = biopan_service._resolve_status_file
handle_run = biopan_service.handle_run


@pytest.fixture
def session_dir(tmp_path, monkeypatch):
    """A valid session dir under a temp SESSION_ROOT, with a fresh registry."""
    root = tmp_path / "biopan"
    sess = root / "abcdef0123456789"
    (sess / "config").mkdir(parents=True)
    (sess / "input").mkdir()
    (sess / "input" / "input.csv").write_text("x\n", encoding="utf-8")
    monkeypatch.setattr(biopan_service, "SESSION_ROOT", root)
    monkeypatch.setattr(biopan_service, "_registry", biopan_service._Registry(8))
    return sess


def _stub_run_session(monkeypatch, recorder=None, output=None, reprocessed=False, raises=None):
    def fake(params, *, dataset=None, exporter=None):
        if recorder is not None:
            recorder.append({"params": params, "dataset": dataset, "exporter": exporter})
        if raises is not None:
            raise raises
        return RunResult(
            dataset=object(),
            exporter=object(),
            reprocessed=reprocessed,
            output=output if output is not None else {"written": {"summary.json": "ok"}},
        )

    monkeypatch.setattr(biopan_service, "run_session", fake)


# --- request validation ------------------------------------------------------

def test_resolve_session_dir_rejects_outside_root(session_dir, monkeypatch):
    with pytest.raises(BioPANRequestError):
        _resolve_session_dir({"session_dir": "/etc"})


def test_resolve_session_dir_accepts_valid(session_dir):
    assert _resolve_session_dir({"session_dir": str(session_dir)}) == session_dir.resolve()


def test_resolve_status_file_must_be_under_session(session_dir):
    good = _resolve_status_file(session_dir.resolve(), str(session_dir / "config" / "reactions_status.json"))
    assert good is not None
    with pytest.raises(BioPANRequestError):
        _resolve_status_file(session_dir.resolve(), "/tmp/evil.json")


def test_params_from_payload_maps_fields(session_dir):
    params = _params_from_payload(
        {
            "build_view": True,
            "scope": "tables",
            "family": "pathway",
            "level": "species",
            "disease_group": "d",
            "control_group": "c",
            "threshold": "0.01",
            "paired": "TRUE",
            "sample_group": ["S1=c", "S2=d"],
        },
        session_dir.resolve(),
    )
    assert isinstance(params, RunParams)
    assert params.build_view and params.scope == "tables" and params.family == "pathway"
    assert params.level == "species" and params.threshold == 0.01 and params.paired is True
    assert params.sample_group == ["S1=c", "S2=d"]


def test_params_from_payload_rejects_bad_enum(session_dir):
    with pytest.raises(BioPANRequestError):
        _params_from_payload({"scope": "bogus"}, session_dir.resolve())


# --- foreground requests -----------------------------------------------------

def test_handle_run_build_view_returns_sorted_keys(session_dir, monkeypatch):
    _stub_run_session(monkeypatch, output={"built_view": {"b.json": 1, "a.json": 2}})
    status, body = handle_run({"session_dir": str(session_dir), "build_view": True,
                               "disease_group": "d", "control_group": "c"})
    assert status == 200
    assert body == {"status": "ok", "built_view": ["a.json", "b.json"]}


def test_warm_registry_reuse(session_dir, monkeypatch):
    calls = []
    _stub_run_session(monkeypatch, recorder=calls)
    payload = {"session_dir": str(session_dir), "summary_only": True}

    handle_run(payload)
    handle_run(payload)

    # First call has a cold registry (no warm objects); second reuses what the
    # first stored.
    assert calls[0]["dataset"] is None and calls[0]["exporter"] is None
    assert calls[1]["dataset"] is not None and calls[1]["exporter"] is not None


def test_reprocess_does_not_reuse_warm(session_dir, monkeypatch):
    calls = []
    _stub_run_session(monkeypatch, recorder=calls, reprocessed=True)
    # Warm the registry first.
    handle_run({"session_dir": str(session_dir), "summary_only": True})
    # A reprocess (csv_path present) must not inject the stale warm objects.
    handle_run({"session_dir": str(session_dir), "csv_path": str(session_dir / "input" / "input.csv")})
    assert calls[1]["dataset"] is None and calls[1]["exporter"] is None


# --- background jobs + status files ------------------------------------------

def _wait_status(path, expected, timeout=5.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            if data.get("status") == expected:
                return True
        except (OSError, ValueError):
            # File may not exist yet or be mid-write; keep polling.
            pass
        time.sleep(0.02)
    return False


def test_background_writes_ready_status(session_dir, monkeypatch):
    _stub_run_session(monkeypatch)
    status_file = session_dir / "config" / "reactions_status.json"
    status, body = handle_run({
        "session_dir": str(session_dir),
        "fetch_reactions": True,
        "summary_only": True,
        "background": True,
        "status_file": str(status_file),
    })
    assert status == 202 and body == {"status": "accepted"}
    assert _wait_status(status_file, "ready")


def test_background_writes_error_status_on_failure(session_dir, monkeypatch):
    _stub_run_session(monkeypatch, raises=RuntimeError("boom"))
    status_file = session_dir / "config" / "reactions_status.json"
    handle_run({
        "session_dir": str(session_dir),
        "background": True,
        "status_file": str(status_file),
    })
    assert _wait_status(status_file, "error")


def test_wait_for_reactions_returns_when_status_ready(session_dir, monkeypatch):
    _stub_run_session(monkeypatch)
    # Stage B already finished.
    (session_dir / "config" / "reactions_status.json").write_text('{"status":"ready"}', encoding="utf-8")
    quant_status = session_dir / "config" / "quantification_status.json"
    handle_run({
        "session_dir": str(session_dir),
        "sample_group": ["S1=c"],
        "lazy_bundle": True,
        "wait_for_reactions": True,
        "background": True,
        "status_file": str(quant_status),
    })
    assert _wait_status(quant_status, "ready")
