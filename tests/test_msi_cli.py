"""Tests for the lipidmaps-msi CLI and molecule-DB provisioning (offline)."""

from pathlib import Path

import pytest

# These modules import cleanly without the optional MSI deps (pyimzml/matplotlib are
# imported lazily inside the code paths that need them).
from lipidmaps import msi_cli
from lipidmaps.data.annotation import db_provision
from lipidmaps.data.models.lmsd import LMSD
from lipidmaps.data.models.refmet import RefMet
from lipidmaps.data.models.sample import LipidDataset

REPO_ROOT = Path(__file__).resolve().parents[1]
DEMO_IMZML = REPO_ROOT / "tests" / "data" / "inputs" / "demo" / "msi" / "demo_brain.imzML"
DB_PATH = REPO_ROOT / "tests" / "data" / "core_metabolome_v3.csv"


@pytest.fixture
def _offline(monkeypatch):
    """Stub all network-backed lookups so the CLI runs fully offline."""
    monkeypatch.setattr(RefMet, "validate_metabolite_names", staticmethod(lambda names, *a, **k: []))
    monkeypatch.setattr(LMSD, "get_lm_ids_by_name", staticmethod(lambda names: []))
    monkeypatch.setattr(LMSD, "get_molecules_by_lm_id", staticmethod(lambda lmids: []))


def test_msi_cli_end_to_end(tmp_path, _offline):
    pytest.importorskip("pyimzml")
    pytest.importorskip("matplotlib")
    out = tmp_path / "out"
    rc = msi_cli.main([
        str(DEMO_IMZML),
        "--database", str(DB_PATH),
        "--no-fetch-reactions",
        "--out", str(out),
    ])
    assert rc == 0

    ds_json = out / "processed_dataset.json"
    assert ds_json.exists()
    dataset = LipidDataset.model_validate_json(ds_json.read_text(encoding="utf-8"))
    assert dataset.is_spatial
    assert len(dataset.spatial_samples()) > 100

    ions = (out / "ions.csv").read_text(encoding="utf-8").splitlines()
    assert ions[0].startswith("input_name,")
    assert len(ions) > 1  # header + ≥1 ion
    # Auto-annotated ions are labeled formula [adduct].
    assert any(" [" in line for line in ions[1:])

    pngs = list((out / "images").glob("*.png"))
    assert len(pngs) >= 1
    assert all(p.stat().st_size > 0 for p in pngs)


def test_msi_cli_missing_imzml(tmp_path, _offline):
    rc = msi_cli.main([str(tmp_path / "nope.imzML"), "--no-fetch-reactions"])
    assert rc == 2


def test_msi_cli_html_output(tmp_path, _offline):
    pytest.importorskip("pyimzml")
    pytest.importorskip("matplotlib")
    pytest.importorskip("plotly")
    out = tmp_path / "out"
    rc = msi_cli.main([
        str(DEMO_IMZML), "--database", str(DB_PATH),
        "--no-fetch-reactions", "--html", "--max-images", "2", "--out", str(out),
    ])
    assert rc == 0
    htmls = list((out / "images").glob("*.html"))
    assert len(htmls) >= 1 and all(h.stat().st_size > 0 for h in htmls)
    # plotly.js is emitted once and shared by the HTML files.
    assert (out / "images" / "plotly.min.js").exists()


def test_pixel_size_parsing():
    from lipidmaps.data.ingestion.imzml_reader import ImzMLIngestion

    class _P:
        imzmldict = {"pixel size (x)": 10.0, "pixel size (y)": 12.0}

    assert ImzMLIngestion._pixel_size(_P()) == (10.0, 12.0)

    class _Q:
        imzmldict = {"max count of pixels x": 5}

    assert ImzMLIngestion._pixel_size(_Q()) is None


def test_resolve_db_prefers_explicit_arg(tmp_path):
    f = tmp_path / "db.csv"
    f.write_text("id\tname\tformula\n", encoding="utf-8")
    assert db_provision.resolve_metabolome_db(str(f)) == f


def test_resolve_db_missing_arg_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        db_provision.resolve_metabolome_db(str(tmp_path / "absent.csv"))


def test_resolve_db_env_and_cache_and_error(tmp_path, monkeypatch):
    # No arg; disable repo fallback and point cache at an empty dir.
    monkeypatch.setattr(db_provision, "_repo_fallback", lambda: None)
    monkeypatch.setattr(db_provision, "_cache_path", lambda: tmp_path / "cache" / "db.csv")
    monkeypatch.delenv(db_provision.ENV_DB_PATH, raising=False)
    monkeypatch.delenv(db_provision.ENV_DB_URL, raising=False)
    # Nothing available -> clear error.
    with pytest.raises(FileNotFoundError):
        db_provision.resolve_metabolome_db(None)

    # Env var pointing at a real file wins.
    envdb = tmp_path / "env.csv"
    envdb.write_text("id\tname\tformula\n", encoding="utf-8")
    monkeypatch.setenv(db_provision.ENV_DB_PATH, str(envdb))
    assert db_provision.resolve_metabolome_db(None) == envdb


def test_resolve_db_repo_fallback_is_real():
    # In the source tree the repo fallback should find the committed DB.
    assert db_provision._repo_fallback() == DB_PATH
