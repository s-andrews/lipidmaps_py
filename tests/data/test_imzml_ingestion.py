"""Tests for imzML ingestion and the import_imzml pipeline (offline).

Uses the committed synthetic MSI fixture under inputs/demo/msi/ and small
in-test imzML files written with pyimzml.ImzMLWriter. Network-dependent
standardization/reactions are disabled (use_refmet/fetch_reactions=False).
"""

from pathlib import Path

import numpy as np
import pytest

from lipidmaps.data.ingestion.imzml_reader import (
    ImzMLIngestion,
    IonAnnotation,
    parse_annotation_csv,
    pixel_name,
)

pytest.importorskip("pyimzml")

FIXTURE_DIR = Path(__file__).resolve().parent / "inputs" / "demo" / "msi"
DEMO_IMZML = FIXTURE_DIR / "demo_brain.imzML"
DEMO_ANNOTATIONS = FIXTURE_DIR / "annotations.csv"


def _write_imzml(tmp_path, mz_axis, spectra):
    """spectra: list of (intensities, (x, y, z)). Returns the .imzML path."""
    from pyimzml.ImzMLWriter import ImzMLWriter

    path = tmp_path / "synth.imzML"
    with ImzMLWriter(str(path)) as writer:
        for intensities, coord in spectra:
            writer.addSpectrum(np.asarray(mz_axis, float), np.asarray(intensities, float), coord)
    return path


def test_pixel_name_is_sortable_and_deterministic():
    assert pixel_name(1, 2, 0) == "px_x0001_y0002_z00"
    assert pixel_name(12, 34, 1) == "px_x0012_y0034_z01"


def test_parse_annotation_csv_fixture():
    anns = parse_annotation_csv(DEMO_ANNOTATIONS)
    assert len(anns) >= 5
    by_name = {a.name: a for a in anns}
    assert "PC 34:1" in by_name
    assert by_name["PC 34:1"].mz > 0


def test_parse_annotation_csv_requires_name_and_mz(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("foo,bar\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError):
        parse_annotation_csv(bad)


def test_reader_extracts_ion_within_tolerance(tmp_path):
    mz_axis = [700.5000, 810.6000]
    # 2x1 pixels: intensities differ by x for the first ion.
    spectra = [
        ([10.0, 5.0], (1, 1, 1)),
        ([20.0, 5.0], (2, 1, 1)),
    ]
    path = _write_imzml(tmp_path, mz_axis, spectra)
    anns = [IonAnnotation(name="A", mz=700.5), IonAnnotation(name="B", mz=810.6)]
    result = ImzMLIngestion().read(path, anns, mz_tolerance_ppm=10)

    assert result.pixel_count == 2
    a_vals = result.ion_values["A"]
    assert a_vals[pixel_name(1, 1, 1)] == 10.0
    assert a_vals[pixel_name(2, 1, 1)] == 20.0


def test_reader_returns_zero_outside_tolerance(tmp_path):
    # Peak sits ~700 ppm away from the annotation target -> no match.
    path = _write_imzml(tmp_path, [700.9], [([99.0], (1, 1, 1))])
    anns = [IonAnnotation(name="A", mz=700.5)]
    result = ImzMLIngestion().read(path, anns, mz_tolerance_ppm=10)
    assert result.ion_values["A"][pixel_name(1, 1, 1)] == 0.0


def test_reader_bbox_crop(tmp_path):
    mz_axis = [700.5]
    spectra = [([1.0], (x, 1, 1)) for x in range(1, 6)]
    path = _write_imzml(tmp_path, mz_axis, spectra)
    anns = [IonAnnotation(name="A", mz=700.5)]
    result = ImzMLIngestion().read(path, anns, bbox=(2, 1, 4, 1))
    assert result.pixel_count == 3  # x in {2,3,4}


def test_import_imzml_builds_spatial_dataset():
    from lipidmaps import import_imzml

    data = import_imzml(
        DEMO_IMZML,
        DEMO_ANNOTATIONS,
        use_refmet=False,
        use_headgroups=False,
        fetch_reactions=False,
    )
    ds = data.dataset
    assert ds.is_spatial is True
    assert len(ds.spatial_samples()) == len(ds.samples) > 100
    # Every sample carries a 3D coordinate.
    sample = ds.spatial_samples()[0]
    assert sample.coordinates is not None
    assert isinstance(sample.coordinates.x, int)
    # A known ion has per-pixel intensities.
    pc = next(lp for lp in ds.lipids if lp.input_name == "PC 34:1")
    assert any(v is not None and v > 0 for v in pc.values.values())
