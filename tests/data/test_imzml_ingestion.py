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
    assert "PC(34:1)" in by_name
    assert by_name["PC(34:1)"].mz > 0


def test_parse_annotation_csv_requires_name_and_mz(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("foo,bar\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError):
        parse_annotation_csv(bad)


def test_parse_metabolights_maf(tmp_path):
    """A tab-delimited MetaboLights MAF parses via its standard column names."""
    maf = tmp_path / "m_MTBLS1_maf.tsv"
    maf.write_text(
        "database_identifier\tchemical_formula\tmetabolite_identification\tmass_to_charge\n"
        "\tC27H42O11\tCortolone-3-glucuronide\t525.3386\n"
        "LMFA00000001\tC4H4N2OS\t2-Thiouracil\t145.9975\n"
        "CHEBI:1234\tC6H6O4\tSomeAcid\t124.9992\n"
        "\t\tNoMz\t\n",  # skipped (no mz)
        encoding="utf-8",
    )
    anns = parse_annotation_csv(maf)
    assert len(anns) == 3
    by_name = {a.name: a for a in anns}
    assert by_name["Cortolone-3-glucuronide"].mz == 525.3386
    assert by_name["Cortolone-3-glucuronide"].formula == "C27H42O11"
    # database_identifier used as lm_id only when it's a LIPID MAPS id.
    assert by_name["2-Thiouracil"].lm_id == "LMFA00000001"
    assert by_name["SomeAcid"].lm_id is None  # ChEBI id ignored


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
    pc = next(lp for lp in ds.lipids if lp.input_name == "PC(34:1)")
    assert any(v is not None and v > 0 for v in pc.values.values())


def test_reader_auto_annotates_from_database(tmp_path):
    """With no annotation file, match observed peaks against a molecule database."""
    from lipidmaps.data.annotation.mz_annotator import (
        MetaboliteDatabase,
        POSITIVE_ADDUCTS,
        formula_monoisotopic_mass,
    )

    glucose_mz = POSITIVE_ADDUCTS[0].mz(formula_monoisotopic_mass("C6H12O6"))
    path = _write_imzml(tmp_path, [glucose_mz], [([100.0], (1, 1, 1)), ([50.0], (2, 1, 1))])

    db = MetaboliteDatabase(
        {"C6H12O6": {"mass": formula_monoisotopic_mass("C6H12O6"),
                     "candidates": [{"name": "Glucose", "id": "HMDB1"}]}}
    )
    result = ImzMLIngestion().read(path, annotations=None, database=db, mz_tolerance_ppm=5)

    assert len(result.annotations) == 1
    ann = result.annotations[0]
    assert ann.name == "C6H12O6 [M+H]+"
    assert ann.formula == "C6H12O6"
    assert {c["name"] for c in ann.candidates} == {"Glucose"}
    assert result.ion_values[ann.name][pixel_name(1, 1, 1)] == 100.0


def test_reader_without_annotation_or_database_errors(tmp_path):
    path = _write_imzml(tmp_path, [700.5], [([1.0], (1, 1, 1))])
    with pytest.raises(ValueError):
        ImzMLIngestion().read(path)  # no annotations, no database


def test_import_imzml_stack_builds_3d(tmp_path):
    """Two imzML files stack into one 3D dataset with distinct z-layers."""
    from lipidmaps import import_imzml_stack

    mz_axis = [700.5000, 810.6000]
    # Two 2x1 'sections', each its own imzML file in its own dir.
    d1 = tmp_path / "a"
    d1.mkdir()
    d2 = tmp_path / "b"
    d2.mkdir()
    p1 = _write_imzml(d1, mz_axis, [([10.0, 5.0], (1, 1, 1)), ([20.0, 5.0], (2, 1, 1))])
    p2 = _write_imzml(d2, mz_axis, [([30.0, 5.0], (1, 1, 1)), ([40.0, 5.0], (2, 1, 1))])
    anns = [IonAnnotation(name="A", mz=700.5), IonAnnotation(name="B", mz=810.6)]

    data = import_imzml_stack(
        [p1, p2], annotation_path=anns, z_spacing_um=15.0,
        use_refmet=False, use_headgroups=False, fetch_reactions=False,
    )
    ds = data.dataset
    assert ds.is_3d and ds.z_slice_count == 2
    assert ds.z_layers() == [0, 1]
    assert len(ds.spatial_samples()) == 4  # 2 pixels × 2 layers
    assert ds.z_spacing_um == 15.0
    # Sample names are unique across layers (encode the layer z).
    assert len({s.sample_name for s in ds.samples}) == 4
    a = next(lp for lp in ds.lipids if lp.input_name == "A")
    # Layer-1 pixel (from file b) carries its value under the layer-encoded name.
    assert a.values[pixel_name(1, 1, 1)] == 30.0


def test_import_imzml_groups_two_stacks(tmp_path):
    """Two labelled stacks merge into one dataset with distinct groups."""
    from lipidmaps import import_imzml_groups

    mz_axis = [700.5000, 810.6000]
    da = tmp_path / "a"
    da.mkdir()
    db = tmp_path / "b"
    db.mkdir()
    pa = _write_imzml(da, mz_axis, [([10.0, 5.0], (1, 1, 1)), ([20.0, 5.0], (2, 1, 1))])
    pb = _write_imzml(db, mz_axis, [([30.0, 5.0], (1, 1, 1)), ([40.0, 5.0], (2, 1, 1))])
    anns = [IonAnnotation(name="A", mz=700.5), IonAnnotation(name="B", mz=810.6)]

    data = import_imzml_groups(
        {"control": [pa], "disease": [pb]}, annotation_path=anns,
        use_refmet=False, use_headgroups=False, fetch_reactions=False,
    )
    ds = data.dataset
    assert ds.is_spatial
    assert {s.group for s in ds.samples} == {"control", "disease"}
    assert len(ds.samples) == 4
    # Group-prefixed sample names keep both stacks' identical coords distinct.
    assert all(s.sample_name.startswith(("control:", "disease:")) for s in ds.samples)
    a = next(lp for lp in ds.lipids if lp.input_name == "A")
    # z in the sample name is the file's layer index within the group (0 here), not the imzML z.
    assert a.values["control:" + pixel_name(1, 1, 0)] == 10.0
    assert a.values["disease:" + pixel_name(1, 1, 0)] == 30.0


def test_import_imzml_auto_annotation_end_to_end():
    from lipidmaps import import_imzml

    data = import_imzml(
        DEMO_IMZML,
        annotation_path=None,
        database=str(FIXTURE_DIR.parents[2] / "core_metabolome_v3.csv"),
        mz_tolerance_ppm=5,
        use_refmet=False,
        use_headgroups=False,
        fetch_reactions=False,
    )
    ds = data.dataset
    assert ds.is_spatial
    # Ions are labeled formula [adduct] and carry candidate molecule names.
    assert all(" [" in lp.input_name and "]" in lp.input_name for lp in ds.lipids)
    assert any(lp.annotation_candidates for lp in ds.lipids)


def test_resolve_candidate_lm_ids(monkeypatch):
    """Candidate molecule names are standardized via RefMet to set lm_id/std name."""
    from lipidmaps.data.data_manager import DataManager
    from lipidmaps.data.models.refmet import RefMet
    from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata, PixelCoordinate

    lipid = QuantifiedLipid(
        input_name="C42H82NO8P [M+H]+",
        values={"px_x0001_y0001_z01": 5.0},
        annotation_candidates=[
            {"name": "PC(20:1/14:0)", "id": "HMDB8"},        # -> std name only
            {"name": "PC(16:0/18:1(11Z))", "id": "HMDB9"},   # -> specific lm_id
        ],
    )
    ds = LipidDataset(
        samples=[SampleMetadata(sample_name="px_x0001_y0001_z01", group="t",
                                coordinates=PixelCoordinate(x=1, y=1, z=1))],
        lipids=[lipid],
    )

    def fake_refmet(names):
        assert "PC(16:0/18:1(11Z))" in names
        return [
            {"input_name": "PC(20:1/14:0)", "standardized_name": "PC 34:1", "lm_id": None},
            {"input_name": "PC(16:0/18:1(11Z))", "standardized_name": "PC 16:0/18:1(11Z)",
             "lm_id": "LMGP01010576"},
        ]

    monkeypatch.setattr(RefMet, "validate_metabolite_names", staticmethod(fake_refmet))

    n = DataManager()._resolve_candidate_lm_ids(ds)
    assert n == 1
    # Prefers the candidate that yields a specific lm_id.
    assert lipid.lm_id == "LMGP01010576"
    assert lipid.lm_id_found_by == "candidate"
    assert lipid.standardized_name == "PC 16:0/18:1(11Z)"
