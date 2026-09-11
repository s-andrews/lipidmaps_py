"""Unit tests for spatial helpers (no network, no pyimzml required)."""

import numpy as np

from lipidmaps.data.models.sample import (
    LipidDataset,
    PixelCoordinate,
    QuantifiedLipid,
    SampleMetadata,
)
from lipidmaps.data.utils.spatial import (
    lipid_spatial_series,
    values_to_grid,
    voronoi_regions,
    z_slices,
)


def _grid_coords(n):
    return [(x, y, 0) for y in range(n) for x in range(n)]


def test_values_to_grid_places_values_by_coordinate():
    coords = _grid_coords(3)
    values = list(range(9))
    grid, mask, xr, yr = values_to_grid(coords, values)
    assert grid.shape == (3, 3)
    assert xr == (0, 2) and yr == (0, 2)
    # grid[y, x]; pixel (x=2, y=0) had value 2.
    assert grid[0, 2] == 2
    assert grid[2, 0] == 6
    assert mask.all()


def test_values_to_grid_marks_absent_pixels():
    coords = [(0, 0, 0), (1, 0, 0)]
    values = [5.0, None]
    grid, mask, _, _ = values_to_grid(coords, values)
    assert grid[0, 0] == 5.0
    assert mask[0, 0] and not mask[0, 1]
    assert np.isnan(grid[0, 1])


def test_values_to_grid_filters_by_z():
    coords = [(0, 0, 0), (0, 0, 1)]
    values = [1.0, 2.0]
    grid, _, _, _ = values_to_grid(coords, values, z=1)
    assert grid.shape == (1, 1)
    assert grid[0, 0] == 2.0


def test_voronoi_regions_one_per_pixel():
    coords = _grid_coords(3)
    values = list(range(9))
    regions = voronoi_regions(coords, values)
    assert len(regions) == 9
    # Every region is a real polygon clipped to the bbox.
    assert all(r["polygon"].shape[0] >= 3 for r in regions)
    assert {(r["x"], r["y"]) for r in regions} == {(c[0], c[1]) for c in coords}


def test_voronoi_regions_too_few_points():
    assert voronoi_regions([(0, 0, 0), (1, 1, 0)]) == []


def test_voronoi_regions_include_z():
    regions = voronoi_regions(_grid_coords(3), list(range(9)), z=0)
    assert all(r["z"] == 0 for r in regions)
    # z filtering + reported z stay consistent for a non-zero slice.
    coords = [(x, y, 2) for y in range(3) for x in range(3)]
    regions2 = voronoi_regions(coords, list(range(9)), z=2)
    assert regions2 and all(r["z"] == 2 for r in regions2)


def test_z_slices():
    coords = [(0, 0, 0), (1, 0, 2), (0, 1, 0)]
    assert z_slices(coords) == [0, 2]


def test_lipid_spatial_series_reads_coordinates_and_values():
    s1 = SampleMetadata(sample_name="px_a", group="t", coordinates=PixelCoordinate(x=0, y=0))
    s2 = SampleMetadata(sample_name="px_b", group="t", coordinates=PixelCoordinate(x=1, y=0))
    plain = SampleMetadata(sample_name="not_spatial", group="t")
    lipid = QuantifiedLipid(input_name="PC 34:1", values={"px_a": 10.0, "px_b": 20.0})
    ds = LipidDataset(samples=[s1, s2, plain], lipids=[lipid])

    coords, values = lipid_spatial_series(ds, "PC 34:1")
    assert coords == [(0, 0, 0), (1, 0, 0)]  # plain sample skipped (no coordinates)
    assert values == [10.0, 20.0]
    assert ds.is_spatial is True
    assert len(ds.spatial_samples()) == 2


def test_dataset_without_coordinates_is_not_spatial():
    ds = LipidDataset(
        samples=[SampleMetadata(sample_name="s1", group="g")],
        lipids=[QuantifiedLipid(input_name="PC 34:1", values={"s1": 1.0})],
    )
    assert ds.is_spatial is False
    assert ds.spatial_samples() == []


def test_single_layer_is_not_3d():
    samples = [SampleMetadata(sample_name=f"px{i}", group="t",
                              coordinates=PixelCoordinate(x=i, y=0, z=1)) for i in range(3)]
    ds = LipidDataset(samples=samples, lipids=[])
    assert ds.z_layers() == [1]
    assert ds.z_slice_count == 1
    assert ds.is_3d is False


def test_multi_layer_is_3d():
    samples = [SampleMetadata(sample_name=f"px{z}_{i}", group="t",
                              coordinates=PixelCoordinate(x=i, y=0, z=z))
               for z in (0, 1, 2) for i in range(2)]
    ds = LipidDataset(samples=samples, lipids=[])
    assert ds.z_layers() == [0, 1, 2]
    assert ds.z_slice_count == 3
    assert ds.is_3d is True


def test_ion_display_name_prefers_resolved_name():
    from lipidmaps.data.utils.spatial import ion_display_full, ion_display_name

    resolved = QuantifiedLipid(input_name="C42H82NO8P [M+H]+", standardized_name="PC 34:1", values={})
    assert ion_display_name(resolved) == "PC 34:1"
    assert ion_display_full(resolved) == "PC 34:1 (C42H82NO8P [M+H]+)"

    only_cand = QuantifiedLipid(input_name="C6H12O6 [M+H]+", values={},
                                annotation_candidates=[{"name": "Glucose", "id": "H1"}])
    assert ion_display_name(only_cand) == "Glucose"

    plain = QuantifiedLipid(input_name="PC 34:1", values={})
    assert ion_display_name(plain) == "PC 34:1"
    assert ion_display_full(plain) == "PC 34:1"
