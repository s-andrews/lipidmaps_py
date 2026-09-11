"""Unit tests for spatial helpers (no network, no pyimzml required)."""

import numpy as np
import pytest

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


def test_reaction_effect_sizes_between_groups():
    from lipidmaps.data.utils.spatial import reaction_effect_sizes
    from lipidmaps.data.models.reaction import ReactionData, CompoundComponent

    def sm(name, group, x):
        return SampleMetadata(sample_name=name, group=group, coordinates=PixelCoordinate(x=x, y=0))

    samples = [sm("a1", "ctrl", 0), sm("a2", "ctrl", 1), sm("b1", "dis", 0), sm("b2", "dis", 1)]
    react = QuantifiedLipid(input_name="R", lm_id="LM_R",
                            values={"a1": 10.0, "a2": 10.0, "b1": 10.0, "b2": 10.0})
    prod = QuantifiedLipid(input_name="P", lm_id="LM_P",
                           values={"a1": 10.0, "a2": 10.0, "b1": 40.0, "b2": 40.0})
    rx = ReactionData(reaction_name="R -> P",
                      reactants=[CompoundComponent(compound_lm_id="LM_R")],
                      products=[CompoundComponent(compound_lm_id="LM_P")])
    ds = LipidDataset(samples=samples, lipids=[react, prod], reactions=[rx])

    rows = reaction_effect_sizes(ds, "ctrl", "dis")
    assert len(rows) == 1
    row = rows[0]
    assert row["mean_log2ratio_ctrl"] == 0.0          # log2(10/10)
    assert row["mean_log2ratio_dis"] == 2.0            # log2(40/10)
    assert row["effect_size_delta"] == 2.0


def test_aggregate_replicates_section_and_tile():
    from lipidmaps.data.utils.spatial import aggregate_replicates

    # Two groups, each two z-sections, two pixels per section.
    samples, values = [], {}
    for g in ("ctrl", "dis"):
        for z in (0, 1):
            for x in (0, 1):
                sn = f"{g}:{z}:{x}"
                samples.append(SampleMetadata(sample_name=sn, group=g,
                                              coordinates=PixelCoordinate(x=x, y=0, z=z)))
                values[sn] = 10.0 if g == "ctrl" else 20.0
    lipid = QuantifiedLipid(input_name="R", lm_id="LM_R", values=values)
    ds = LipidDataset(samples=samples, lipids=[lipid])

    # pixel → unchanged
    assert aggregate_replicates(ds, "pixel") is ds

    # section → one unit per (group, z): 2 per group = 4 units, values = mean (unchanged here)
    sec = aggregate_replicates(ds, "section")
    assert len(sec.samples) == 4
    counts = {}
    for s in sec.samples:
        counts[s.group] = counts.get(s.group, 0) + 1
    assert counts == {"ctrl": 2, "dis": 2}
    r = sec.lipids[0]
    assert all(v == 10.0 for k, v in r.values.items() if k.startswith("ctrl"))

    # tile with 4 tiles: units are a subset (few distinct x,y here) but grouped correctly
    tile = aggregate_replicates(ds, "tile", tiles=4)
    assert {s.group for s in tile.samples} == {"ctrl", "dis"}
    assert tile.reactions == ds.reactions  # reactions carried over


def test_reaction_gene_labels():
    from lipidmaps.data.utils.spatial import reaction_gene_labels
    from lipidmaps.data.models.reaction import ReactionData

    rx = ReactionData(
        reaction_name="R -> P",
        genes=[{"gene_name": "COMT", "uniprot_id": "P21964"}, {"uniprot_id": "Q9Y2"}],
        proteins=[{"ec_number": "2.7.8.27"}],
    )
    labels = reaction_gene_labels(rx)
    assert "COMT" in labels
    assert "Q9Y2" in labels          # falls back to uniprot when no gene_name
    assert "EC 2.7.8.27" in labels


def test_prominent_reactions_and_cloud():
    from lipidmaps.data.utils.spatial import prominent_reactions, reaction_ratio_cloud
    from lipidmaps.data.models.reaction import ReactionData, CompoundComponent

    samples = [SampleMetadata(sample_name=f"px{i}", group="t", coordinates=PixelCoordinate(x=i, y=0))
               for i in range(3)]
    react = QuantifiedLipid(input_name="R", lm_id="LM_R", values={"px0": 10.0, "px1": 10.0, "px2": 10.0})
    prod = QuantifiedLipid(input_name="P", lm_id="LM_P", values={"px0": 20.0, "px1": 40.0, "px2": 0.0})
    rx = ReactionData(reaction_name="R -> P",
                      reactants=[CompoundComponent(compound_lm_id="LM_R")],
                      products=[CompoundComponent(compound_lm_id="LM_P")])
    ds = LipidDataset(samples=samples, lipids=[react, prod], reactions=[rx])

    ranked = prominent_reactions(ds)
    assert ranked and ranked[0][0] is rx
    assert ranked[0][3] == 30.0 + 60.0  # sum reactant + product signal

    coords, vals = reaction_ratio_cloud(ds, "R", "P")
    assert len(coords) == 3
    assert vals[0] == 1.0            # log2(20/10)
    assert vals[2] is None           # product 0 -> undefined


def test_normalize_per_section():
    from lipidmaps.data.utils.spatial import normalize_per_section

    coords = [(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1)]
    values = [10.0, 20.0, 100.0, 300.0]
    norm = normalize_per_section(coords, values)
    # z=0: 10..20 → 0,1 ; z=1: 100..300 → 0,1  (each section scaled to itself)
    assert norm == [0.0, 1.0, 0.0, 1.0]
    # None preserved
    assert normalize_per_section([(0, 0, 0)], [None]) == [None]


def test_tile_aligned_delta():
    from lipidmaps.data.utils.spatial import tile_aligned_delta
    from lipidmaps.data.models.reaction import ReactionData

    # Independent 2x1 grids per group; tile grid n=1 → whole slice is one tile.
    samples, rvals, pvals = [], {}, {}
    for g, pfac in (("ctrl", 1.0), ("dis", 4.0)):
        for x in (0, 1):
            sn = f"{g}:{x}"
            samples.append(SampleMetadata(sample_name=sn, group=g, coordinates=PixelCoordinate(x=x, y=0, z=0)))
            rvals[sn] = 10.0
            pvals[sn] = 10.0 * pfac
    react = QuantifiedLipid(input_name="R", lm_id="LM_R", values=rvals)
    prod = QuantifiedLipid(input_name="P", lm_id="LM_P", values=pvals)
    ds = LipidDataset(samples=samples, lipids=[react, prod],
                      reactions=[ReactionData(reaction_name="R->P")])

    res = tile_aligned_delta(ds, "R", "P", "ctrl", "dis", n=1)
    # ctrl log2(1)=0, dis log2(4)=2 → Δ=2 for the single aligned tile (z0,row0,col0)
    assert res["delta"][(0, 0, 0)] == 2.0


def test_tile_delta_figure_single_cube():
    pytest.importorskip("plotly")
    from lipidmaps.data.utils.spatial_html import tile_delta_figure

    fig = tile_delta_figure({"delta": {(0, 0, 0): 2.0, (0, 0, 1): -1.0}, "n": 2},
                            "R -> P Δ", gene_text="COMT")
    assert fig is not None
    assert len(fig.data) == 1        # one cube of aligned tiles


def test_reaction_volume_figure_tiled():
    pytest.importorskip("plotly")
    from lipidmaps.data.utils.spatial_html import reaction_volume_figure

    coords = [(0, 0, 0), (1, 0, 1)]
    fig = reaction_volume_figure(
        [("control", coords, [0.0, 1.0]), ("disease", coords, [1.0, 2.0])],
        "R -> P", gene_text="COMT, EC 2.7.8.27", z_spacing_um=20.0,
    )
    assert fig is not None
    assert len(fig.data) == 2        # one 3D scene per stack (tiled)


def test_compare_stack_reactions_structure():
    from lipidmaps.data.utils.spatial import compare_stack_reactions
    from lipidmaps.data.models.reaction import ReactionData, CompoundComponent

    samples = [
        SampleMetadata(sample_name=f"{g}:{i}", group=g, coordinates=PixelCoordinate(x=i, y=0))
        for g in ("ctrl", "dis") for i in range(3)
    ]
    react = QuantifiedLipid(input_name="R", lm_id="LM_R",
                            values={s.sample_name: 10.0 for s in samples})
    prod = QuantifiedLipid(input_name="P", lm_id="LM_P",
                           values={s.sample_name: (30.0 if s.group == "dis" else 10.0) for s in samples})
    rx = ReactionData(reaction_name="R -> P",
                      reactants=[CompoundComponent(compound_lm_id="LM_R")],
                      products=[CompoundComponent(compound_lm_id="LM_P")])
    ds = LipidDataset(samples=samples, lipids=[react, prod], reactions=[rx])

    result = compare_stack_reactions(ds, "ctrl", "dis", threshold=0.05)
    assert set(result) == {"rows", "graph"}
    assert "nodes" in result["graph"] and "edges" in result["graph"]
    for r in result["rows"]:
        assert {"reaction", "z_score", "direction", "significant"} <= set(r)


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
