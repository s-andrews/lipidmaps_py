"""Spatial helpers for mass-spectrometry-imaging (MSI) datasets.

Pure, Streamlit-free functions so they can be unit-tested without a UI:

- :func:`lipid_spatial_series` pulls a lipid's per-pixel values + coordinates out
  of a :class:`~lipidmaps.data.models.sample.LipidDataset`.
- :func:`values_to_grid` rasterizes one z-slice into a dense 2D array for an
  image/heatmap renderer (e.g. Plotly ``px.imshow``).
- :func:`voronoi_regions` builds per-pixel Voronoi polygons (clipped to the data
  bounding box) for irregular/sparse pixel grids, using ``scipy.spatial.Voronoi``
  (already a core dependency).

Coordinates are ``(x, y, z)``; z defaults to 0 (2D acquisition).
"""

from __future__ import annotations

import logging
from typing import List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

Coord = Tuple[int, int, int]


def lipid_spatial_series(
    dataset, lipid, z: Optional[int] = None
) -> Tuple[List[Coord], List[Optional[float]]]:
    """Return ``(coords, values)`` for one lipid across the dataset's MSI pixels.

    ``lipid`` may be a ``QuantifiedLipid`` or an ``input_name`` string. Only samples
    that carry coordinates are considered; when ``z`` is given, only that slice.
    """
    name = getattr(lipid, "input_name", lipid)
    target = None
    for lp in dataset.lipids:
        if lp is lipid or lp.input_name == name:
            target = lp
            break
    if target is None:
        return [], []

    coords: List[Coord] = []
    values: List[Optional[float]] = []
    for sample in dataset.samples:
        c = sample.coordinates
        if c is None:
            continue
        if z is not None and c.z != z:
            continue
        coords.append((c.x, c.y, c.z))
        values.append(target.values.get(sample.sample_name))
    return coords, values


def z_slices(coords: Sequence[Coord]) -> List[int]:
    """Sorted unique z values present in ``coords``."""
    return sorted({int(c[2]) if len(c) > 2 else 0 for c in coords})


def values_to_grid(
    coords: Sequence[Coord],
    values: Sequence[Optional[float]],
    z: int = 0,
    fill: float = np.nan,
) -> Tuple[np.ndarray, np.ndarray, Tuple[int, int], Tuple[int, int]]:
    """Rasterize a z-slice into a dense 2D array indexed ``grid[y, x]``.

    Returns ``(grid, mask, (x_min, x_max), (y_min, y_max))`` where ``mask`` marks
    pixels that had a (non-None) value. Absent pixels are ``fill`` (NaN by default),
    which renders as blank in most image renderers.
    """
    pts = [
        (c, v)
        for c, v in zip(coords, values)
        if (len(c) < 3 or int(c[2]) == z)
    ]
    if not pts:
        empty = np.zeros((0, 0))
        return empty, empty.astype(bool), (0, 0), (0, 0)

    xs = [int(c[0]) for c, _ in pts]
    ys = [int(c[1]) for c, _ in pts]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    width = x_max - x_min + 1
    height = y_max - y_min + 1

    grid = np.full((height, width), fill, dtype=float)
    mask = np.zeros((height, width), dtype=bool)
    for c, v in pts:
        row = int(c[1]) - y_min
        col = int(c[0]) - x_min
        if v is not None:
            grid[row, col] = float(v)
            mask[row, col] = True
    return grid, mask, (x_min, x_max), (y_min, y_max)


def _voronoi_finite_polygons_2d(vor, radius: Optional[float] = None):
    """Reconstruct finite Voronoi regions in 2D (standard recipe).

    Returns ``(regions, vertices)`` where each region is a list of indices into
    ``vertices``. Infinite ridges are projected out to ``radius`` so boundary cells
    become finite polygons (later clipped to the bounding box).
    """
    if vor.points.shape[1] != 2:
        raise ValueError("Voronoi reconstruction requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()
    center = vor.points.mean(axis=0)
    if radius is None:
        radius = np.ptp(vor.points, axis=0).max() * 2

    # Map each point to the ridges around it.
    all_ridges: dict = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p1, region_idx in enumerate(vor.point_region):
        vertices = vor.regions[region_idx]
        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue

        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]
        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                continue
            # Compute the missing endpoint of an infinite ridge.
            t = vor.points[p2] - vor.points[p1]
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])
            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius
            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # Order region vertices counter-clockwise.
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = [new_region[i] for i in np.argsort(angles)]
        new_regions.append(new_region)

    return new_regions, np.asarray(new_vertices)


def _clip_polygon(polygon: np.ndarray, bbox: Tuple[float, float, float, float]) -> np.ndarray:
    """Clip a convex polygon to an axis-aligned rectangle (Sutherland-Hodgman)."""
    x_min, y_min, x_max, y_max = bbox
    # Each edge: (inside-test, intersect) for the four rectangle boundaries.
    edges = [
        (lambda p: p[0] >= x_min, 0, x_min),
        (lambda p: p[0] <= x_max, 0, x_max),
        (lambda p: p[1] >= y_min, 1, y_min),
        (lambda p: p[1] <= y_max, 1, y_max),
    ]
    output = [np.asarray(p, dtype=float) for p in polygon]
    for inside, axis, bound in edges:
        if not output:
            break
        cur = output
        output = []
        for i in range(len(cur)):
            a = cur[i - 1]
            b = cur[i]
            a_in = inside(a)
            b_in = inside(b)
            if b_in:
                if not a_in:
                    output.append(_intersect(a, b, axis, bound))
                output.append(b)
            elif a_in:
                output.append(_intersect(a, b, axis, bound))
    return np.asarray(output) if output else np.zeros((0, 2))


def _intersect(a: np.ndarray, b: np.ndarray, axis: int, bound: float) -> np.ndarray:
    """Point where segment a-b crosses the line ``coord[axis] == bound``."""
    other = 1 - axis
    denom = b[axis] - a[axis]
    if denom == 0:
        return a.copy()
    t = (bound - a[axis]) / denom
    p = np.empty(2)
    p[axis] = bound
    p[other] = a[other] + t * (b[other] - a[other])
    return p


def voronoi_regions(
    coords: Sequence[Coord],
    values: Optional[Sequence[Optional[float]]] = None,
    z: int = 0,
) -> List[dict]:
    """Voronoi polygons (clipped to the data bbox) for a z-slice.

    Returns one dict per pixel: ``{"x", "y", "value", "polygon"}`` where
    ``polygon`` is an ``(N, 2)`` array of clipped vertices. Falls back to an empty
    list when there are too few distinct points for a tessellation (< 4).
    """
    from scipy.spatial import Voronoi

    if values is None:
        values = [None] * len(coords)
    pts = [
        (int(c[0]), int(c[1]), v)
        for c, v in zip(coords, values)
        if (len(c) < 3 or int(c[2]) == z)
    ]
    if len(pts) < 4:
        logger.debug("voronoi_regions: only %d points; skipping tessellation", len(pts))
        return []

    points = np.array([[p[0], p[1]] for p in pts], dtype=float)
    vor = Voronoi(points)
    regions, vertices = _voronoi_finite_polygons_2d(vor)

    # Pad the bounding box slightly so edge cells have visible area.
    x_min, y_min = points.min(axis=0) - 0.5
    x_max, y_max = points.max(axis=0) + 0.5
    bbox = (x_min, y_min, x_max, y_max)

    out: List[dict] = []
    for (px, py, val), region in zip(pts, regions):
        poly = vertices[region]
        clipped = _clip_polygon(poly, bbox)
        out.append({"x": px, "y": py, "value": val, "polygon": clipped})
    return out


def aggregate_by_region(dataset, lipid, regions) -> dict:
    """Scaffold: mean value of ``lipid`` per named pixel region.

    ``regions`` maps region_name -> set/list of pixel sample_names. Returns
    region_name -> mean intensity. Region/ROI-based reaction comparison builds on
    this; it is intentionally minimal for now.
    """
    name = getattr(lipid, "input_name", lipid)
    target = next((lp for lp in dataset.lipids if lp.input_name == name), None)
    if target is None:
        return {}
    out: dict = {}
    for region_name, pixel_names in regions.items():
        vals = [
            target.values.get(pn)
            for pn in pixel_names
            if target.values.get(pn) is not None
        ]
        out[region_name] = float(np.mean(vals)) if vals else None
    return out
