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


def nice_length(value: float) -> float:
    """Round a length up to a 1/2/5×10ⁿ 'nice' value (for scale bars)."""
    import math

    if value <= 0:
        return 1.0
    exp = math.floor(math.log10(value))
    base = value / (10 ** exp)
    nice = 1 if base < 1.5 else 2 if base < 3 else 5 if base < 7 else 10
    return nice * (10 ** exp)


def ion_display_name(lipid) -> str:
    """Friendliest name for an ion: standardized_name > first candidate > input_name.

    Auto-annotated MSI ions are labeled ``formula [adduct]`` in ``input_name``; this
    surfaces the resolved molecule name when one exists.
    """
    name = getattr(lipid, "standardized_name", None)
    if name:
        return name
    candidates = getattr(lipid, "annotation_candidates", None)
    if candidates and candidates[0].get("name"):
        return candidates[0]["name"]
    return lipid.input_name


def ion_display_full(lipid) -> str:
    """``name (formula [adduct])`` when a molecule name is known, else the label."""
    name = ion_display_name(lipid)
    return name if name == lipid.input_name else f"{name} ({lipid.input_name})"


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
        out.append({"x": px, "y": py, "z": z, "value": val, "polygon": clipped})
    return out


def ion_names_for_component_ids(dataset, lm_ids: set) -> List[str]:
    """Measured ion names whose lm_id/generic_lm_id is in ``lm_ids``."""
    if not lm_ids:
        return []
    names = []
    for lp in dataset.lipids:
        if lp.lm_id in lm_ids or lp.generic_lm_id in lm_ids:
            if any(v is not None for v in lp.values.values()):
                names.append(lp.input_name)
    return names


def component_ids(components) -> set:
    """Collect lm_id/generic_lm_id from a reaction's reactant/product components."""
    ids = set()
    for comp in components or []:
        for attr in ("compound_lm_id", "compound_generic_lm_id"):
            val = getattr(comp, attr, None)
            if val:
                ids.add(val)
    return ids


def spatial_reactions(dataset):
    """Reactions (deduped by name) with a measured reactant AND product ion.

    Returns ``(reaction, reactant_ion_names, product_ion_names)`` tuples — the
    reactions meaningful to map over tissue. Shared by the Streamlit UI and the CLI.
    """
    out = []
    seen = set()
    for rx in getattr(dataset, "reactions", None) or []:
        react = ion_names_for_component_ids(dataset, component_ids(rx.reactants))
        prod = ion_names_for_component_ids(dataset, component_ids(rx.products))
        if not (react and prod):
            continue
        label = rx.reaction_name or f"reaction {getattr(rx, 'reaction_id', '')}"
        if label in seen:
            continue
        seen.add(label)
        out.append((rx, react, prod))
    return out


def ratio_series(dataset, reactant_name, product_name, z: Optional[int] = None):
    """Per-pixel ``log2(product / reactant)`` series (None where undefined)."""
    rc, rv = lipid_spatial_series(dataset, reactant_name, z=z)
    pc, pv = lipid_spatial_series(dataset, product_name, z=z)
    r_by_coord = {c: v for c, v in zip(rc, rv)}
    coords: List[Coord] = []
    vals: List[Optional[float]] = []
    for c, p in zip(pc, pv):
        r = r_by_coord.get(c)
        coords.append(c)
        if r and r > 0 and p is not None and p > 0:
            vals.append(float(np.log2(p / r)))
        else:
            vals.append(None)
    return coords, vals


def _samples_in_group(dataset, group: str) -> set:
    return {s.sample_name for s in dataset.samples if s.group == group and s.sample_name}


def _mean_log2_ratio(dataset, reactant_name, product_name, sample_names) -> Optional[float]:
    """Mean per-pixel log2(product/reactant) over the given samples, or None."""
    r = next((lp for lp in dataset.lipids if lp.input_name == reactant_name), None)
    p = next((lp for lp in dataset.lipids if lp.input_name == product_name), None)
    if r is None or p is None:
        return None
    ratios = []
    for sn in sample_names:
        rv = r.values.get(sn)
        pv = p.values.get(sn)
        if rv and rv > 0 and pv and pv > 0:
            ratios.append(float(np.log2(pv / rv)))
    return float(np.mean(ratios)) if ratios else None


def reaction_effect_sizes(dataset, group_a: str, group_b: str) -> List[dict]:
    """Per-reaction effect size between two stacks: Δ of mean log2(product/reactant).

    A magnitude/direction measure that (unlike the pixel-level t-test z-score) is not
    inflated by pixel count. One row per reaction with a measured reactant+product ion.
    """
    a_names = _samples_in_group(dataset, group_a)
    b_names = _samples_in_group(dataset, group_b)
    out = []
    for rx, react_ions, prod_ions in spatial_reactions(dataset):
        r, p = react_ions[0], prod_ions[0]
        ma = _mean_log2_ratio(dataset, r, p, a_names)
        mb = _mean_log2_ratio(dataset, r, p, b_names)
        delta = (mb - ma) if (ma is not None and mb is not None) else None
        out.append({
            "reaction": rx.reaction_name or f"reaction {getattr(rx, 'reaction_id', '')}",
            "reactant": r, "product": p,
            f"mean_log2ratio_{group_a}": ma, f"mean_log2ratio_{group_b}": mb,
            "effect_size_delta": delta,
            "genes": ", ".join(reaction_gene_labels(rx)),
        })
    out.sort(key=lambda d: abs(d["effect_size_delta"]) if d["effect_size_delta"] is not None else -1,
             reverse=True)
    return out


def reaction_gene_labels(reaction, limit: int = 12) -> List[str]:
    """Gene symbols (+ EC numbers) involved in a reaction, from its genes/proteins."""
    out, seen = [], set()
    for g in getattr(reaction, "genes", None) or []:
        name = g.get("gene_name") or g.get("gene") or g.get("symbol") or g.get("uniprot_id")
        if name and name not in seen:
            seen.add(name)
            out.append(str(name))
    for pr in getattr(reaction, "proteins", None) or []:
        ec = pr.get("ec_number") or pr.get("ec")
        label = f"EC {ec}" if ec else None
        if label and label not in seen:
            seen.add(label)
            out.append(label)
    return out[:limit]


def prominent_reactions(dataset, top: Optional[int] = None) -> List[tuple]:
    """Reactions ranked by spatial prominence = total reactant+product signal.

    Returns ``(reaction, reactant_ions, product_ions, score)`` tuples, most prominent
    first — used to pick which reactions get a 3D render.
    """
    by_name = {lp.input_name: lp for lp in dataset.lipids}
    scored = []
    for rx, react_ions, prod_ions in spatial_reactions(dataset):
        r = by_name.get(react_ions[0])
        p = by_name.get(prod_ions[0])
        score = 0.0
        for lp in (r, p):
            if lp is not None:
                score += float(sum(v for v in lp.values.values() if v))
        scored.append((rx, react_ions, prod_ions, score))
    scored.sort(key=lambda t: t[3], reverse=True)
    return scored[:top] if top else scored


def reaction_ratio_cloud(dataset, reactant_name, product_name, group: Optional[str] = None):
    """3D point cloud of per-pixel log2(product/reactant), optionally within one group.

    Returns ``(coords, values)`` over every spatial pixel (of ``group`` if given);
    value is ``log2(product/reactant)`` where both are positive, else ``None``.
    """
    r = next((lp for lp in dataset.lipids if lp.input_name == reactant_name), None)
    p = next((lp for lp in dataset.lipids if lp.input_name == product_name), None)
    coords: List[Coord] = []
    vals: List[Optional[float]] = []
    if r is None or p is None:
        return coords, vals
    for s in dataset.samples:
        if s.coordinates is None or (group is not None and s.group != group):
            continue
        rv = r.values.get(s.sample_name)
        pv = p.values.get(s.sample_name)
        coords.append((s.coordinates.x, s.coordinates.y, s.coordinates.z))
        vals.append(float(np.log2(pv / rv)) if (rv and rv > 0 and pv and pv > 0) else None)
    return coords, vals


def normalize_per_section(coords, values):
    """Min-max normalize values to [0,1] **within each z-section** (None stays None).

    So each slice is coloured relative to its own points, revealing within-section
    structure instead of being dominated by section-to-section scale differences.
    """
    by_z: dict = {}
    for i, c in enumerate(coords):
        z = c[2] if len(c) > 2 else 0
        by_z.setdefault(z, []).append(i)
    out = [None] * len(values)
    for idxs in by_z.values():
        present = [values[i] for i in idxs if values[i] is not None]
        if not present:
            continue
        lo, hi = min(present), max(present)
        span = (hi - lo) or 1.0
        for i in idxs:
            v = values[i]
            out[i] = None if v is None else (v - lo) / span
    return out


def _extent(coords):
    xs = [p[0] for p in coords]
    ys = [p[1] for p in coords]
    return (min(xs), max(xs), min(ys), max(ys)) if coords else (0, 0, 0, 0)


def _tile_of(p, ext, n):
    x0, x1, y0, y1 = ext
    wx = ((x1 - x0 + 1) / n) or 1.0
    wy = ((y1 - y0 + 1) / n) or 1.0
    col = min(n - 1, int((p[0] - x0) / wx))
    row = min(n - 1, int((p[1] - y0) / wy))
    z = p[2] if len(p) > 2 else 0
    return (z, row, col)


def tile_aligned_delta(dataset, reactant_name, product_name, group_a, group_b, n: int = 4) -> dict:
    """Compare two stacks by **aligned tiles** (valid even for independent pixel grids).

    Each stack's (x, y) extent is partitioned into an n×n grid per z-slice; for each
    aligned ``(z, row, col)`` tile we take the mean log2(product/reactant) in each stack
    and their difference ``Δ = mean_b − mean_a``. Returns
    ``{"delta": {(z,row,col): Δ}, "n": n}`` — a tile-grid comparison for one cube.
    """
    ca, va = reaction_ratio_cloud(dataset, reactant_name, product_name, group=group_a)
    cb, vb = reaction_ratio_cloud(dataset, reactant_name, product_name, group=group_b)
    exa, exb = _extent(ca), _extent(cb)

    def tile_means(coords, vals, ext):
        acc: dict = {}
        for p, v in zip(coords, vals):
            if v is None:
                continue
            acc.setdefault(_tile_of(p, ext, n), []).append(v)
        return {k: float(np.mean(v)) for k, v in acc.items()}

    ma, mb = tile_means(ca, va, exa), tile_means(cb, vb, exb)
    delta = {k: (mb[k] - ma[k]) for k in (set(ma) & set(mb))}
    return {"delta": delta, "n": n, "mean_a": ma, "mean_b": mb}


def aggregate_replicates(dataset, by: str = "section", tiles: int = 16):
    """Collapse pixels into replicate units so stats treat each unit (not each pixel)
    as a replicate — countering pixel pseudo-replication in the between-stack t-test.

    ``by``:
      - ``"pixel"`` — no change (each pixel is a replicate; returns the dataset as-is).
      - ``"section"`` — one unit per (group, z-slice); good for multi-section 3D stacks.
      - ``"tile"`` — partition each group's (x, y) extent into a ~√tiles × √tiles grid;
        each occupied tile is a unit (gives real replicates even for one 2D section).

    Returns a NEW ``LipidDataset`` whose samples are the units (group preserved,
    ``values`` = mean ion intensity over the unit's pixels), carrying the original
    ``reactions`` so the BioPAN comparison works unchanged.
    """
    import math
    from ..models.sample import LipidDataset, QuantifiedLipid, SampleMetadata

    by = (by or "pixel").lower()
    if by == "pixel":
        return dataset

    per_group: dict = {}
    for s in dataset.samples:
        if s.coordinates is not None:
            per_group.setdefault(s.group, []).append(s)

    unit_of: dict = {}
    for group, samples in per_group.items():
        if by == "section":
            for s in samples:
                unit_of[s.sample_name] = f"{group}#z{s.coordinates.z}"
        elif by == "tile":
            xs = [s.coordinates.x for s in samples]
            ys = [s.coordinates.y for s in samples]
            x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
            n = max(1, int(math.ceil(math.sqrt(max(1, tiles)))))
            wx = (x1 - x0 + 1) / n
            wy = (y1 - y0 + 1) / n
            for s in samples:
                cx = min(n - 1, int((s.coordinates.x - x0) / wx)) if wx > 0 else 0
                cy = min(n - 1, int((s.coordinates.y - y0) / wy)) if wy > 0 else 0
                unit_of[s.sample_name] = f"{group}#t{cy}_{cx}"
        else:
            raise ValueError(f"unknown replicate unit: {by!r} (use pixel/section/tile)")

    unit_group: dict = {}
    unit_members: dict = {}
    for s in dataset.samples:
        u = unit_of.get(s.sample_name)
        if u is None:
            continue
        unit_group[u] = s.group
        unit_members.setdefault(u, []).append(s.sample_name)
    units = sorted(unit_group)

    new_samples = [SampleMetadata(sample_name=u, group=unit_group[u]) for u in units]
    new_lipids = []
    for lp in dataset.lipids:
        vals = {}
        for u in units:
            xs = [lp.values.get(m) for m in unit_members[u]]
            xs = [v for v in xs if v is not None]
            if xs:
                vals[u] = float(np.mean(xs))
        new_lipids.append(QuantifiedLipid(
            input_name=lp.input_name, lm_id=lp.lm_id, generic_lm_id=lp.generic_lm_id,
            standardized_name=lp.standardized_name, reactions=lp.reactions,
            annotation_candidates=lp.annotation_candidates, values=vals,
        ))
    return LipidDataset(samples=new_samples, lipids=new_lipids, reactions=dataset.reactions)


def compare_stack_reactions(dataset, group_a: str, group_b: str,
                            threshold: float = 0.05, level: str = "species") -> dict:
    """Compare reactions between two stacks (groups) via BioPAN's z-score machinery.

    ``group_a`` is the control, ``group_b`` the condition of interest. Returns
    ``{"rows": [...], "graph": {nodes, edges}}`` where each row is
    ``{reaction, source, target, z_score, direction, significant}`` ranked by |z|.
    Pure reuse of ``BioPANPathwayExporter.build_reaction_graph`` — the per-reaction
    z-score is a t-test on per-sample (per-pixel) substrate→product ratios between the
    two groups. See :func:`reaction_effect_sizes` for a pixel-count-robust magnitude.
    """
    from scipy.stats import norm
    from ..biopan_pathway_exporter import BioPANPathwayExporter

    exporter = BioPANPathwayExporter(dataset=dataset)
    graph = exporter.build_reaction_graph(
        disease_group=group_b, control_group=group_a, level=level,
    )
    labels = {n["data"]["id"]: (n["data"].get("label") or n["data"].get("name") or n["data"]["id"])
              for n in graph.get("nodes", [])}
    crit = float(norm.ppf(1.0 - threshold))
    rows = []
    for edge in graph.get("edges", []):
        d = edge.get("data", {})
        z = d.get("weight")
        if z is None:
            continue
        src, tgt = labels.get(d.get("source"), d.get("source")), labels.get(d.get("target"), d.get("target"))
        rows.append({
            "reaction": f"{src} → {tgt}",
            "source": src, "target": tgt, "z_score": z,
            "direction": f"up in {group_b}" if z > 0 else f"up in {group_a}",
            "significant": abs(z) >= crit,
        })
    rows.sort(key=lambda r: abs(r["z_score"]), reverse=True)
    return {"rows": rows, "graph": graph}


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
