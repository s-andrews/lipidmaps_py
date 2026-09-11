"""Plotly spatial-map builders shared by the Streamlit app and the CLI HTML output.

Returns ``plotly.graph_objects.Figure`` objects with rich hover (pixel x, y, z +
sample_name + value), clear axis labels, z in the title, and an optional µm scale bar
when the pixel size is known. plotly is an optional ([msi]/[dev]) dependency, imported
lazily so the core package never needs it.
"""

from __future__ import annotations

import logging
from typing import Optional, Sequence, Tuple

import numpy as np

from .spatial import Coord, nice_length, values_to_grid, voronoi_regions

logger = logging.getLogger(__name__)


def _plotly():
    import plotly.express as px
    import plotly.graph_objects as go
    return px, go


def _sample_name(x: int, y: int, z: int) -> str:
    return f"px_x{int(x):04d}_y{int(y):04d}_z{int(z):02d}"


def _add_scale_bar(go_fig, x_min, x_max, y_min, y_max, pixel_size_um):
    """Draw a horizontal µm scale bar near the bottom-left, if pixel size is known."""
    if not pixel_size_um:
        return
    psx = pixel_size_um[0]
    if not psx or psx <= 0:
        return
    span_px = max(1.0, x_max - x_min)
    length_um = nice_length(span_px * psx * 0.25)
    length_px = length_um / psx
    x0 = x_min + 0.04 * span_px
    y0 = y_min + 0.06 * max(1.0, (y_max - y_min))
    go_fig.add_shape(type="line", x0=x0, x1=x0 + length_px, y0=y0, y1=y0,
                     line=dict(color="white", width=3))
    label = f"{length_um:g} µm"
    go_fig.add_annotation(x=x0 + length_px / 2, y=y0, text=label, showarrow=False,
                          yshift=10, font=dict(color="white", size=11))


def _title_with_z(title: str, z: int) -> str:
    return title if f"z={z}" in title else f"{title} (z={z})"


def grid_figure(coords: Sequence[Coord], values, z: int, title: str,
                pixel_size_um: Optional[Tuple[float, float]] = None,
                cmap: str = "Viridis", midpoint: Optional[float] = None):
    """Grid ion image (one z-slice) as a plotly figure with x/y/z/value hover."""
    px, go = _plotly()
    grid, _mask, (x_min, x_max), (y_min, y_max) = values_to_grid(coords, values, z=z)
    fig = px.imshow(
        grid, origin="lower", color_continuous_scale=cmap, color_continuous_midpoint=midpoint,
        aspect="equal", labels={"color": "value"}, title=_title_with_z(title, z),
        x=list(range(x_min, x_min + grid.shape[1])) if grid.size else None,
        y=list(range(y_min, y_min + grid.shape[0])) if grid.size else None,
    )
    fig.update_traces(
        hovertemplate="pixel x=%{x}, y=%{y}, z=" + str(z) + "<br>value=%{z}<extra></extra>"
    )
    fig.update_layout(margin=dict(l=10, r=10, t=40, b=10),
                      xaxis_title="pixel x (col)", yaxis_title="pixel y (row)")
    if grid.size:
        _add_scale_bar(fig, x_min, x_max, y_min, y_max, pixel_size_um)
    return fig


def scatter3d_figure(coords: Sequence[Coord], values, title: str,
                     cmap: str = "Viridis", midpoint: Optional[float] = None,
                     z_spacing_um: Optional[float] = None):
    """Volumetric point cloud across all z-slices, coloured by value (plotly Scatter3d).

    Only pixels with a value are plotted. Use for stacked 3D datasets; returns None if
    there is nothing to show.
    """
    _px, go = _plotly()
    xs, ys, zs, vs, texts = [], [], [], [], []
    for c, v in zip(coords, values):
        if v is None:
            continue
        x, y, z = (tuple(c) + (0, 0, 0))[:3]
        xs.append(x)
        ys.append(y)
        zs.append(z * z_spacing_um if z_spacing_um else z)
        vs.append(float(v))
        texts.append(f"pixel ({x}, {y}, {z})<br>{_sample_name(x, y, z)}<br>value={v}")
    if not xs:
        return None
    marker = dict(size=3, color=vs, colorscale=cmap, showscale=True,
                  colorbar=dict(title="value"), opacity=0.8)
    if midpoint is not None:
        extent = max(abs(min(vs) - midpoint), abs(max(vs) - midpoint)) or 1.0
        marker.update(cmin=midpoint - extent, cmax=midpoint + extent)
    fig = go.Figure(go.Scatter3d(x=xs, y=ys, z=zs, mode="markers", marker=marker,
                                 hoverinfo="text", text=texts))
    z_title = "z (µm)" if z_spacing_um else "z-slice"
    fig.update_layout(title=title, margin=dict(l=0, r=0, t=40, b=0),
                      scene=dict(xaxis_title="pixel x", yaxis_title="pixel y", zaxis_title=z_title))
    return fig


def reaction_volume_figure(panels, title: str, gene_text: Optional[str] = None,
                           z_spacing_um: Optional[float] = None, cmap: str = "RdBu_r",
                           midpoint: float = 0.0, color_mode: str = "global"):
    """Tiled 3D reaction figure: one Scatter3d scene per panel, coloured by log2 ratio.

    ``panels`` is a list of ``(panel_title, coords, values)``. One panel → a single 3D
    volume (like the ion 3D view); two panels → stacks tiled side by side.
    ``color_mode``: ``"global"`` colours all points on one diverging scale around
    ``midpoint``; ``"section"`` colours each point **relative to the other points in its
    own z-slice** (min-max within section, Viridis) so within-section structure shows.
    ``gene_text`` (genes involved) is shown as a caption. Returns None if nothing to plot.
    """
    from plotly.subplots import make_subplots
    from .spatial import normalize_per_section

    section = color_mode == "section"
    prepared = []
    for ptitle, coords, values in panels:
        color_vals = normalize_per_section(coords, values) if section else values
        xs, ys, zs, cs, texts = [], [], [], [], []
        for c, v, cv in zip(coords, values, color_vals):
            if v is None or cv is None:
                continue
            x, y, z = (tuple(c) + (0, 0, 0))[:3]
            xs.append(x)
            ys.append(y)
            zs.append(z * z_spacing_um if z_spacing_um else z)
            cs.append(float(cv))
            texts.append(f"pixel ({x}, {y}, {z})<br>log2(prod/react)={v:.3f}")
        prepared.append((ptitle, xs, ys, zs, cs, texts))
    if not any(p[4] for p in prepared):
        return None

    if section:
        cmin, cmax, scale, bar = 0.0, 1.0, "Viridis", "rel. (per section)"
    else:
        allc = [c for p in prepared for c in p[4]]
        extent = max((abs(c - midpoint) for c in allc), default=1.0) or 1.0
        cmin, cmax, scale, bar = midpoint - extent, midpoint + extent, cmap, "log2(p/r)"
    z_title = "z (µm)" if z_spacing_um else "z-slice"

    n = len(prepared)
    fig = make_subplots(rows=1, cols=n, specs=[[{"type": "scene"}] * n],
                        subplot_titles=[p[0] for p in prepared])
    for i, (ptitle, xs, ys, zs, cs, texts) in enumerate(prepared):
        fig.add_trace(_go().Scatter3d(
            x=xs, y=ys, z=zs, mode="markers",
            marker=dict(size=3, color=cs, colorscale=scale, cmin=cmin, cmax=cmax,
                        opacity=0.85, showscale=(i == n - 1), colorbar=dict(title=bar)),
            hoverinfo="text", text=texts, name=ptitle,
        ), row=1, col=i + 1)

    scene = dict(xaxis_title="pixel x", yaxis_title="pixel y", zaxis_title=z_title)
    layout = {"title": title, "margin": dict(l=0, r=0, t=60, b=40)}
    for i in range(n):
        layout["scene" if i == 0 else f"scene{i + 1}"] = scene
    if gene_text:
        layout["annotations"] = [dict(text=f"Genes: {gene_text}", showarrow=False,
                                      xref="paper", yref="paper", x=0, y=-0.02,
                                      xanchor="left", font=dict(size=12))]
    fig.update_layout(**layout)
    return fig


def tile_delta_figure(delta_result: dict, title: str, gene_text: Optional[str] = None,
                      z_spacing_um: Optional[float] = None):
    """One cube of aligned-tile Δ: a marker per (z, tile) coloured by Δ log2(product/reactant).

    Consumes :func:`lipidmaps.data.utils.spatial.tile_aligned_delta` output. Valid for
    independent pixel grids because tiles align by relative grid position, not by pixel.
    """
    go = _go()
    delta = delta_result.get("delta", {})
    if not delta:
        return None
    xs, ys, zs, vs, texts = [], [], [], [], []
    for (z, row, col), d in delta.items():
        xs.append(col)
        ys.append(row)
        zs.append(z * z_spacing_um if z_spacing_um else z)
        vs.append(float(d))
        texts.append(f"tile (row {row}, col {col}) z={z}<br>Δ log2(p/r)={d:+.3f}")
    extent = max((abs(v) for v in vs), default=1.0) or 1.0
    fig = go.Figure(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="markers",
        marker=dict(size=8, color=vs, colorscale="RdBu_r", cmin=-extent, cmax=extent,
                    opacity=0.9, showscale=True, colorbar=dict(title="Δ log2(p/r)")),
        hoverinfo="text", text=texts,
    ))
    layout = {"title": title, "margin": dict(l=0, r=0, t=60, b=40),
              "scene": dict(xaxis_title="tile col", yaxis_title="tile row",
                            zaxis_title="z (µm)" if z_spacing_um else "z-slice")}
    if gene_text:
        layout["annotations"] = [dict(text=f"Genes: {gene_text}", showarrow=False,
                                      xref="paper", yref="paper", x=0, y=-0.02,
                                      xanchor="left", font=dict(size=12))]
    fig.update_layout(**layout)
    return fig


def _go():
    import plotly.graph_objects as go
    return go


def voronoi_figure(coords: Sequence[Coord], values, z: int, title: str,
                   pixel_size_um: Optional[Tuple[float, float]] = None,
                   cmap: str = "Viridis", midpoint: Optional[float] = None):
    """Voronoi map as a plotly figure; each tile hovers pixel x,y,z + sample_name + value."""
    px, go = _plotly()
    from plotly.colors import sample_colorscale

    regions = voronoi_regions(coords, values, z=z)
    if not regions:
        return None

    vals = np.array([r["value"] if r["value"] is not None else np.nan for r in regions], dtype=float)
    finite = vals[np.isfinite(vals)]
    if midpoint is not None:
        extent = float(np.max(np.abs(finite - midpoint))) if finite.size else 1.0
        extent = extent or 1.0
        norm = [0.5 if not np.isfinite(v) else max(0.0, min(1.0, 0.5 + 0.5 * (v - midpoint) / extent)) for v in vals]
    else:
        vmin, vmax = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 1.0)
        span = (vmax - vmin) or 1.0
        norm = [0.0 if not np.isfinite(v) else max(0.0, min(1.0, (v - vmin) / span)) for v in vals]
    colors = sample_colorscale(cmap, norm)

    xs_all, ys_all = [], []
    fig = go.Figure()
    for region, color in zip(regions, colors):
        poly = region["polygon"]
        if poly.shape[0] < 3:
            continue
        rx, ry, rz = region["x"], region["y"], region.get("z", z)
        xs = list(poly[:, 0]) + [poly[0, 0]]
        ys = list(poly[:, 1]) + [poly[0, 1]]
        xs_all += list(poly[:, 0])
        ys_all += list(poly[:, 1])
        fig.add_trace(go.Scatter(
            x=xs, y=ys, fill="toself", mode="lines", fillcolor=color,
            line=dict(width=0.3, color="rgba(0,0,0,0.25)"), hoverinfo="text",
            text=f"pixel ({rx}, {ry}, {rz})<br>{_sample_name(rx, ry, rz)}<br>value={region['value']}",
            showlegend=False,
        ))
    fig.update_layout(title=_title_with_z(title, z), xaxis_title="pixel x (col)",
                      yaxis_title="pixel y (row)", margin=dict(l=10, r=10, t=40, b=10))
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    if xs_all:
        _add_scale_bar(fig, min(xs_all), max(xs_all), min(ys_all), max(ys_all), pixel_size_um)
    return fig
