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
