"""Static (matplotlib) renderers for MSI spatial maps -- used by the CLI.

Keeps Plotly confined to the Streamlit UI; these write PNGs for headless/offline
runs on large datasets. matplotlib is an optional (``[msi]``) dependency, imported
lazily with a non-interactive backend so importing this module never requires a
display and the core package stays matplotlib-free.

All functions reuse the pure helpers in :mod:`lipidmaps.data.utils.spatial`.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from .spatial import Coord, nice_length, values_to_grid, voronoi_regions

logger = logging.getLogger(__name__)

_XLABEL = "pixel x (col)"
_YLABEL = "pixel y (row)"


def _title_z(title: str, z: int) -> str:
    return title if f"z={z}" in title else f"{title} (z={z})"


def _scale_bar(ax, x_min, x_max, y_min, y_max, pixel_size_um) -> None:
    """Draw a µm scale bar (bottom-left) when the pixel size is known."""
    if not pixel_size_um:
        return
    psx = pixel_size_um[0]
    if not psx or psx <= 0:
        return
    span = max(1.0, x_max - x_min)
    length_um = nice_length(span * psx * 0.25)
    length_px = length_um / psx
    x0 = x_min + 0.04 * span
    y0 = y_min + 0.06 * max(1.0, (y_max - y_min))
    ax.plot([x0, x0 + length_px], [y0, y0], color="white", lw=3, solid_capstyle="butt")
    ax.text(x0 + length_px / 2, y0, f"{length_um:g} µm", color="white",
            ha="center", va="bottom", fontsize=8)


def _plt():
    """Return the pyplot module with the Agg backend, or raise a helpful error."""
    try:
        import matplotlib
        matplotlib.use("Agg")  # headless; no display needed
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - only without the extra
        raise ImportError(
            "Rendering MSI PNGs requires matplotlib (optional). "
            "Install it with: pip install lipidmaps_py[msi]"
        ) from exc
    return plt


def save_ion_grid_png(
    coords: Sequence[Coord],
    values: Sequence[Optional[float]],
    out_path: str | Path,
    title: str,
    z: int = 0,
    cmap: str = "viridis",
    pixel_size_um=None,
) -> Path:
    """Save a grid ion image (one z-slice) as a PNG."""
    plt = _plt()
    grid, _mask, (x_min, x_max), (y_min, y_max) = values_to_grid(coords, values, z=z)
    fig, ax = plt.subplots(figsize=(5, 4))
    x_hi = x_min - 0.5 + grid.shape[1]
    y_hi = y_min - 0.5 + grid.shape[0]
    im = ax.imshow(
        grid, origin="lower", cmap=cmap, aspect="equal",
        extent=[x_min - 0.5, x_hi, y_min - 0.5, y_hi],
    )
    ax.set_title(_title_z(title, z), fontsize=10)
    ax.set_xlabel(_XLABEL)
    ax.set_ylabel(_YLABEL)
    if grid.size:
        _scale_bar(ax, x_min - 0.5, x_hi, y_min - 0.5, y_hi, pixel_size_um)
    fig.colorbar(im, ax=ax, shrink=0.85, label="value")
    return _save(fig, out_path)


def save_voronoi_png(
    coords: Sequence[Coord],
    values: Sequence[Optional[float]],
    out_path: str | Path,
    title: str,
    z: int = 0,
    cmap: str = "viridis",
    midpoint: Optional[float] = None,
    pixel_size_um=None,
) -> Optional[Path]:
    """Save a Voronoi map as a PNG, or return None if too few pixels."""
    plt = _plt()
    from matplotlib.collections import PolyCollection

    regions = voronoi_regions(coords, values, z=z)
    if not regions:
        return None
    polys = [r["polygon"] for r in regions if r["polygon"].shape[0] >= 3]
    vals = np.array(
        [r["value"] if r["value"] is not None else np.nan
         for r in regions if r["polygon"].shape[0] >= 3],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(5, 4))
    coll = PolyCollection(polys, array=vals, cmap=cmap, edgecolors="none")
    if midpoint is not None:
        finite = vals[np.isfinite(vals)]
        extent = float(np.max(np.abs(finite - midpoint))) if finite.size else 1.0
        coll.set_clim(midpoint - extent, midpoint + extent)
    ax.add_collection(coll)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_title(_title_z(title, z), fontsize=10)
    ax.set_xlabel(_XLABEL)
    ax.set_ylabel(_YLABEL)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    _scale_bar(ax, x0, x1, y0, y1, pixel_size_um)
    fig.colorbar(coll, ax=ax, shrink=0.85, label="value")
    return _save(fig, out_path)


def save_ratio_png(
    coords: Sequence[Coord],
    ratio_values: Sequence[Optional[float]],
    out_path: str | Path,
    title: str,
    z: int = 0,
) -> Path:
    """Save a diverging log2(product/reactant) map centered at 0."""
    plt = _plt()
    grid, _mask, (x_min, _x), (y_min, _y) = values_to_grid(coords, ratio_values, z=z)
    finite = grid[np.isfinite(grid)]
    extent = float(np.max(np.abs(finite))) if finite.size else 1.0
    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(
        grid, origin="lower", cmap="RdBu_r", aspect="equal", vmin=-extent, vmax=extent,
        extent=[x_min - 0.5, x_min - 0.5 + grid.shape[1], y_min - 0.5, y_min - 0.5 + grid.shape[0]],
    )
    ax.set_title(_title_z(title, z), fontsize=10)
    ax.set_xlabel(_XLABEL)
    ax.set_ylabel(_YLABEL)
    fig.colorbar(im, ax=ax, shrink=0.85, label="log2(product/reactant)")
    return _save(fig, out_path)


def _save(fig, out_path: str | Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    import matplotlib.pyplot as plt
    plt.close(fig)
    return out_path
