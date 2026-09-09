"""Streamlit UI for spatial (MSI / imzML) lipidomics exploration.

Renders per-lipid ion images (regular grid) and Voronoi maps (irregular pixels),
plus a minimal reaction overlay showing reactant-vs-product ion maps side by side.
Kept thin: all spatial math lives in ``lipidmaps.data.utils.spatial``.
"""

from __future__ import annotations

import logging
from typing import List

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from plotly.colors import sample_colorscale

from lipidmaps.data.utils.spatial import (
    lipid_spatial_series,
    values_to_grid,
    voronoi_regions,
    z_slices,
)

logger = logging.getLogger(__name__)


def _ion_options(dataset) -> List[str]:
    """Ion names that have at least one non-null pixel value."""
    return [
        lp.input_name
        for lp in dataset.lipids
        if any(v is not None for v in lp.values.values())
    ]


def _ion_image_figure(dataset, lipid_name: str, z: int):
    coords, values = lipid_spatial_series(dataset, lipid_name, z=z)
    grid, _mask, _xr, _yr = values_to_grid(coords, values, z=z)
    fig = px.imshow(
        grid,
        origin="lower",
        color_continuous_scale="Viridis",
        aspect="equal",
        labels={"color": "intensity", "x": "x", "y": "y"},
        title=f"{lipid_name} (z={z})",
    )
    fig.update_layout(margin=dict(l=10, r=10, t=40, b=10))
    return fig


def _voronoi_figure(dataset, lipid_name: str, z: int):
    coords, values = lipid_spatial_series(dataset, lipid_name, z=z)
    regions = voronoi_regions(coords, values, z=z)
    if not regions:
        return None

    vals = np.array(
        [r["value"] if r["value"] is not None else np.nan for r in regions], dtype=float
    )
    finite = vals[np.isfinite(vals)]
    vmin, vmax = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 1.0)
    span = (vmax - vmin) or 1.0
    norm = [0.0 if not np.isfinite(v) else max(0.0, min(1.0, (v - vmin) / span)) for v in vals]
    colors = sample_colorscale("Viridis", norm)

    fig = go.Figure()
    for region, color in zip(regions, colors):
        poly = region["polygon"]
        if poly.shape[0] < 3:
            continue
        xs = list(poly[:, 0]) + [poly[0, 0]]
        ys = list(poly[:, 1]) + [poly[0, 1]]
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                fill="toself",
                mode="lines",
                fillcolor=color,
                line=dict(width=0.3, color="rgba(0,0,0,0.25)"),
                hoverinfo="text",
                text=f"({region['x']},{region['y']}) = {region['value']}",
                showlegend=False,
            )
        )
    fig.update_layout(
        title=f"{lipid_name} — Voronoi (z={z})",
        xaxis_title="x",
        yaxis_title="y",
        margin=dict(l=10, r=10, t=40, b=10),
    )
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


def _ion_names_for_component_ids(dataset, lm_ids: set) -> List[str]:
    """Ion names in the dataset whose lm_id/generic_lm_id is in ``lm_ids``."""
    if not lm_ids:
        return []
    names = []
    for lp in dataset.lipids:
        if lp.lm_id in lm_ids or lp.generic_lm_id in lm_ids:
            if any(v is not None for v in lp.values.values()):
                names.append(lp.input_name)
    return names


def _component_ids(components) -> set:
    ids = set()
    for comp in components or []:
        for attr in ("compound_lm_id", "compound_generic_lm_id"):
            val = getattr(comp, attr, None)
            if val:
                ids.add(val)
    return ids


def _render_side(dataset, title: str, ion_names: List[str], z: int, key: str):
    st.markdown(f"**{title}**")
    if not ion_names:
        st.caption("not measured in this dataset")
        return
    name = ion_names[0] if len(ion_names) == 1 else st.selectbox(
        title, ion_names, key=key
    )
    st.plotly_chart(_ion_image_figure(dataset, name, z), use_container_width=True, key=key + "_img")


def render_spatial_explorer(dataset, tab_key_prefix: str = "spatial") -> None:
    st.subheader("Spatial ion maps")
    st.caption(
        f"{len(dataset.spatial_samples())} pixels · {len(dataset.lipids)} annotated ions"
    )

    options = _ion_options(dataset)
    if not options:
        st.info("No annotated ions with intensity were found in this dataset.")
        return

    lipid_name = st.selectbox("Lipid / ion", options, key=f"{tab_key_prefix}_lipid")

    coords, _ = lipid_spatial_series(dataset, lipid_name)
    slices = z_slices(coords) or [0]
    z = slices[0] if len(slices) == 1 else st.selectbox(
        "z-slice", slices, key=f"{tab_key_prefix}_z"
    )

    mode = st.radio(
        "Render mode",
        ["Ion image (grid)", "Voronoi (irregular pixels)"],
        horizontal=True,
        key=f"{tab_key_prefix}_mode",
    )
    if mode.startswith("Ion"):
        st.plotly_chart(
            _ion_image_figure(dataset, lipid_name, z),
            use_container_width=True,
            key=f"{tab_key_prefix}_img",
        )
    else:
        fig = _voronoi_figure(dataset, lipid_name, z)
        if fig is None:
            st.info("Too few pixels for a Voronoi tessellation.")
        else:
            st.plotly_chart(fig, use_container_width=True, key=f"{tab_key_prefix}_vor")

    # --- Reaction overlay (v1: reactant vs product ion maps side by side) ---
    st.markdown("---")
    st.subheader("Reaction overlay")
    lipid_obj = next((lp for lp in dataset.lipids if lp.input_name == lipid_name), None)
    reactions = getattr(lipid_obj, "reactions", None) or []
    if not reactions:
        st.caption(
            "No reactions attached to this lipid. Enable 'fetch reactions' when "
            "processing to see reactant/product spatial comparisons."
        )
        return

    labels = [r.reaction_name or f"reaction {r.reaction_id}" for r in reactions]
    picked = st.selectbox("Reaction", labels, key=f"{tab_key_prefix}_rx")
    rx = reactions[labels.index(picked)]

    reactant_ions = _ion_names_for_component_ids(dataset, _component_ids(rx.reactants))
    product_ions = _ion_names_for_component_ids(dataset, _component_ids(rx.products))

    left, right = st.columns(2)
    with left:
        _render_side(dataset, "Reactants", reactant_ions, z, key=f"{tab_key_prefix}_rxn_reactant")
    with right:
        _render_side(dataset, "Products", product_ions, z, key=f"{tab_key_prefix}_rxn_product")
