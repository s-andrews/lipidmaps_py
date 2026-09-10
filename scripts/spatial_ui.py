"""Streamlit UI for spatial (MSI / imzML) lipidomics exploration.

Renders per-lipid ion images (regular grid) and Voronoi maps (irregular pixels),
plus a minimal reaction overlay showing reactant-vs-product ion maps side by side.
Kept thin: all spatial math lives in ``lipidmaps.data.utils.spatial``.
"""

from __future__ import annotations

import logging
from typing import List

import streamlit as st

from lipidmaps.data.utils.spatial import (
    lipid_spatial_series,
    ratio_series,
    spatial_reactions,
    z_slices,
)
# Plotly figure builders live in the package so the Streamlit app and the CLI's HTML
# output render identically (rich hover, z labels, µm scale bar).
from lipidmaps.data.utils.spatial_html import grid_figure, voronoi_figure

logger = logging.getLogger(__name__)


def _ion_options(dataset) -> List[str]:
    """Ion names that have at least one non-null pixel value."""
    return [
        lp.input_name
        for lp in dataset.lipids
        if any(v is not None for v in lp.values.values())
    ]


def _render_series(coords, values, z, mode_is_grid, title, key,
                   cmap="Viridis", midpoint=None, pixel_size_um=None):
    """Render one (coords, values) series as a grid or Voronoi map into this column."""
    if mode_is_grid:
        st.plotly_chart(
            grid_figure(coords, values, z, title, pixel_size_um=pixel_size_um,
                        cmap=cmap, midpoint=midpoint),
            use_container_width=True, key=key,
        )
    else:
        fig = voronoi_figure(coords, values, z, title, pixel_size_um=pixel_size_um,
                             cmap=cmap, midpoint=midpoint)
        if fig is None:
            st.caption("too few pixels for Voronoi")
        else:
            st.plotly_chart(fig, use_container_width=True, key=key)


# Reaction/ion spatial helpers live in the package (lipidmaps.data.utils.spatial) so
# the CLI and this UI share one implementation: spatial_reactions(), ratio_series().


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

    # For auto-annotated (formula+adduct) ions, surface the candidate molecule names.
    selected = next((lp for lp in dataset.lipids if lp.input_name == lipid_name), None)
    candidates = getattr(selected, "annotation_candidates", None) if selected else None
    if candidates:
        lm = getattr(selected, "lm_id", None)
        header = f"{len(candidates)} candidate molecule(s)"
        if lm:
            header += f" · resolved LM ID: {lm} ({selected.standardized_name})"
        with st.expander(header):
            st.write(", ".join(c.get("name", "") for c in candidates[:40]))

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
    mode_is_grid = mode.startswith("Ion")
    ps = getattr(dataset, "pixel_size_um", None)
    if ps:
        st.caption(f"Pixel size: {ps[0]:g} × {ps[1]:g} µm · hover a tile for pixel (x,y,z) + value")
    else:
        st.caption("Hover a tile for pixel (x, y, z) + value (no pixel size in file metadata)")
    _render_series(
        *lipid_spatial_series(dataset, lipid_name, z=z), z, mode_is_grid,
        f"{lipid_name}", key=f"{tab_key_prefix}_single", pixel_size_um=ps,
    )

    # --- Reaction spatial explorer: click a found reaction -> reactant vs product
    # vs product/reactant ratio, mapped over the tissue (grid or Voronoi). ---
    st.markdown("---")
    st.subheader("Reactions over tissue")
    spatial_rx = spatial_reactions(dataset)
    if not spatial_rx:
        st.caption(
            "No found reaction has both a measured reactant and product ion. "
            "Enable 'fetch reactions' when processing an annotated/auto-annotated imzML."
        )
        return

    labels = [
        (rx.reaction_name or f"reaction {rx.reaction_id}") for rx, _, _ in spatial_rx
    ]
    idx = st.selectbox(
        "Found reactions", range(len(labels)),
        format_func=lambda i: labels[i], key=f"{tab_key_prefix}_rx",
    )
    rx, reactant_ions, product_ions = spatial_rx[idx]

    reactant_name = reactant_ions[0] if len(reactant_ions) == 1 else st.selectbox(
        "Reactant ion", reactant_ions, key=f"{tab_key_prefix}_rxn_reactant_sel")
    product_name = product_ions[0] if len(product_ions) == 1 else st.selectbox(
        "Product ion", product_ions, key=f"{tab_key_prefix}_rxn_product_sel")
    st.caption(f"{reactant_name}  →  {product_name}")

    c_react, c_prod, c_ratio = st.columns(3)
    with c_react:
        st.markdown("**Reactant**")
        _render_series(
            *lipid_spatial_series(dataset, reactant_name, z=z), z, mode_is_grid,
            f"{reactant_name}", key=f"{tab_key_prefix}_rxn_react_map", pixel_size_um=ps,
        )
    with c_prod:
        st.markdown("**Product**")
        _render_series(
            *lipid_spatial_series(dataset, product_name, z=z), z, mode_is_grid,
            f"{product_name}", key=f"{tab_key_prefix}_rxn_prod_map", pixel_size_um=ps,
        )
    with c_ratio:
        st.markdown("**log2(product / reactant)**")
        rc, rvals = ratio_series(dataset, reactant_name, product_name, z)
        _render_series(
            rc, rvals, z, mode_is_grid, "reaction activity",
            key=f"{tab_key_prefix}_rxn_ratio_map", cmap="RdBu_r", midpoint=0.0,
            pixel_size_um=ps,
        )
