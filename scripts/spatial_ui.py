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
    ion_display_full,
    ion_display_name,
    lipid_spatial_series,
    ratio_series,
    spatial_reactions,
    z_slices,
)
# Plotly figure builders live in the package so the Streamlit app and the CLI's HTML
# output render identically (rich hover, z labels, µm scale bar).
from lipidmaps.data.utils.spatial_html import grid_figure, scatter3d_figure, voronoi_figure

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
    dims = (f"3D: {dataset.z_slice_count} z-slices" if dataset.is_3d
            else "2D (single z-slice)")
    st.caption(
        f"{dims} · {len(dataset.spatial_samples())} pixels · {len(dataset.lipids)} annotated ions"
    )
    groups = sorted({s.group for s in dataset.samples if s.group})
    if len(groups) > 1:
        st.info(
            f"Multiple stacks/groups detected ({', '.join(groups)}). Open the **BioPAN** tab "
            "to compare reactions between them (ranked z-score table + network)."
        )

    options = _ion_options(dataset)
    if not options:
        st.info("No annotated ions with intensity were found in this dataset.")
        return

    # Show the friendliest name (resolved molecule name when known), keep input_name
    # as the value so lookups still work.
    by_name = {lp.input_name: lp for lp in dataset.lipids}
    label_map = {n: ion_display_full(by_name[n]) for n in options}
    lipid_name = st.selectbox(
        "Lipid / ion", options, key=f"{tab_key_prefix}_lipid",
        format_func=lambda n: label_map.get(n, n),
    )
    selected = by_name.get(lipid_name)
    title = ion_display_full(selected) if selected else lipid_name

    # For auto-annotated (formula+adduct) ions, surface the candidate molecule names.
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

    modes = ["Ion image (grid)", "Voronoi (irregular pixels)"]
    if dataset.is_3d:
        modes.append("3D volume")
    mode = st.radio("Render mode", modes, horizontal=True, key=f"{tab_key_prefix}_mode")
    ps = getattr(dataset, "pixel_size_um", None)
    if ps:
        st.caption(f"Pixel size: {ps[0]:g} × {ps[1]:g} µm · hover a tile for pixel (x,y,z) + value")
    else:
        st.caption("Hover a tile for pixel (x, y, z) + value (no pixel size in file metadata)")

    if mode == "3D volume":
        coords_all, values_all = lipid_spatial_series(dataset, lipid_name)  # all z
        fig = scatter3d_figure(coords_all, values_all, f"{title} (3D)",
                               z_spacing_um=getattr(dataset, "z_spacing_um", None))
        if fig is None:
            st.caption("no values to plot in 3D")
        else:
            st.plotly_chart(fig, use_container_width=True, key=f"{tab_key_prefix}_single3d")
    else:
        _render_series(
            *lipid_spatial_series(dataset, lipid_name, z=z), z, mode.startswith("Ion"),
            title, key=f"{tab_key_prefix}_single", pixel_size_um=ps,
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

    def _lbl(name):
        lp = by_name.get(name)
        return ion_display_name(lp) if lp else name

    st.caption(f"{_lbl(reactant_name)}  →  {_lbl(product_name)}")
    # The reaction maps use grid unless the user explicitly picked Voronoi.
    rxn_grid = not mode.startswith("Voronoi")

    c_react, c_prod, c_ratio = st.columns(3)
    with c_react:
        st.markdown("**Reactant**")
        _render_series(
            *lipid_spatial_series(dataset, reactant_name, z=z), z, rxn_grid,
            _lbl(reactant_name), key=f"{tab_key_prefix}_rxn_react_map", pixel_size_um=ps,
        )
    with c_prod:
        st.markdown("**Product**")
        _render_series(
            *lipid_spatial_series(dataset, product_name, z=z), z, rxn_grid,
            _lbl(product_name), key=f"{tab_key_prefix}_rxn_prod_map", pixel_size_um=ps,
        )
    with c_ratio:
        st.markdown("**log2(product / reactant)**")
        rc, rvals = ratio_series(dataset, reactant_name, product_name, z)
        _render_series(
            rc, rvals, z, rxn_grid, "reaction activity",
            key=f"{tab_key_prefix}_rxn_ratio_map", cmap="RdBu_r", midpoint=0.0,
            pixel_size_um=ps,
        )
