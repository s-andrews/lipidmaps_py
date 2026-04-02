import json

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from lipidmaps.data.biopan_pathway_exporter import BioPANPathwayExporter


def _dataset_groups(dataset):
    groups = []
    for sample in getattr(dataset, "samples", []) or []:
        group_name = getattr(sample, "group", None) or "Unknown"
        if group_name not in groups:
            groups.append(group_name)
    return groups


def _flatten_biopan_table(table_payload):
    rows = []
    for row in table_payload.get("pathways", []) or []:
        data = row.get("data", {})
        rows.append(
            {
                "transition": (data.get("pathway") or "").replace("&#8594;", " -> "),
                "score": data.get("score"),
                "gene": data.get("gene") or "",
                "nodes": len(row.get("node_ids", []) or []),
                "edges": len(row.get("edge_ids", []) or []),
            }
        )
    return pd.DataFrame(rows)


def _highlight_ids(highlight_payload, key):
    values = []
    for raw_value in (highlight_payload.get(key, "") or "").split(","):
        cleaned = raw_value.strip().lstrip("#")
        if cleaned:
            values.append(cleaned)
    return values


def _render_cytoscape_graph(graph_payload, highlight_payload, component_id, height=540):
    if not graph_payload.get("nodes") and not graph_payload.get("edges"):
        st.info("No graph elements available for the current selection.")
        return

    payload = {
        "nodes": graph_payload.get("nodes", []),
        "edges": graph_payload.get("edges", []),
        "highlight_nodes": _highlight_ids(highlight_payload, "nodes"),
        "highlight_edges": _highlight_ids(highlight_payload, "edges"),
    }
    graph_dom_id = f"cy_{component_id}"
    html = f"""
    <div style=\"font-family: sans-serif;\">
        <div style=\"display:flex; gap:16px; flex-wrap:wrap; margin-bottom:10px; font-size:13px;\">
            <span><strong>Graph</strong>: Cytoscape.js layout</span>
            <span><span style=\"display:inline-block;width:10px;height:10px;background:#24a19c;border-radius:50%;margin-right:6px;\"></span>Active</span>
            <span><span style=\"display:inline-block;width:10px;height:10px;background:#A634C7;border-radius:50%;margin-right:6px;\"></span>Suppressed</span>
            <span><span style=\"display:inline-block;width:10px;height:10px;border:2px solid #ff8c42;border-radius:50%;margin-right:6px;\"></span>Highlighted by threshold</span>
        </div>
        <div id=\"{graph_dom_id}\" style=\"width:100%; height:{height}px; border:1px solid #d9d9d9; border-radius:10px; background:#ffffff;\"></div>
    </div>
    <script src=\"https://unpkg.com/cytoscape@3.30.4/dist/cytoscape.min.js\"></script>
    <script>
        const payload = {json.dumps(payload)};
        const graphElement = document.getElementById({json.dumps(graph_dom_id)});

        if (typeof cytoscape === 'undefined') {{
            graphElement.innerHTML = '<div style="padding:16px;font-family:sans-serif;color:#666;">Cytoscape.js could not be loaded. Check network access for the embedded graph.</div>';
        }} else {{
            const highlightNodes = new Set(payload.highlight_nodes || []);
            const highlightEdges = new Set(payload.highlight_edges || []);
            const elements = [...(payload.nodes || []), ...(payload.edges || [])];

            const cy = cytoscape({{
                container: graphElement,
                elements,
                layout: {{ name: 'cose', animate: false, fit: true, padding: 24, nodeRepulsion: 9000, idealEdgeLength: 120 }},
                style: [
                    {{
                        selector: 'node',
                        style: {{
                            'background-color': '#f6f3ea',
                            'border-color': '#333333',
                            'border-width': 1.2,
                            'shape': 'ellipse',
                            'label': 'data(label)',
                            'font-size': 12,
                            'text-wrap': 'wrap',
                            'text-max-width': 110,
                            'text-valign': 'center',
                            'color': '#222222',
                            'width': 54,
                            'height': 54
                        }}
                    }},
                    {{
                        selector: 'edge',
                        style: {{
                            'curve-style': 'bezier',
                            'target-arrow-shape': 'triangle',
                            'arrow-scale': 1.1,
                            'line-color': 'data(color)',
                            'target-arrow-color': 'data(color)',
                            'width': 3,
                            'opacity': 0.8
                        }}
                    }},
                    {{
                        selector: 'node.highlighted',
                        style: {{
                            'border-color': '#ff8c42',
                            'border-width': 3,
                            'opacity': 1
                        }}
                    }},
                    {{
                        selector: 'edge.highlighted',
                        style: {{
                            'line-color': '#ff8c42',
                            'target-arrow-color': '#ff8c42',
                            'width': 5,
                            'opacity': 1
                        }}
                    }}
                ]
            }});

            cy.nodes().forEach((node) => {{
                if (highlightNodes.has(node.id())) {{
                    node.addClass('highlighted');
                }}
            }});
            cy.edges().forEach((edge) => {{
                if (highlightEdges.has(edge.id())) {{
                    edge.addClass('highlighted');
                }}
            }});
        }}
    </script>
    """
    components.html(html, height=height + 60, scrolling=False)


def render_biopan_explorer(dataset, tab_key_prefix="biopan"):
    groups = _dataset_groups(dataset)
    hero = st.container(border=True)
    hero.subheader("Reaction Explorer")
    hero.caption("Compare sample groups, inspect ranked reaction/pathway tables, and view Cytoscape reaction graphs generated from the lipidmaps_py BioPAN exporter.")

    if len(groups) < 2:
        hero.info("This view needs at least two sample groups in the processed dataset before it can render BioPAN comparisons.")
        return

    control_col, disease_col, threshold_col, paired_col = hero.columns(4)
    with control_col:
        control_group = st.selectbox(
            "Control group",
            groups,
            index=0,
            key=f"{tab_key_prefix}_control_group",
        )
    disease_options = [group for group in groups if group != control_group] or groups
    with disease_col:
        disease_group = st.selectbox(
            "Disease group",
            disease_options,
            index=0,
            key=f"{tab_key_prefix}_disease_group",
        )
    with threshold_col:
        biopan_threshold = st.number_input(
            "P-value threshold",
            min_value=0.0001,
            max_value=0.5,
            value=0.05,
            step=0.01,
            format="%.4f",
            key=f"{tab_key_prefix}_threshold",
        )
    with paired_col:
        paired = st.checkbox("Paired comparison", value=False, key=f"{tab_key_prefix}_paired")

    graph_family_col, graph_mode_col, table_family_col, table_mode_col, level_col = hero.columns(5)
    with graph_family_col:
        graph_family = st.selectbox(
            "Graph family",
            ["reaction", "pathway"],
            key=f"{tab_key_prefix}_graph_family",
        )
    with graph_mode_col:
        graph_mode = st.selectbox(
            "Graph mode",
            ["active", "suppressed"],
            key=f"{tab_key_prefix}_graph_mode",
        )
    with table_family_col:
        table_family = st.selectbox(
            "Table family",
            ["reaction", "pathway"],
            key=f"{tab_key_prefix}_table_family",
        )
    with table_mode_col:
        table_mode = st.selectbox(
            "Table mode",
            ["active", "suppressed", "most_active", "most_suppressed"],
            key=f"{tab_key_prefix}_table_mode",
        )
    with level_col:
        level = st.selectbox(
            "Level",
            ["class", "species"],
            key=f"{tab_key_prefix}_level",
        )

    with st.spinner("Building BioPAN reaction tables and graph..."):
        exporter = BioPANPathwayExporter(dataset=dataset)
        result_set, reaction_lookup = exporter.build_reaction_match_set(dataset)
        if graph_family == "reaction":
            graph_payload = exporter.build_reaction_graph(
                disease_group=disease_group,
                control_group=control_group,
                level=level,
                mode=graph_mode,
                paired=paired,
                dataset=dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
            )
            graph_highlight = exporter.build_reaction_highlight(
                disease_group=disease_group,
                control_group=control_group,
                threshold=biopan_threshold,
                level=level,
                mode=graph_mode,
                paired=paired,
                dataset=dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
            )
        else:
            graph_payload = exporter.build_pathway_graph(
                disease_group=disease_group,
                control_group=control_group,
                level=level,
                mode=graph_mode,
                paired=paired,
                dataset=dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
            )
            graph_highlight = exporter.build_pathway_highlight(
                disease_group=disease_group,
                control_group=control_group,
                threshold=biopan_threshold,
                level=level,
                mode=graph_mode,
                paired=paired,
                dataset=dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
            )

        if table_family == "reaction":
            biopan_table = exporter.build_reaction_table(
                disease_group=disease_group,
                control_group=control_group,
                threshold=biopan_threshold,
                level=level,
                mode=table_mode,
                paired=paired,
                dataset=dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
            )
        else:
            biopan_table = exporter.build_pathway_table(
                disease_group=disease_group,
                control_group=control_group,
                threshold=biopan_threshold,
                level=level,
                mode=table_mode,
                paired=paired,
                dataset=dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
            )

    graph_stats_col, edge_stats_col, table_stats_col = hero.columns(3)
    with graph_stats_col:
        st.metric("Graph nodes", len(graph_payload.get("nodes", [])))
    with edge_stats_col:
        st.metric("Graph edges", len(graph_payload.get("edges", [])))
    with table_stats_col:
        st.metric("Ranked rows", len(biopan_table.get("pathways", []) or []))

    hero.caption(
        f"Comparison: {disease_group} vs {control_group} | level={level} | graph={graph_family}:{graph_mode} | table={table_family}:{table_mode}"
    )
    _render_cytoscape_graph(
        graph_payload,
        graph_highlight,
        component_id=f"{tab_key_prefix}_{graph_family}_{graph_mode}_{level}_{disease_group}_{control_group}",
    )

    table_df = _flatten_biopan_table(biopan_table)
    if table_df.empty:
        st.info("No BioPAN rows passed the current threshold for this selection.")
    else:
        st.dataframe(table_df, hide_index=True, use_container_width=True)
        csv_bytes = table_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download BioPAN table CSV",
            data=csv_bytes,
            file_name=f"biopan_{table_family}_{table_mode}_{level}_{disease_group}_vs_{control_group}.csv",
            mime="text/csv",
            key=f"{tab_key_prefix}_download_csv",
        )