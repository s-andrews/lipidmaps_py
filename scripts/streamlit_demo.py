import streamlit as st
import os
import sys
import pandas as pd
import plotly.express as px
from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.models import reaction

# ---------------------- BOOTSTRAP ----------------------
dir_path = os.path.dirname(os.path.realpath(__file__))
src_path = os.path.abspath(os.path.join(dir_path, '../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

def main():
    st.set_page_config(
        page_title="LIPID MAPS Quantitative Data Demo",
        page_icon="LM",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    st.header("Quantitative Data Demo")
    
    # Display file being analyzed if processed
    if st.session_state.get("processed") and st.session_state.get("_last_file_used"):
        file_name = st.session_state.get("_last_file_used")
        st.markdown(f"""
        **Analyzing:** `{file_name}`  
        View the results in the tabs below.
        """)
    else:
        st.markdown("""
        Use the sidebar to select a file or upload your CSV.  
        Process the file to see standardized lipid annotations.
        """)

    # ---------------------- SESSION DEFAULTS ----------------------
    defaults = {
        "file_to_use": None,
        "dataset": None,
        "processed": False,
        "validation_issues": [],
        "generic_lm_id_assigned": False,
        "reactions_fetched": False,
        "validation_passed": None,
        "validation_summary": None,
        "has_validation_report": False,
        "show_all_issues": False,
        "show_validation_section": True,
        "reactions": [],              # persistent reactions
        "taxonomy_group": "all",
    }

    for k, v in defaults.items():
        st.session_state.setdefault(k, v)

    # ---------------------- SIDEBAR ----------------------
    processed = False
    generic_lm_id_button = False
    fetch_reactions_button = False

    with st.sidebar:
        st.title("LIPID MAPS API")

        # ---- FILE selection ----
        with st.expander("File", expanded=not bool(st.session_state["processed"])):
            test_data_dir = os.path.abspath(os.path.join(dir_path, '../tests/data/inputs/demo'))
            try:
                test_files = [f for f in os.listdir(test_data_dir)
                            if f.endswith((".tsv", ".csv"))]
            except:
                test_files = []
                st.warning(f"Test data directory not found: {test_data_dir}")

            selected_file = st.selectbox("Select test data", ["(none)"] + test_files)
            uploaded_file = st.file_uploader("Or upload CSV", type=["csv", "tsv"])

            file_to_use = None
            file_name = None
            if uploaded_file:
                import tempfile
                with tempfile.NamedTemporaryFile(delete=False) as tmp:
                    tmp.write(uploaded_file.read())
                    file_to_use = tmp.name
                    file_name = uploaded_file.name
                selected_file = "(none)"
            elif selected_file != "(none)":
                file_to_use = os.path.join(test_data_dir, selected_file)
                file_name = selected_file
                uploaded_file = None
            
            prev_file_name = st.session_state.get("_last_file_used")
            st.session_state["file_to_use"] = file_to_use

            # If the chosen file changed, reset dependent session state to sensible defaults
            if file_name != prev_file_name:
                reset_keys = [
                    "dataset",
                    "processed",
                    "validation_issues",
                    "generic_lm_id_assigned",
                    "reactions_fetched",
                    "validation_passed",
                    "validation_summary",
                    "has_validation_report",
                    "show_all_issues",
                    "show_validation_section",
                    "reactions",
                ]

                for key in reset_keys:
                    if key in defaults:
                        st.session_state[key] = defaults[key]
                    else:
                        # fallback sensible values
                        st.session_state[key] = [] if key.endswith("issues") or key == "reactions" else False

                # clear any UI selection tied to previous dataset
                if "selected_lipid_idx" in st.session_state:
                    st.session_state.pop("selected_lipid_idx", None)

                st.session_state["_last_file_used"] = file_name

        # ---- OPTIONS ----
        with st.expander("Options", expanded=True):
            file_chosen = bool(st.session_state["file_to_use"])
            use_verification = st.checkbox(
                "Use Verification",
                value=True,
                disabled=not file_chosen,
            )
            validate_data = use_verification if file_chosen else False

            use_refmet = st.checkbox(
                "Use Refmet",
                value=True,
                disabled=not file_chosen,
            )

            use_headgroups = st.checkbox(
                "Assign Generic LMIDs from headgroups",
                value=True,
                disabled=not file_chosen,
            )

            fetch_reactions = st.checkbox(
                "Fetch reactions by LM ID",
                value=True,
                disabled=not file_chosen,
            )

            # Taxonomy group selection for reactions fetching
            allowed_taxonomy_groups = ['all', 'bacteria', 'archaea', 'fungi', 'viridiplantae', 'mammalia', 'arthropoda', 'eukaryota']
            # allow choosing taxonomy only once dataset is processed
            st.session_state.setdefault("taxonomy_group", "all")
            taxonomy_group = st.selectbox(
                "Taxonomy group for reactions",
                allowed_taxonomy_groups,
                index=allowed_taxonomy_groups.index(st.session_state.get("taxonomy_group", "all")),
                disabled=not (fetch_reactions and file_chosen),
                help="Restrict reactions lookup to this taxonomy group when fetching reactions. Choose 'all' to omit taxonomy filter.",
            )
            st.session_state["taxonomy_group"] = taxonomy_group

            processed = st.button("Process Data", disabled=not file_chosen)
            if st.session_state["processed"] and file_chosen:
                st.badge("Success! Please use the Processed tab to view results", icon=":material/check:", color="green")

        # # ---- TOOLS ----
        # with st.expander("Tools", expanded=True):
        #     processed_flag = bool(st.session_state["processed"])

        #     generic_lm_id_button = st.button("Assign Generic LMIDs",
        #                                     disabled=(not processed_flag or use_headgroups))
        #     if st.session_state["generic_lm_id_assigned"]:
        #         st.badge("Generic LMIDs assigned", icon=":material/check:", color="green")

        #     fetch_reactions_button = st.button("Fetch reactions by LM ID",
        #                                     disabled=(not processed_flag or fetch_reactions))
        #     if st.session_state["reactions_fetched"]:
        #         st.badge("Reactions fetched", icon=":material/check:", color="green")

        # ---- PAGES ----
        with st.expander("All", expanded=True):
            if st.button("All Reactions"):
                st.session_state["show_all_reactions"] = True
            if st.button("Dashboard"):
                st.session_state["show_all_reactions"] = False

        # ---- VIEW ----
        with st.expander("View", expanded=True):
            st.session_state["show_validation_section"] = st.checkbox(
                "Show Validation Report",
                value=st.session_state["show_validation_section"]
            )

    # ---------------------- TABS ----------------------
    # If user requested All Reactions, load and run `scripts/reactions.py` as main content
    if st.session_state.get("show_all_reactions"):
        try:
            # Try to import the module by path and call its `main()` if present
            import importlib.util
            mod_path = os.path.join(dir_path, "reactions.py")
            if os.path.exists(mod_path):
                spec = importlib.util.spec_from_file_location("scripts.reactions", mod_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)  # type: ignore
                if hasattr(module, "main") and callable(module.main):
                    module.main()
                else:
                    # executed file may include top-level Streamlit code
                    pass
            else:
                st.error(f"Reactions file not found: {mod_path}")
        except Exception as e:
            st.error(f"Failed to load All Reactions page: {e}")
        return
    tab_labels = ["Preview", "Processed", "Reactions", "Validation"]
    tabs = st.tabs(tab_labels)

    tab_index = {name.lower(): i for i, name in enumerate(tab_labels)}

    # --------------------------------------------------------------
    # PREVIEW TAB
    # --------------------------------------------------------------
    with tabs[tab_index["preview"]]:
        if not st.session_state["file_to_use"]:
            st.info("Please select or upload a file to preview.")
        else:
            st.subheader("Preview of Selected File")

            fp = st.session_state["file_to_use"]
            try:
                _, ext = os.path.splitext(fp)
                if ext.lower() in [".tsv", ".txt"]:
                    try:
                        df = pd.read_csv(fp, sep="\t")
                    except pd.errors.ParserError:
                        df = pd.read_csv(fp, sep="\t", skiprows=[1])
                else:
                    try:
                        df = pd.read_csv(fp)
                    except pd.errors.ParserError:
                        df = pd.read_csv(fp, skiprows=[1])

                st.write(f"Rows: {df.shape[0]}, Columns: {df.shape[1]}")
                st.dataframe(df, hide_index=True)
            except Exception as e:
                st.error(f"Error reading file: {e}")

    # --------------------------------------------------------------
    # PROCESSING ACTION
    # --------------------------------------------------------------
    if processed and st.session_state["file_to_use"]:
        try:
            fp = st.session_state["file_to_use"]

            from lipidmaps import process_csv
            dataset = process_csv(fp, validate_data=validate_data, use_refmet=use_refmet, use_headgroups=use_headgroups, taxonomy_group=taxonomy_group)
            st.session_state["dataset"] = dataset
            st.session_state["processed"] = True
            st.session_state["reactions"] = []  # clear old reactions
            st.write(f"validate_data: {validate_data}")
            # Validation report handling
            if validate_data and getattr(dataset, "validation_report", None):
                vr = dataset.validation_report
                st.session_state["validation_passed"] = vr.passed
                st.session_state["validation_issues"] = vr.issues or []
                st.session_state["validation_summary"] = vr.summary or None
                st.session_state["has_validation_report"] = True
            else:
                st.session_state["has_validation_report"] = False
                st.session_state["validation_issues"] = []

            st.rerun()   # important

        except Exception as e:
            st.error(f"Error processing file: {e}")

    # --------------------------------------------------------------
    # PROCESSED TAB — ALWAYS RENDER
    # --------------------------------------------------------------
    with tabs[tab_index["processed"]]:
        dataset = st.session_state.get("dataset")
        processed_table_container = st.container(border=True)
        if not dataset:
            processed_table_container.info("No processed dataset yet.")
        else:
            processed_table_container.subheader("Processed Lipid Annotations")

            # Build DataFrame and keep mapping to lipid objects
            lipid_rows = []
            for idx, lipid in enumerate(dataset.lipids):
                lipid_rows.append({
                    "index": idx,
                    "input_name": getattr(lipid, "input_name", None),
                    "standardized_name": getattr(lipid, "standardized_name", None),
                    "lm_id": getattr(lipid, "lm_id", None),
                    "generic_lm_id": getattr(lipid, "generic_lm_id", None),
                    "refmet": getattr(lipid, "refmet_id", None),
                    "mass": getattr(lipid, "mass", None),
                    "main_class": getattr(lipid, "main_class", None),
                    "sub_class": getattr(lipid, "sub_class", None),
                })
            df_proc = pd.DataFrame(lipid_rows)

            processed_table_container.write(f"Rows: {df_proc.shape[0]}, Columns: {df_proc.shape[1]}")
            # Show the table (read-only) and provide a selection control underneath.
            processed_table_container.dataframe(df_proc, hide_index=True)

            search_container = st.container(border=True)  # reset container for additional content below the table
            # ---- Search field for quick query testing ----
            search_container.subheader("Search / Quick Query")
            search_q = search_container.text_input("Search lipids (name, standardized, LMID, class)", key="processed_search_query")
            if search_q:
                q = search_q.lower()
                mask = (
                    df_proc["input_name"].fillna("").str.lower().str.contains(q)
                    | df_proc["standardized_name"].fillna("").str.lower().str.contains(q)
                    | df_proc["lm_id"].fillna("").astype(str).str.lower().str.contains(q)
                    | df_proc["generic_lm_id"].fillna("").astype(str).str.lower().str.contains(q)
                    | df_proc["main_class"].fillna("").str.lower().str.contains(q)
                    | df_proc["sub_class"].fillna("").str.lower().str.contains(q)
                )
                df_matches = df_proc[mask]
                search_container.write(f"Matches: {df_matches.shape[0]}")
                if df_matches.empty:
                    search_container.info("No matches found for your query.")
                else:
                    search_container.dataframe(df_matches, hide_index=True)

                    # allow quick selection of one of the matches
                    opts = [f"{r['index']}: {r['input_name']}" for _, r in df_matches.iterrows()]
                    sel = search_container.selectbox("Select matched lipid to inspect", opts, key="search_match_select")
                    try:
                        sel_idx = int(sel.split(":", 1)[0])
                    except Exception:
                        sel_idx = None
                    if sel_idx is not None:
                        st.session_state["selected_lipid_idx"] = sel_idx
                        lipid = dataset.lipids[sel_idx]
                        search_container.markdown(f"**Selected:** {lipid.input_name} — LMID: {getattr(lipid, 'lm_id', None)}")
                        # Show simple per-sample bar chart for the selected lipid
                        vals = getattr(lipid, "values", {}) or {}
                        if vals:
                            try:
                                sample_order = dataset.sample_names
                            except Exception:
                                sample_order = list(vals.keys())
                            rows = [{"sample": sid, "value": vals[sid]} for sid in sample_order if sid in vals]
                            df_vals = pd.DataFrame(rows)
                            fig = px.bar(
                                df_vals,
                                x="sample",
                                y="value",
                                title=f"Quantitation for {lipid.input_name}",
                                color="sample",
                                color_discrete_sequence=px.colors.qualitative.Plotly,
                            )
                            search_container.plotly_chart(fig, use_container_width=True, key="search_selected_lipid_chart")

                    # show suggested Query code for reproducing the search
                    try:
                        from lipidmaps.data.models.query import attr_contains
                        q_code = " | ".join([
                            "attr_contains('input_name', '%s')" % search_q,
                            "attr_contains('standardized_name', '%s')" % search_q,
                            "attr_contains('lm_id', '%s')" % search_q,
                            "attr_contains('generic_lm_id', '%s')" % search_q,
                            "attr_contains('main_class', '%s')" % search_q,
                            "attr_contains('sub_class', '%s')" % search_q,
                        ])
                        search_container.code(f"# Equivalent dataset Query\nq = {q_code}\nresults = dataset.query_lipids(q, combine='or')\nlen(results)", language="python")
                    except Exception:
                        pass

            # ---- Mass range query ----
            with search_container.expander("Search by mass range", expanded=False):
                try:
                    min_mass = st.number_input("Min mass (inclusive)", value=0.0, format="%.6f", key="mass_min")
                    max_mass = st.number_input("Max mass (inclusive)", value=1000.0, format="%.6f", key="mass_max")
                except Exception:
                    min_mass = 0.0
                    max_mass = 1000.0

                if st.button("Find lipids by mass", key="mass_search"):
                    try:
                        def mass_pred(l):
                            m = getattr(l, "mass", None)
                            if m is None:
                                return False
                            try:
                                return (m >= float(min_mass)) and (m <= float(max_mass))
                            except Exception:
                                return False

                        results = dataset.query_lipids(mass_pred)
                        st.write(f"Matches: {len(results)}")
                        if not results:
                            st.info("No lipids found in that mass range.")
                        else:
                            rows = [{"input_name": r.input_name, "mass": getattr(r, "mass", None), "lm_id": getattr(r, "lm_id", None)} for r in results]
                            df_mass = pd.DataFrame(rows)
                            st.dataframe(df_mass, hide_index=True)

                            opts = [f"{i}: {r.input_name} (mass={getattr(r,'mass',None)})" for i, r in enumerate(results)]
                            sel = st.selectbox("Select lipid to inspect", opts, key="mass_select")
                            try:
                                sel_idx = int(sel.split(":", 1)[0])
                            except Exception:
                                sel_idx = None
                            if sel_idx is not None:
                                lipid = results[sel_idx]
                                st.session_state["selected_lipid_idx"] = next((i for i, L in enumerate(dataset.lipids) if L.input_name == lipid.input_name), None)
                                st.markdown(f"**Selected:** {lipid.input_name} — mass: {getattr(lipid, 'mass', None)} — LMID: {getattr(lipid, 'lm_id', None)}")
                    except Exception as e:
                        st.error(f"Mass query failed: {e}")

            # Simple per-sample listing using dataset helper
            st.subheader("Per-sample lipid values")
            sample_opts = dataset.list_sample_names() if getattr(dataset, 'samples', None) else []
            if sample_opts:
                sample_sel = st.selectbox("Select sample to list lipid values", sample_opts, key="sample_list_select")
                data = dataset.get_lipid_values_for_samples(sample_sel)
                # display array of objects as bar chart and table
                try:
                    df_vals = pd.DataFrame([d for d in data if d.get("value") is not None])
                except Exception:
                    df_vals = pd.DataFrame(data)

                    if df_vals.empty:
                        st.info(f"No lipid values for sample {sample_sel}.")
                    else:
                        fig = px.bar(
                            df_vals,
                            x="input_name",
                            y="value",
                            title=f"Lipid values for sample {sample_sel}",
                            color="input_name",
                            color_discrete_sequence=px.colors.qualitative.Plotly,
                        )
                        fig.update_layout(xaxis_title="Lipid", yaxis_title="Value", xaxis_tickangle=-45)
                        st.plotly_chart(fig, use_container_width=True, key="sample_lipid_bar_chart")
            else:
                st.info("No samples available in dataset.")

            # Mean value per main class for selected sample
            if sample_opts:
                st.subheader("Mean value per Main Class")
                try:
                    main_classes = df_proc["main_class"].dropna().unique().tolist()
                except Exception:
                    main_classes = []

                class_rows = []
                for mc in main_classes:
                    class_lipids = [l for l in dataset.lipids if getattr(l, "main_class", None) == mc]
                    if not class_lipids:
                        continue
                    mean_val = dataset.mean_value_for_lipids(sample_sel, class_lipids, skip_missing=True)
                    class_rows.append({"main_class": mc, "mean_value": mean_val})
                # After collecting mean values per class, render a single chart
                if class_rows:
                    df_class = pd.DataFrame(class_rows).sort_values(by="mean_value", ascending=False)
                    fig = px.bar(
                        df_class,
                        x="main_class",
                        y="mean_value",
                        title=f"Mean per main class for sample {sample_sel}",
                        color="main_class",
                        color_discrete_sequence=px.colors.qualitative.Set2,
                    )
                    fig.update_layout(xaxis_title="Main class", yaxis_title="Mean value", xaxis_tickangle=-45)
                    st.plotly_chart(fig, use_container_width=True, key="mean_value_per_main_class")
                else:
                    st.info("No main class information available to compute means.")

            # Provide a stable selection UI (selectbox) for choosing a lipid to view sample values bar chart
            options = [f"{i}: {getattr(l, 'input_name', '')}" for i, l in enumerate(dataset.lipids)]
            if options:
                sel = st.selectbox("Select lipid to view sample values bar chart", options, key="lipid_select")
                try:
                    selected_idx = int(sel.split(":", 1)[0])
                except Exception:
                    selected_idx = None
                if selected_idx is not None:
                    st.session_state["selected_lipid_idx"] = selected_idx
                    lipid = dataset.lipids[selected_idx]
                    st.subheader(f"Sample values bar chart for: {lipid.input_name}")

                    # Build bar chart of quantitation values
                    values = getattr(lipid, "values", {}) or {}
                    if values:
                        try:
                            sample_order = dataset.sample_names
                        except Exception:
                            sample_order = list(values.keys())

                        rows = []
                        for sid in sample_order:
                            if sid in values:
                                rows.append({"sample": sid, "value": values[sid]})
                        if not rows:
                            rows = [{"sample": k, "value": v} for k, v in values.items()]

                        df_vals = pd.DataFrame(rows)
                        fig = px.bar(
                            df_vals,
                            x="sample",
                            y="value",
                            title=f"Quantitation for {lipid.input_name}",
                            color="sample",
                            color_discrete_sequence=px.colors.sequential.Viridis,
                        )
                        fig.update_layout(xaxis_title="Sample", yaxis_title="Value")
                        st.plotly_chart(fig, use_container_width=True, key="sample_lipid_bar_chart_2")
                    else:
                        st.info("No quantitation values available for this lipid.")
            else:
                st.info("No lipids available to select.")

            # Note: Query examples are shown in the README; see README.md for usage of dataset.query_lipids

            # ---- Pie Charts ----

            # Main class
            if "main_class" in df_proc:
                counts = df_proc["main_class"].value_counts()
                if len(counts) > 0:
                    st.subheader("Main Class Distribution")
                    fig = px.pie(values=counts.values, names=counts.index)
                    st.plotly_chart(fig, use_container_width=True, key="main_class_distribution")

            # LMID Found
            if "lm_id" in df_proc:
                lm_counts = df_proc["lm_id"].notna().value_counts()
                if len(lm_counts) > 0:
                    st.subheader("LM ID Found Distribution")
                    labels = ["LM ID Found" if x else "Not Found" for x in lm_counts.index]
                    fig = px.pie(values=lm_counts.values,
                                names=labels,
                                color=labels,
                                color_discrete_map={"LM ID Found": "#4caf50", "Not Found": "#ffb3b3"})
                    st.plotly_chart(fig, use_container_width=True, key="lm_id_found_distribution")
            # Generic LMID
            if "generic_lm_id" in df_proc:
                gen_counts = df_proc["generic_lm_id"].notna().value_counts()
                if len(gen_counts) > 0:
                    st.subheader("Generic LM ID Found Distribution")
                    labels = ["Found" if x else "Not Found" for x in gen_counts.index]
                    # Use green for found, red for not found
                    color_map = {"Found": "#4caf50", "Not Found": "#ffb3b3"}
                    fig = px.pie(values=gen_counts.values,
                                names=labels,
                                color=labels,
                                color_discrete_map=color_map)
                    st.plotly_chart(fig, use_container_width=True, key="generic_lm_id_found_distribution")

            # Neither LMID nor Generic LMID
            if "lm_id" in df_proc and "generic_lm_id" in df_proc:
                neither = (~df_proc["lm_id"].notna()) & (~df_proc["generic_lm_id"].notna())
                neither_counts = neither.value_counts()
                if len(neither_counts) > 0:
                    st.subheader("Neither LM ID nor Generic LM ID Found")
                    labels = ["Neither Found" if x else "At Least One Found" for x in neither_counts.index]
                    color_map = {"Neither Found": "#ffb3b3", "At Least One Found": "#4caf50"}
                    fig = px.pie(values=neither_counts.values,
                                names=labels,
                                color=labels,
                                color_discrete_map=color_map)
                    st.plotly_chart(fig, use_container_width=True, key="neither_lm_id_found_distribution")


    # --------------------------------------------------------------
    # GENERIC LMID ASSIGNMENT
    # --------------------------------------------------------------
    if generic_lm_id_button and st.session_state["dataset"] is not None:
        ds = st.session_state["dataset"]
        updated = ds.fill_headgroups_from_name()
        updated = ds.fill_generic_lm_ids_from_headgroups()
        st.session_state["generic_lm_id_assigned"] = True
        st.success(f"Updated {updated} lipids using headgroup mapping.")
        st.rerun()  # refresh processed page


    # --------------------------------------------------------------
    # REACTIONS TAB
    # --------------------------------------------------------------

    if fetch_reactions_button and st.session_state["dataset"] is not None:
        ds = st.session_state["dataset"]

        try:
            tg = st.session_state.get("taxonomy_group", None)   
            if tg and tg != "all":
                reactions = ds.fetch_reactions_by_lm_id(
                    reaction_type="species-level",
                    only_lipid_components=False,
                    taxonomy_group=tg,
                )
            else:
                # 'all' selected — omit taxonomy_group from request
                reactions = ds.fetch_reactions_by_lm_id(
                    reaction_type="species-level",
                    only_lipid_components=False,
                )
            st.session_state["reactions"] = reactions

            st.success(f"Fetched {len(reactions)} reactions.")
            st.session_state["reactions_fetched"] = True
            st.rerun()  # IMPORTANT: stable, never clears processed page
        except Exception as e:
            st.error(f"Error fetching reactions: {e}")

    with tabs[tab_index["reactions"]]:
        st.subheader("Reactions for LM IDs")
        st.write(f"Taxonomy group filter: {st.session_state.get('taxonomy_group', 'all')}")

        if not getattr(st.session_state["dataset"], "reactions", None):
            st.info("No reactions fetched yet. Use Tools → Fetch reactions by LM ID.")
        else:
            # Build reactions table, handling pathway dicts or objects
            rxn_rows = []
            for r in getattr(st.session_state["dataset"], "reactions", []):
                # Create clickable LM ID links for reactants
                reactants_links = []
                for c in getattr(r, "reactants", []):
                    lm_id = getattr(c, "compound_lm_id", None)
                    if lm_id:
                        reactants_links.append(f'<a href="https://lipidmaps.org/databases/lmsd/{lm_id}" target="_blank">{lm_id}</a>')
                reactants_str = ", ".join(reactants_links) if reactants_links else "N/A"
                
                # Create clickable LM ID links for products
                products_links = []
                for c in getattr(r, "products", []):
                    lm_id = getattr(c, "compound_lm_id", None)
                    if lm_id:
                        products_links.append(f'<a href="https://lipidmaps.org/databases/lmsd/{lm_id}" target="_blank">{lm_id}</a>')
                products_str = ", ".join(products_links) if products_links else "N/A"
                
                pathways_str = ", ".join([
                    s for s in ((p.get("name") if isinstance(p, dict) else getattr(p, "pathway_name", None)) for p in getattr(r, "pathways", [])) if s
                ])
                # Create clickable EC number links to BRENDA
                ec_links = []
                for p in getattr(r, "proteins", []):
                    ec = p.get("ec_number") if isinstance(p, dict) else getattr(p, "ec_number", None)
                    if ec:
                        ec_links.append(f'<a href="https://www.brenda-enzymes.org/enzyme.php?ecno={ec}" target="_blank">{ec}</a>')
                ec_str = ", ".join(ec_links) if ec_links else "N/A"
                genes_str = ", ".join([
                    s for s in ((p.get("gene_name") if isinstance(p, dict) else getattr(p, "gene_name", None)) for p in getattr(r, "genes", [])) if s
                ])
                organisms = ", ".join(r.organisms) if r.organisms else "N/A"
                rxn_rows.append({
                    "reaction_id": getattr(r, "reaction_id", None),
                    "reaction_name": getattr(r, "reaction_name", None),
                    "reactants": reactants_str,
                    "products": products_str,
                    "pathways": pathways_str,
                    "ec_number": ec_str,
                    "organisms": organisms,
                })

            rxn_df = pd.DataFrame(rxn_rows)
            st.write(f"Rows: {rxn_df.shape[0]}, Columns: {rxn_df.shape[1]}")
            # Display with HTML rendering for clickable EC number links
            st.write(rxn_df.to_html(escape=False, index=False), unsafe_allow_html=True)
        
            # For each reaction, show a small graph and metadata (reactants -> reaction -> products)
            def reaction_to_dot(reaction):
                def label_for(component):
                    if component is None:
                        return ""
                    if isinstance(component, dict):
                        return component.get("compound_lm_id") or component.get("compound_name")
                    return getattr(component, "compound_lm_id", None) or getattr(component, "compound_name", None)

                reactants = [label_for(c) for c in getattr(reaction, "reactants", []) or []]
                products = [label_for(c) for c in getattr(reaction, "products", []) or []]

                rxn_label = (getattr(reaction, "reaction_name", "") or "").replace('"', '\\"')
                enzymes = []
                if getattr(reaction, "enzyme_ids", None):
                    enzymes = list(getattr(reaction, "enzyme_ids") or [])
                elif getattr(reaction, "proteins", None):
                    enzymes = [ (p.get("ec_number") if isinstance(p, dict) else getattr(p, "ec_number", None)) for p in getattr(reaction, "proteins", []) ]
                enzyme_label = ", ".join([e for e in enzymes if e])

                lines = ["digraph reaction {", "  rankdir=LR;", "  node [shape=box, style=filled, fillcolor=\"#EFEFEF\"];"]

                for i, r in enumerate(reactants):
                    safe = str(r).replace('"', '\\"')
                    lines.append(f'  r{i} [label="{safe}"];')

                for j, p in enumerate(products):
                    safe = str(p).replace('"', '\\"')
                    lines.append(f'  p{j} [label="{safe}"];')

                # reaction centre
                rxn_safe = rxn_label if rxn_label else getattr(reaction, "reaction_id", "reaction")
                rxn_safe = rxn_safe.replace('"', '\\"')
                lines.append(f'  rxn [label="{rxn_safe}\n{enzyme_label}", shape=diamond, style=filled, fillcolor=\"#FFDDAA\"];')

                for i in range(len(reactants)):
                    lines.append(f'  r{i} -> rxn;')
                for j in range(len(products)):
                    lines.append(f'  rxn -> p{j};')

                lines.append('}')
                return "\n".join(lines)

            # Provide a dropdown to select a reaction and view its graph + metadata
            reactions = getattr(st.session_state["dataset"], "reactions", [])
            if reactions:
                reaction_options = [f"{i}: {getattr(r, 'reaction_name', '')} ({getattr(r, 'reaction_id', '')})" for i, r in enumerate(reactions)]
                sel = st.selectbox("Select reaction to view", ["(none)"] + reaction_options, key="reaction_select")
                if sel and sel != "(none)":
                    try:
                        selected_idx = int(sel.split(":", 1)[0])
                    except Exception:
                        selected_idx = None

                    if selected_idx is not None:
                        r = reactions[selected_idx]
                        try:
                            dot = reaction_to_dot(r)
                            st.graphviz_chart(dot)
                        except Exception as e:
                            st.write(f"Could not render graph: {e}")

                        # metadata: enzymes and pathways
                        enzymes = getattr(r, "enzyme_ids", None) or []
                        if not enzymes and getattr(r, "proteins", None):
                            enzymes = [ (p.get("ec_number") if isinstance(p, dict) else getattr(p, "ec_number", None)) for p in getattr(r, "proteins", []) ]
                        st.write("**Enzymes / EC numbers:**", ", ".join([e for e in enzymes if e]))

                        pathways_list = []
                        for p in getattr(r, "pathways", []) or []:
                            if isinstance(p, dict):
                                pathways_list.append(p.get("name") or p.get("pathway_name"))
                            else:
                                pathways_list.append(getattr(p, "pathway_name", None) or getattr(p, "name", None))
                        st.write("**Pathways:**", ", ".join([p for p in pathways_list if p]))

            # Build pathways table, handling pathway dicts or objects and multiple pathways per reaction
            pathway_rows = []
            for r in getattr(st.session_state["dataset"], "reactions", []):
                for p in getattr(r, "pathways", []) or []:
                    if isinstance(p, dict):
                        pid = p.get("id")
                        name = p.get("name")
                        desc = p.get("wikipathways_description") or p.get("description")
                        org = p.get("organism")
                    else:
                        pid = getattr(p, "id", None)
                        name = getattr(p, "name", None) 
                        desc = getattr(p, "wikipathways_description", None)
                        org = getattr(p, "organism", None)

                    pathway_rows.append({
                        "pathway_name": name,
                        "Description": desc,
                        "organism": org,
                    })

            if pathway_rows:
                # Deduplicate pathways by (pathway_id, pathway_name)
                unique_rows = []
                seen = set()
                for row in pathway_rows:
                    key = (row.get("pathway_id"), row.get("pathway_name"))
                    if key in seen:
                        continue
                    seen.add(key)
                    unique_rows.append(row)

                pathway_df = pd.DataFrame(unique_rows)
                st.subheader("Pathways for Reactions")
                st.dataframe(pathway_df, hide_index=True)
            else:
                st.info("No pathway entries available for fetched reactions.")
            # Show lipids annotated with reactions
            dataset = st.session_state.get("dataset")
            if dataset:
                lipids_with_rxn = dataset.get_lipids_with_reactions()
                if lipids_with_rxn:
                    st.subheader("Lipids Annotated with Reactions")

                    lip_rows = []
                    for lip in lipids_with_rxn:
                        reactions_str = ", ".join([
                            s for s in (getattr(r, "reaction_name", None) for r in (lip.reactions or [])) if s
                        ])
                        lip_rows.append({
                            "input_name": getattr(lip, "input_name", None),
                            "lm_id": getattr(lip, "lm_id", None),
                            "generic_lm_id": getattr(lip, "generic_lm_id", None),
                            "reactions": reactions_str,
                        })
                    lip_df = pd.DataFrame(lip_rows)

                    st.dataframe(lip_df, hide_index=True)
                else:
                    st.info("No lipids in the dataset have reactions.")


    # --------------------------------------------------------------
    # VALIDATION TAB
    # --------------------------------------------------------------
    with tabs[tab_index["validation"]]:
        if not st.session_state["show_validation_section"]:
            st.info("Validation report is hidden.")
        elif not st.session_state["has_validation_report"]:
            st.info("Run processing with verification enabled.")
        else:
            st.subheader("Validation Report")
            st.write(f"Passed: {st.session_state['validation_passed']}")

            issues = getattr(st.session_state, "validation_issues", [])
            st.write(f"Issues: {len(issues)}")

            # summary = getattr(st.session_state, "validation_summary", None)
            # st.write(f"Summary: {summary}")

            show_all = st.checkbox("Show all issues", key="show_all_issues")
            to_show = issues if show_all else issues[:5]

            for issue in to_show:
                st.write("-", getattr(issue, "message", str(issue)))

            if len(issues) > 5 and not show_all:
                st.write(f"...and {len(issues) - 5} more.")

if __name__ == "__main__":
    main()