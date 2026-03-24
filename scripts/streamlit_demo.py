import streamlit as st
import os
import sys
import pandas as pd
import plotly.express as px
from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.models import reaction
from lipidmaps.data.quantitation import QuantitationAnalyzer, NormalizationMethod

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
        "refmet_failed": False,
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

            transpose_file = st.checkbox(
                "Transpose CSV (lipids as columns)",
                value=False,
                disabled=not file_chosen,
            )

            has_labels = st.checkbox(
                "CSV has header labels",
                value= False,
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

        # ---- RUN / DEBUG ----
        with st.expander("Run demo locally", expanded=False):
            st.markdown("Run the Streamlit demo using your virtual environment:")
            st.code("source venv/bin/activate\nstreamlit run scripts/streamlit_demo.py", language="bash")

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
    tab_labels = ["Preview", "Processed", "Reactions", "Validation", "Parser"]
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
            dataset = process_csv(fp, validate_data=validate_data, use_refmet=use_refmet, use_headgroups=use_headgroups, taxonomy_group=taxonomy_group, transpose_file=transpose_file, has_labels=has_labels)
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

            # Track RefMet API status
            st.session_state["refmet_failed"] = getattr(dataset, "refmet_failed", False)

            # Auto-refresh headgroup -> lipid-index mapping for UI convenience
            try:
                hg_map_idx = {}
                for i, lipid in enumerate(dataset.lipids):
                    try:
                        s = lipid.structure
                    except Exception:
                        s = None
                    if s is None:
                        continue
                    hg = getattr(s, "headgroup", None) or "Unknown"
                    hg_map_idx.setdefault(hg, []).append(i)
                st.session_state["hg_map_idx"] = hg_map_idx
            except Exception:
                # Non-fatal; mapping is only a convenience for the UI
                pass

            # If labels were supplied in the CSV, prefer them for sample groups
            try:
                if has_labels:
                    for s in getattr(dataset, "samples", []) or []:
                        lbl = getattr(s, "label", None)
                        if lbl:
                            s.group = lbl
            except Exception:
                # non-fatal UI convenience
                pass

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

            # Show warning if RefMet API failed
            if st.session_state.get("refmet_failed"):
                processed_table_container.warning(
                    "RefMet API request failed. Lipid standardization and classification data may be incomplete. "
                    "This can happen due to network issues or API unavailability. You can try reprocessing later."
                )

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
            # Show experimental conditions (groups) and allow filtering samples by condition
            try:
                conds = sorted({(s.group or "Unknown") for s in (dataset.samples or [])})
            except Exception:
                conds = []

            cond_sel = None
            if conds:
                cond_sel = st.selectbox("Filter by condition (group)", ["(all)"] + conds, key="condition_select")

            # Build sample options, optionally filtered by selected condition
            try:
                if cond_sel and cond_sel != "(all)":
                    sample_opts = [s.sample_name for s in dataset.samples if (s.group or "Unknown") == cond_sel]
                else:
                    sample_opts = dataset.list_sample_names() if getattr(dataset, 'samples', None) else []
            except Exception:
                sample_opts = dataset.list_sample_names() if getattr(dataset, 'samples', None) else []

            # If labels were read from CSV, display them alongside sample names
            try:
                sample_rows = []
                for s in (dataset.samples or []):
                    if sample_opts and s.sample_name not in sample_opts:
                        continue
                    sample_rows.append({"sample_name": s.sample_name, "group": getattr(s, 'group', None), "label": getattr(s, 'label', None)})
                if sample_rows:
                    st.write("Sample metadata")
                    st.dataframe(pd.DataFrame(sample_rows), hide_index=True)
            except Exception:
                pass

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

            # ---------------------- NORMALIZATION ----------------------
            st.subheader("Normalization")
            try:
                from lipidmaps.data.quantitation import NormalizationMethod, QuantitationAnalyzer
                norm_methods = [m for m in NormalizationMethod]
                method_labels = {m: m.value for m in norm_methods}

                sel_method = st.selectbox("Normalization method", options=norm_methods, format_func=lambda m: m.value, key="norm_method_select")

                # Explanatory text for each normalization method
                method_descriptions = {
                    NormalizationMethod.NONE: (
                        "No normalization. Values are left unchanged so you can inspect raw measurements "
                        "and compare downstream effects of applying other methods."
                    ),
                    NormalizationMethod.TOTAL_LIPID: (
                        "Total-lipid normalization divides each lipid value by the sum of all lipid values "
                        "in the same sample, then optionally multiplies by a scale factor (e.g. 1e6 for ppm). "
                        "This compensates for differences in total loading or sample amount between samples."
                    ),
                    NormalizationMethod.INTERNAL_STANDARD: (
                        "Internal-standard normalization divides each lipid's value by the value of a chosen "
                        "internal standard lipid measured in the same sample. This corrects for instrument "
                        "variation and sample-specific recovery, but requires a consistent internal standard "
                        "present across samples."
                    ),
                    NormalizationMethod.LOG2: (
                        "Log2 transform: computes log2(value + offset) to reduce skew and stabilize variance. "
                        "An offset (default 1) prevents taking the log of zero. Useful when data are heavy-tailed."
                    ),
                    NormalizationMethod.LOG10: (
                        "Log10 transform: computes log10(value + offset) to reduce skew and stabilize variance. "
                        "An offset (default 1) prevents taking the log of zero. Use based on interpretation preference."
                    ),
                    NormalizationMethod.MEDIAN_CENTER: (
                        "Median-centering subtracts the median value (across all lipids) for each sample from each lipid, "
                        "centering sample distributions and removing sample-specific offsets. Useful for comparability."
                    ),
                    NormalizationMethod.ZSCORE: (
                        "Per-lipid z-score: for each lipid, compute (value - mean) / std across samples. "
                        "This centers each lipid to mean 0 with unit variance, highlighting relative changes across samples."
                    ),
                    NormalizationMethod.QUANTILE: (
                        "Quantile normalization forces the distribution of values across samples to be identical by averaging "
                        "ranked values. It is useful when sample distributions should be comparable, but assumes similar global composition."
                    ),
                }

                try:
                    desc = method_descriptions.get(sel_method)
                    if desc:
                        st.info(desc)
                except Exception:
                    pass

    # Parser tab moved to avoid interfering with surrounding try/except blocks.
    # The interactive parser UI will be rendered just before the Validation tab.

                internal_std = None
                if sel_method == NormalizationMethod.INTERNAL_STANDARD:
                    # allow selecting from lipid names or entering free text
                    lipid_names = [getattr(l, 'input_name', '') for l in dataset.lipids]
                    internal_std = st.selectbox("Internal standard (choose lipid)", options=["(none)"] + lipid_names, key="internal_std_select")
                    if internal_std == "(none)":
                        internal_std = st.text_input("Or type internal standard name", value="", key="internal_std_text") or None

                # normalized values are stored separately on each lipid; no in-place overwrite
                scale_factor = None
                if sel_method == NormalizationMethod.TOTAL_LIPID:
                    scale_factor = st.number_input("Scale factor (e.g. 1e6 for ppm)", value=1e6, format="%.0f", key="norm_scale_factor")

                if st.button("Apply normalization", key="apply_normalization"):
                    try:
                        qa = QuantitationAnalyzer(dataset=dataset)
                        # build a reproducible method key encoding method + params
                        method_key = sel_method.value
                        if sel_method == NormalizationMethod.TOTAL_LIPID and scale_factor is not None:
                            method_key = f"{method_key}:scale={float(scale_factor)}"
                        if sel_method == NormalizationMethod.INTERNAL_STANDARD:
                            std_name = internal_std or st.session_state.get("internal_std_text")
                            method_key = f"{method_key}:std={std_name}"
                        if sel_method == NormalizationMethod.NONE:
                            norm_res = {l.input_name: l.values.copy() for l in dataset.lipids}
                        elif sel_method == NormalizationMethod.TOTAL_LIPID:
                            norm_res = qa.normalize_total_lipid(scale_factor=scale_factor)
                        elif sel_method == NormalizationMethod.INTERNAL_STANDARD:
                            std_name = internal_std or st.session_state.get("internal_std_text")
                            if not std_name:
                                st.error("Please provide an internal standard name or select one.")
                                norm_res = None
                            else:
                                norm_res = qa.normalize_internal_standard(std_name)
                        elif sel_method == NormalizationMethod.LOG2:
                            norm_res = qa.normalize_log(base=2)
                        elif sel_method == NormalizationMethod.LOG10:
                            norm_res = qa.normalize_log(base=10)
                        elif sel_method == NormalizationMethod.MEDIAN_CENTER:
                            norm_res = qa.normalize_median_center()
                        elif sel_method == NormalizationMethod.ZSCORE:
                            # zscore per-lipid helper on lipid objects may exist
                            norm_res = {l.input_name: l.zscore() for l in dataset.lipids}
                        else:
                            st.error(f"Normalization {sel_method} not implemented in demo")
                            norm_res = None

                        if norm_res is not None:
                            # convert to DataFrame: rows = lipids, columns = samples
                            df_norm = pd.DataFrame.from_dict(norm_res, orient='index')
                            # order columns by dataset sample names if available
                            try:
                                cols = dataset.list_sample_names() if hasattr(dataset, 'list_sample_names') else getattr(dataset, 'sample_names', None)
                                if cols:
                                    df_norm = df_norm[cols]
                            except Exception:
                                pass

                            # store normalized results on each lipid under the method_key
                            per_lipid_failed = False
                            for l in dataset.lipids:
                                if l.input_name in norm_res:
                                    try:
                                        l.set_normalized(method_key, norm_res[l.input_name])
                                    except Exception:
                                        per_lipid_failed = True

                            # If per-lipid storage failed (some models may forbid new fields),
                            # fall back to storing the full mapping on the dataset object.
                            if per_lipid_failed:
                                try:
                                    store = getattr(dataset, '_normalized_store', None)
                                    if store is None:
                                        setattr(dataset, '_normalized_store', {})
                                        store = dataset._normalized_store
                                    store[method_key] = norm_res
                                except Exception:
                                    # last-resort: attach attribute directly
                                    try:
                                        dataset._normalized_store = {method_key: norm_res}
                                    except Exception:
                                        pass

                            st.subheader("Normalized values")
                            st.write(f"Rows: {df_norm.shape[0]}, Columns: {df_norm.shape[1]}")
                            st.dataframe(df_norm, hide_index=False)
                            # offer CSV download
                            try:
                                csv_bytes = df_norm.to_csv(index=True).encode('utf-8')
                                st.download_button("Download normalized CSV", data=csv_bytes, file_name="normalized.csv", mime="text/csv")
                            except Exception:
                                pass
                    except Exception as e:
                        st.error(f"Normalization failed: {e}")
            except Exception:
                st.info("Normalization tools unavailable for this dataset.")

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
        updated = ds.fill_headgroups_from_names()
        compound_headgroups_updated = ds.fill_compound_headgroups_from_lipids
        generic_lm_ids_updated = ds.fill_generic_lm_ids_from_headgroups()
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
            # Show plausible reactions filtered by headgroup rules
            try:
                plausible = []
                try:
                    plausible = dataset.possible_reactions(getattr(st.session_state, "reactions", []))
                except Exception:
                    from lipidmaps.data.utils.lipid_reaction_rules import reactions_possible_in_dataset
                    plausible = reactions_possible_in_dataset(dataset, getattr(st.session_state, "reactions", []))

                if plausible:
                    st.subheader("Plausible reactions (headgroup-based filter)")
                    pr_rows = []
                    for r in plausible:
                        pr_rows.append({
                            "reaction_id": getattr(r, "reaction_id", None),
                            "reaction_name": getattr(r, "reaction_name", None),
                        })
                    try:
                        st.dataframe(pd.DataFrame(pr_rows), hide_index=True)
                    except Exception:
                        st.write(pr_rows)
            except Exception:
                pass

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
                    # "pathways": pathways_str,
                    "ec_number": ec_str,
                    "organisms": organisms,
                    "possible": getattr(r, "possible", None),
                    "possible_explanation": getattr(r, "possible_explanation", None),
                })

            rxn_df = pd.DataFrame(rxn_rows)[0:20]
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
                        # Show rule evaluation result with detailed diagnostics
                        possible = getattr(r, "possible", None)
                        explanation = getattr(r, "possible_explanation", None)
                        evaluation = getattr(r, "evaluation", None) or {}
                        if possible is not None:
                            color = "green" if possible else "red"
                            st.markdown(f"**Possible:** <span style='color:{color}; font-weight:bold'>{possible}</span>", unsafe_allow_html=True)

                        # human-readable explanation (summary)
                        if explanation:
                            with st.expander("Explanation", expanded=False):
                                st.write(explanation)

                        # structured evaluation details
                        if evaluation:
                            with st.expander("Evaluation details", expanded=False):
                                # show the explanation from the evaluator if present
                                if evaluation.get("explanation"):
                                    st.markdown("**Summary:**")
                                    st.write(evaluation.get("explanation"))

                                details = evaluation.get("details") or []
                                if details:
                                    rows = []
                                    for pdict in details:
                                        rinfo = pdict.get("reactant", {})
                                        pinfo = pdict.get("product", {})
                                        rule = pdict.get("rule", {})
                                        rows.append({
                                            "reactant_hg": rinfo.get("headgroup"),
                                            "reactant_linkage": rinfo.get("linkage"),
                                            "reactant_chain_count": rinfo.get("chain_count"),
                                            "reactant_total_carbons": rinfo.get("total_carbons"),
                                            "reactant_total_DB": rinfo.get("total_double_bonds"),
                                            "product_hg": pinfo.get("headgroup"),
                                            "product_linkage": pinfo.get("linkage"),
                                            "product_chain_count": pinfo.get("chain_count"),
                                            "product_total_carbons": pinfo.get("total_carbons"),
                                            "product_total_DB": pinfo.get("total_double_bonds"),
                                            "require_same_linkage": rule.get("require_same_linkage"),
                                            "required_acyl_chains": rule.get("required_acyl_chains"),
                                            "can_convert_to": ",".join(rule.get("can_convert_to") or []) if rule.get("can_convert_to") else None,
                                        })

                                    try:
                                        df_eval = pd.DataFrame(rows)
                                        st.dataframe(df_eval, hide_index=True)
                                    except Exception:
                                        st.write(rows)

                                # raw JSON view
                                try:
                                    st.subheader("Raw evaluation JSON")
                                    st.json(evaluation)
                                except Exception:
                                    pass

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

            # Reaction-level quantitation: z-scores between two groups
            try:
                if dataset and getattr(dataset, 'reactions', None):
                    st.subheader("Reaction-level comparison (z-scores)")
                    from lipidmaps.data.quantitation import QuantitationAnalyzer
                    analyzer = QuantitationAnalyzer(dataset=dataset)
                    groups = sorted({analyzer._sample_group_name(s) for s in dataset.samples})
                    if len(groups) >= 2:
                        g1 = st.selectbox("Group 1 (numerator)", groups, key="rx_g1")
                        g2 = st.selectbox("Group 2 (denominator)", [g for g in groups if g != g1], key="rx_g2")
                        method = st.selectbox("Flux method", ["ratio", "difference"], index=0, key="rx_method")
                        if st.button("Compute reaction z-scores", key="compute_rx_z"):
                            try:
                                rz = analyzer.reaction_zscores(g1, g2, method=method)
                                if not rz:
                                    st.info("No reaction flux data available for selected groups.")
                                else:
                                    rows = []
                                    for rid, info in rz.items():
                                        rows.append({
                                            "reaction_id": rid,
                                            "reaction_name": info.get("reaction_name"),
                                            "group1_mean": info.get("group1_mean"),
                                            "group2_mean": info.get("group2_mean"),
                                            "zscore": info.get("zscore"),
                                            "n1": info.get("n1"),
                                            "n2": info.get("n2"),
                                        })
                                    df_rz = pd.DataFrame(rows)
                                    try:
                                        df_rz = df_rz.sort_values(by='zscore', key=lambda s: s.abs(), ascending=False)
                                    except Exception:
                                        pass
                                    st.dataframe(df_rz, hide_index=True)
                                    try:
                                        csv_bytes = df_rz.to_csv(index=False).encode('utf-8')
                                        st.download_button(
                                            "Download z-scores CSV",
                                            data=csv_bytes,
                                            file_name=f"reaction_zscores_{g1}_vs_{g2}_{method}.csv",
                                            mime="text/csv",
                                        )
                                    except Exception:
                                        pass
                            except Exception as e:
                                st.error(f"Failed to compute reaction z-scores: {e}")
                    else:
                        st.info("Need at least two groups in dataset to compute reaction comparisons.")
            except Exception:
                pass


    # --------------------------------------------------------------
    # VALIDATION TAB
    # --------------------------------------------------------------
    # Insert Parser tab just before Validation to avoid nesting issues
    with tabs[tab_index["parser"]]:
        st.subheader("Lipid Parser Tester")
        st.write("Enter a lipid species name to parse using the ChainParser.")
        parser_input = st.text_input("Lipid name", value="PC 16:0_18:1", key="parser_input")
        parse_now = st.button("Parse", key="parser_parse")
        if parse_now:
            try:
                from lipidmaps.data.utils.chain_parser import parse_lipid, ChainParser

                parsed = parse_lipid(parser_input)
                # Prefer pydantic v2 model_dump, fall back to dict()
                try:
                    parsed_dict = parsed.model_dump()
                except Exception:
                    try:
                        parsed_dict = parsed.dict()
                    except Exception:
                        parsed_dict = str(parsed)

                st.json(parsed_dict)

                # show available headgroups for reference
                try:
                    hgs = list(ChainParser.HEADGROUPS)
                    st.expander("Available headgroups (from headgroups.py)", expanded=False).write(hgs)
                except Exception:
                    pass
            except Exception as e:
                st.error(f"Parser failed: {e}")

        # If a dataset is loaded, offer the option to summarize parsed structures
        if st.session_state.get("dataset") is not None:
            ds = st.session_state.get("dataset")
            # Build and store a lightweight mapping: headgroup -> list of lipid indices
            if st.button("Refresh structures by headgroup", key="get_struct_by_hg"):
                try:
                    hg_map_idx = {}
                    for i, lipid in enumerate(ds.lipids):
                        try:
                            s = lipid.structure
                        except Exception:
                            s = None
                        if s is None:
                            continue
                        hg = getattr(s, "headgroup", None) or "Unknown"
                        hg_map_idx.setdefault(hg, []).append(i)

                    st.session_state["hg_map_idx"] = hg_map_idx
                except Exception as e:
                    st.error(f"Failed to build headgroup mapping: {e}")

            # Allow clearing the mapping explicitly
            if st.button("Clear headgroup mapping", key="clear_hg_map"):
                if "hg_map_idx" in st.session_state:
                    st.session_state.pop("hg_map_idx", None)
                    st.success("Headgroup mapping cleared.")
                else:
                    st.info("No headgroup mapping to clear.")

            # If a mapping exists in session state, render it persistently so selections survive reruns
            hg_map_idx = st.session_state.get("hg_map_idx")
            if hg_map_idx:
                try:
                    summary = {hg: len(v) for hg, v in hg_map_idx.items()}
                    df_summary = pd.DataFrame(list(summary.items()), columns=["headgroup", "count"]).sort_values("count", ascending=False)
                    st.subheader("Structures by headgroup")
                    st.dataframe(df_summary, hide_index=True)

                    sel = st.selectbox("Select headgroup to list structures", ["(all)"] + sorted(list(hg_map_idx.keys())), key="hg_map_select")
                    if sel and sel != "(all)":
                        indices = hg_map_idx.get(sel, [])
                        liprows = []
                        for idx in indices:
                            try:
                                lip = ds.lipids[idx]
                            except Exception:
                                continue
                            reactions = []
                            for r in getattr(lip, "reactions", []) or []:
                                reactions.append(getattr(r, "reaction_name", None) or getattr(r, "reaction_id", None))
                            liprows.append({
                                "index": idx,
                                "input_name": getattr(lip, "input_name", None),
                                "lm_id": getattr(lip, "lm_id", None),
                                "generic_lm_id": getattr(lip, "generic_lm_id", None),
                                "reactions_count": len(reactions),
                                "reactions": ", ".join([r for r in reactions if r])
                            })
                        if liprows:
                            st.subheader(f"Lipids in headgroup: {sel}")
                            st.dataframe(pd.DataFrame(liprows), hide_index=True)
                        else:
                            st.info(f"No lipids or reactions for headgroup {sel}.")
                    else:
                        # show brief sample of structures for all headgroups
                        sample_rows = []
                        for hg, indices in list(hg_map_idx.items())[:50]:
                            examples = []
                            for idx in indices[:3]:
                                try:
                                    s = ds.lipids[idx].standardized_name or ds.lipids[idx].input_name
                                except Exception:
                                    s = None
                                if s:
                                    examples.append(s)
                            sample_rows.append({"headgroup": hg, "examples": examples})
                        st.dataframe(pd.DataFrame(sample_rows), hide_index=True)
                except Exception as e:
                    st.error(f"Failed to display headgroup mapping: {e}")

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