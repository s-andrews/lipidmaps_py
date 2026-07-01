import streamlit as st
import pandas as pd
from math import ceil


from lipidmaps.data.models.reaction import ReactionChecker, ReactionData
from lipidmaps.config import LMSD_REACTIONS_BASE_URL



# ---- Example: Your ReactionData objects ----
# reaction_list = [ReactionData(...), ReactionData(...), ...]

def reaction_to_row(r):
    """Flatten reaction row for the table."""
    return {
        "reaction_id": r.reaction_id,
        "reaction_name": r.reaction_name,
        "num_reactants": len(r.reactants),
        "num_products": len(r.products),
        "num_proteins": len(r.proteins),
        "num_pathways": len(r.pathways),
        "has_lipid_components": r.has_lm_main_components(),
        "type": r.reaction_type or "—",
    }


def display_compound_card(compound):
    """Display a single compound as a card."""
    with st.container(border=True):
        # Use display_name() method for best available name
        st.markdown(f"**{compound.display_name()}**")
        if compound.compound_lm_id:
            st.caption(f"LM ID: `{compound.compound_lm_id}`")
        
        # System name
        if compound.compound_sys_name:
            st.markdown(f"<small>*{compound.compound_sys_name}*</small>", unsafe_allow_html=True)

        # Additional details in expander
        with st.expander("Details", expanded=False):
            details = {}
            for key in ["compound_synonyms", "compound_generic_lm_id", "compound_abbrev",
                        "compound_abbrev_chains", "compound_headgroup", "compound_full_struct", "compound_smiles"]:
                value = getattr(compound, key, None)
                if value:
                    details[key.replace("compound_", "").replace("_", " ").title()] = value
            for key, value in details.items():
                st.text(f"{key}: {value}")


def show_full_reaction(r: "ReactionData"):
    st.subheader(f"Reaction {r.reaction_id}: {r.reaction_name}")
    
    # Show reaction metadata
    metadata_col1, metadata_col2 = st.columns(2)
    with metadata_col1:
        if r.reaction_type:
            st.caption(f"Type: **{r.reaction_type}**")
    with metadata_col2:
        organisms = r.organisms
        if organisms:
            st.caption(f"Organisms: **{', '.join(organisms)}**")

    # Display Reactants and Products side-by-side
    st.markdown("### Reactants & Products")
    reactant_col, product_col = st.columns(2)
    
    with reactant_col:
        st.markdown("#### Reactants")
        if r.reactants:
            for comp in r.reactants:
                display_compound_card(comp)
        else:
            st.info("No reactants")
    
    with product_col:
        st.markdown("#### Products")
        if r.products:
            for comp in r.products:
                display_compound_card(comp)
        else:
            st.info("No products")

    # Additional reaction details
    st.markdown("---")
    
    st.markdown("### Proteins")
    if r.proteins:
        proteins_col1, proteins_col2 = st.columns(2)
        with proteins_col1:
            st.write(f"**Count:** {len(r.proteins)}")
        with proteins_col2:
            organisms = r.organisms
            if organisms:
                st.write(f"**Organisms:** {', '.join(organisms)}")
        with st.expander("Protein Details"):
            st.json(r.proteins)
    else:
        st.info("No proteins")

    st.markdown("### Curations")
    if r.curations:
        st.write(f"**Count:** {len(r.curations)}")
        with st.expander("Curation Details"):
            st.json(r.curations)
    else:
        st.info("No curations")

    st.markdown("### Pathways")
    if r.pathways:
        st.write(f"**Count:** {len(r.pathways)}")
        with st.expander("Pathway Details"):
            st.json(r.pathways)
    else:
        st.info("No pathways")
    
    st.markdown("### Genes")
    if r.genes:
        st.write(f"**Count:** {len(r.genes)}")
        with st.expander("Gene Details"):
            st.json(r.genes)
    else:
        st.info("No genes")

# @st.cache_data(ttl=60 * 60 * 2)
def fetch_all_reactions():
	"""Fetch all reactions and cache result for 2 hours."""
	reaction_checker = ReactionChecker(base_url=LMSD_REACTIONS_BASE_URL)
	return reaction_checker.check_reactions(lm_ids="all", only_lipid_components=False, generic_reactions=False)


# Read selection from main app via st.session_state (direct access option)
selected = st.session_state.get("selected_reaction", None)
reaction_list = []

# Load reactions: prefer a `ReactionData` passed via session state, otherwise fetch all
if selected:
    st.write("Selected reaction from main app:", selected)
    # If the session value is a ReactionData instance, use it directly
    if isinstance(selected, ReactionData):
        reaction_list = [selected]
    else:
        # Attempt to interpret the selected value as an ID and find it in the full set
        all_reactions = fetch_all_reactions()
        try:
            sel_id = int(selected)
        except Exception:
            sel_id = None
        if sel_id is not None:
            reaction_list = [r for r in all_reactions.reactions if r.reaction_id == sel_id]
        else:
            reaction_list = []
else:
    all_reactions = fetch_all_reactions()
    reaction_list = all_reactions.reactions

# Default: sort reactions by reaction_id descending
def rid_key(r: "ReactionData") -> int:
    val = getattr(r, "reaction_id", None)
    if val is None:
        return -1
    try:
        return int(val)
    except Exception:
        return -1

reaction_list = sorted(reaction_list, key=rid_key, reverse=True)

# ---- Quick search using Query helper (searches id, name, and compound names/LMIDs)
search_container = st.container(border=True)
search_q = search_container.text_input("Search reactions (id, name, compound)", key="reactions_search_query")
if search_q:
    try:
        from lipidmaps.data.models.query import from_callable

        q_lower = search_q.lower()

        def match_fn(r: ReactionData) -> bool:
            # reaction id
            try:
                if r.reaction_id is not None and q_lower in str(r.reaction_id):
                    return True
            except Exception:
                pass
            # reaction name
            try:
                if r.reaction_name and q_lower in (r.reaction_name or "").lower():
                    return True
            except Exception:
                pass
            # components (reactants/products)
            try:
                for comp in (r.reactants or []) + (r.products or []):
                    for attr in ("compound_name", "compound_lm_id", "compound_generic_lm_id", "compound_sys_name"):
                        val = getattr(comp, attr, None) or ""
                        if q_lower in str(val).lower():
                            return True
            except Exception:
                pass
            return False

        q = from_callable(match_fn)
        filtered = [r for r in reaction_list if q(r)]
        search_container.write(f"Matches: {len(filtered)}")
        if not filtered:
            search_container.info("No matches found for your query.")
        else:
            # keep the same descending sort for matched results
            reaction_list = sorted(filtered, key=rid_key, reverse=True)
    except Exception:
        # If query helpers aren't available for some reason, skip search
        pass
     
# ---- Build Table ----
rows = [reaction_to_row(r) for r in reaction_list]
df = pd.DataFrame(rows)

# ---- Pagination ----
PAGE_SIZE = 20
total_pages = ceil(len(df) / PAGE_SIZE)

if "page" not in st.session_state:
    st.session_state.page = 1

prev_col, next_col = st.columns(2)
with prev_col:
    if st.button("⬅ Previous") and st.session_state.page > 1:
        st.session_state.page -= 1
with next_col:
    if st.button("Next ➡") and st.session_state.page < total_pages:
        st.session_state.page += 1

start = (st.session_state.page - 1) * PAGE_SIZE
end = start + PAGE_SIZE
page_df = df.iloc[start:end]

st.dataframe(page_df, use_container_width=True, hide_index=True)

# ---- Row selection ----
selected_id = st.selectbox(
    "Select a reaction to view full details",
    page_df["reaction_id"],
)

# ---- Display full reaction ----
selected_reaction = next(r for r in reaction_list if r.reaction_id == selected_id)
show_full_reaction(selected_reaction)
