from .data_manager import DataManager
from .models.reaction import Reaction


def print_annotated_lipids_with_reactions(manager, n=10):
    """
    Print the first n lipids with their annotated reactions.
    """
    if not manager.dataset or not hasattr(manager.dataset, "lipids"):
        print("No dataset or lipids available.")
        return
    for i, lipid in enumerate(manager.dataset.lipids[:n]):
        print(f"Lipid {i+1}: {getattr(lipid, 'input_name')} -> {getattr(lipid, 'standardized_name')} : LMID: {getattr(lipid, 'lm_id', 'Unknown')}")
        if hasattr(lipid, 'reactions') and lipid.reactions:
            for rxn in lipid.reactions:
                # rxn is a SampleReactionInfo (Pydantic model), not a dict
                reaction_name = getattr(rxn, 'reaction_name', None) or 'N/A'
                print(f"Reaction: {reaction_name} ID: {getattr(rxn, 'reaction_id', 'Unknown')} Role: {getattr(rxn, 'role', 'N/A')}")
        else:
            pass
            # print("  No reactions annotated.")

