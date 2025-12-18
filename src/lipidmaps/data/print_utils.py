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
        print(f"Lipid {i+1}: {getattr(lipid, 'input_name', getattr(lipid, 'lm_id', 'Unknown'))}")
        if hasattr(lipid, 'reactions') and lipid.reactions:
            for rxn in lipid.reactions:
                print(f"  Reaction: {rxn.get('reaction_id', 'N/A')} - {rxn.get('reaction_name', 'N/A')}")
        else:
            print("  No reactions annotated.")
        print()
