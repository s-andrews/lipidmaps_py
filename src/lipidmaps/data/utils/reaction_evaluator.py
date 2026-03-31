import re
from typing import Any, Dict, List, Optional, Set

from ..models.reaction import ReactionData, CompoundComponent
from .lipid_reaction_rules import lipid_reaction_rules
from .headgroups import lm_id_to_headgroup


class ReactionEvaluator:
    """Evaluate ReactionData objects against rule metadata and optional dataset evidence.

    Minimal, self-contained evaluator that annotates ReactionData with
    `possible` (bool) and `possible_explanation` (str).
    """

    def __init__(self, rules: Optional[Dict[str, Any]] = None):
        self.rules = rules or lipid_reaction_rules

    def annotate_reactions(self, reactions: List[ReactionData], dataset: Optional[Any] = None) -> None:
        for r in reactions:
            try:
                res = self.evaluate_reaction(r, dataset=dataset)
                if res['possible']:
                    print(f"    Reaction {r.reaction_id} evaluation: {res['possible']} - {res['explanation']}   ")
                r.evaluation = res

            except Exception as e:
                # ensure reactions are always annotated even if evaluator errors
                # print(f"Error evaluating reaction {r.reaction_id}: {str(e)}")
                r.possible = False
    
    def evaluate_reaction(self, reaction: ReactionData, dataset: Optional[Any] = None) -> Dict[str, Any]:
        """
        For each reactant and product generic_lm_id, check which species in the dataset can participate in the reaction.
        Return a dict with possible species for each reactant/product.
        """
        reasons: List[str] = []
        possible = False
        pairs_info: List[str] = []
        possible_species: Dict[str, List[str]] = {"reactants": [], "products": []}

        #STEP 1: Build lookup for dataset species by generic_lm_id
        reaction_gids = set(reaction.list_generic_lm_ids())
        reaction_lmids = set(reaction.list_lm_ids())
        #TODO: acyl CoA species currently not included in generic_lm_id list, so we need to also check lm_ids for those. We should update the generic_lm_id assignment logic to include acyl CoA species as well.
        reaction_gids_except_fa_coa = {gid for gid in reaction_gids if not (gid.startswith("LMFA0705"))}   
        generic_lm_ids_exist = dataset.generic_lm_ids_exist(list(reaction_gids_except_fa_coa)) if dataset else False
        if not generic_lm_ids_exist:
            return {"possible": False, "explanation": f"Missing generic_lm_id for reactants/products in dataset.{reaction_gids}{reaction_lmids}"}
        else:
            ...
            print(f"Evaluating reaction {reaction.reaction_id} - {reaction.reaction_name} with generic_lm_ids: ")
            print(f"    {lm_id_to_headgroup.get(reaction.list_nonfa_noncoa_reactant_lm_ids()[0], 'Unknown')} -> {lm_id_to_headgroup.get(reaction.list_nonfa_noncoa_product_lm_ids()[0], 'Unknown')}")
            # print(f"    LMIDS:   {reaction_lmids}")
            # print(f"    Generic: {reaction_gids}")

        #STEP 2: For reactions with 1 reactant lipid and 1 product lipid
        if len(reaction.list_nonfa_noncoa_reactant_lm_ids()) == 1 and len(reaction.list_nonfa_noncoa_product_lm_ids()) == 1:
            main_reactant = next((compound for compound in reaction.reactants if getattr(compound, "compound_type", None) == "lm_main" and compound.compound_lm_id != "LMFA01010000" and compound.compound_lm_id != "LMFA07050000"), None)
            reactant_lipids = dataset.get_lipids_for_component(main_reactant) if dataset and main_reactant else []
            # if reactant_lipids:
            #     print(f" Reactant lipids: {', '.join([lipid.standardized_name for lipid in reactant_lipids[:5]])}")

            main_product = next((compound for compound in reaction.products if getattr(compound, "compound_type", None) == "lm_main" and compound.compound_lm_id != "LMFA01010000" and compound.compound_lm_id != "LMFA07050000"), None)
            product_lipids = dataset.get_lipids_for_component(main_product) if dataset and main_product else []
            # if product_lipids:
            #     print(f" Product lipids: {', '.join([lipid.standardized_name for lipid in product_lipids[:5]])}")


            for reactant_lipid in reactant_lipids:
                reactant_headgroup = getattr(reactant_lipid.structure, "headgroup", None)
                reactant_input_name = getattr(reactant_lipid.structure, "input_name", None) 
                reactant_rule = self.rules.get("headgroup_reactions", {}).get(reactant_headgroup, {})
                reactant_chains = getattr(reactant_lipid.structure, "chains", [])
                reactant_total_carbons = getattr(reactant_lipid.structure, "total_carbons", 0)
                reactant_total_double_bonds = getattr(reactant_lipid.structure, "total_double_bonds", 0)
                reactant_linkage_type = getattr(reactant_lipid.structure, "linkage_type", None)

                for product_lipid in product_lipids:
                    product_headgroup = getattr(product_lipid.structure, "headgroup", None)
                    if product_headgroup in reactant_rule.get("can_convert_to", []):
                        product_input_name = getattr(product_lipid.structure, "input_name", None)
                        product_rule = reactant_rule["conversion_rules"].get(product_headgroup, None) 

                        product_chains = getattr(product_lipid.structure, "chains", [])
                        product_total_carbons = getattr(product_lipid.structure, "total_carbons", 0)
                        product_total_double_bonds = getattr(product_lipid.structure, "total_double_bonds", 0)
                        product_linkage_type = getattr(product_lipid.structure, "linkage_type", None)

                        #CONDITION 1 IF BOTH SAME LINKAGE: Check if acyl chain count matches rule requirement (if specified)
                        if reactant_linkage_type == product_linkage_type:
                            if reactant_headgroup == "LPC" and product_headgroup == "PC":
                                print(f"    Same linkage - Evaluating possible conversion: {reactant_input_name} -> {product_input_name} ")
                            if product_rule is None or product_rule.get("acyl_chain_change") == 0:
                                if reactant_total_carbons == product_total_carbons and reactant_total_double_bonds == product_total_double_bonds:
                                    possible = True
                                    reasons.append(f"{reactant_input_name} can convert to {product_input_name} based on headgroup reaction rules.")
                                    pairs_info.append(f"{reactant_lipid.standardized_name} -> {product_lipid.standardized_name}")
                                
                            elif product_rule.get("acyl_chain_change")==-1:
                                if product_rule.get("reaction_requirements").get("external_compounds") == []:
                                    # print(f"For {reactant_input_name}->{product_input_name}, we need FA {reactant_total_carbons-product_total_carbons}:{reactant_total_double_bonds-product_total_double_bonds}")
                                    if dataset.check_structure_exist(headgroup="FA", total_carbons = {reactant_total_carbons-product_total_carbons}, total_double_bonds = {reactant_total_double_bonds-product_total_double_bonds}):
                                        possible = True
                                        reasons.append(f"{reactant_input_name} can convert to {product_input_name} based on headgroup reaction rules.")
                                        pairs_info.append(f"{reactant_lipid.standardized_name} -> {product_lipid.standardized_name}")
                                        
                                        # print(f"    {reactant_input_name} can convert to {product_input_name} based on headgroup reaction rules.")
                                reasons.append(f"{product_input_name} does not have enough acyl chains to be a possible product of {reactant_input_name}.")

                            elif product_rule.get("acyl_chain_change")==1:
                                if product_rule.get("reaction_requirements").get("external_compounds")[0] == "acylcoa":
                                    # print("***********************")
                                    # print(f"For {reactant_input_name}->{product_input_name}, we need CoA {product_total_carbons - reactant_total_carbons}:{reactant_total_double_bonds-product_total_double_bonds}")
                                    if dataset.check_structure_exist(headgroup="acyl CoA", total_carbons = {product_total_carbons - reactant_total_carbons}, total_double_bonds = {reactant_total_double_bonds-product_total_double_bonds}):
                                        possible = True
                                        reasons.append(f"{reactant_input_name} can convert to {product_input_name} based on headgroup reaction rules.")
                                        pairs_info.append(f"{reactant_lipid.standardized_name} -> {product_lipid.standardized_name}")

                                        # print(f"    {reactant_input_name} can convert to {product_input_name} based on headgroup reaction rules.")
                                reasons.append(f"{product_input_name} does not have enough acyl chains to be a possible product of {reactant_input_name}.")

                            elif product_rule.get("required_compound") == "fa" and not any(chain.get("is_fa", False) for chain in product_chains):
                                reasons.append(f"{product_input_name} does not have the required fatty acid chain to be a possible product of {reactant_input_name}.")
                            else:
                                reasons.append(f"{reactant_input_name} and {product_input_name} do not match in acyl chain composition, but could still be possible based on headgroup conversion rules.")
        else: # Here for if we have multiple lipids in reactants or products
            pass

        evaluation: Dict[str, Any] = {"possible": possible, "reasons": reasons, "pairs_info": pairs_info}
        evaluation["possible_species"] = possible_species

        if not evaluation["possible"]:
            evaluation["explanation"] = "No matching species in dataset for this reaction."
        else:
            evaluation["explanation"] = f"Possible reactant species: {possible_species['reactants']}; Possible product species: {possible_species['products']}"
        return evaluation


