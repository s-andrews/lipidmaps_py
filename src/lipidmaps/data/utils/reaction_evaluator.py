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
            except Exception as e:
                # ensure reactions are always annotated even if evaluator errors
                setattr(r, "possible", False)
                setattr(r, "possible_explanation", f"error_evaluating:{e}")
                setattr(r, "evaluation", {"possible": False, "explanation": f"error_evaluating:{e}"})
            else:
                setattr(r, "possible", res.get("possible", False))
                setattr(r, "possible_explanation", res.get("explanation", ""))
                # also attach the full evaluation dict for UI selection / detailed view
                try:
                    setattr(r, "evaluation", res)
                except Exception:
                    # best-effort: evaluation may be non-serializable; fall back to basic fields
                    setattr(r, "evaluation", {"possible": res.get("possible", False), "explanation": res.get("explanation", "")})

    def evaluate_reaction(self, reaction: ReactionData, dataset: Optional[Any] = None) -> Dict[str, Any]:
        """
        For each reactant and product generic_lm_id, check which species in the dataset can participate in the reaction.
        Return a dict with possible species for each reactant/product.
        """
        reasons: List[str] = []
        possible = True
        pairs_info: List[Dict[str, Any]] = []
        possible_species: Dict[str, List[str]] = {"reactants": [], "products": []}

        #STEP 1: Build lookup for dataset species by generic_lm_id
        
        reaction_gids = set(reaction.list_generic_lm_ids())
        generic_lm_ids_exist = dataset.generic_lm_ids_exist(list(reaction_gids)) if dataset else False
        if not generic_lm_ids_exist:
            return {"possible": False, "explanation": "One or more generic_lm_id in reaction not found in dataset."}

        #STEP 2: For reactions with 1 reactant lipid and 1 product lipid
        if len(reaction.list_reactant_lm_ids()) == 1 and len(reaction.list_product_lm_ids()) == 1:
            print("Reaction has 1 reactant and 1 product with generic_lm_id.")
            print(f"Reactant : {reaction.reaction_id} - {reaction.reaction_name}")
            # Find the first reactant with compound_type == 'lm_main'
            main_reactant = next((compound for compound in reaction.reactants if getattr(compound, "compound_type", None) == "lm_main"), None)
            component_lipids = dataset.get_lipids_for_component(main_reactant) if dataset and main_reactant else []
            if component_lipids:
                print(f"Dataset has {len(component_lipids)} lipids for component {getattr(main_reactant, 'compound_name', 'None')}.")
                print([lipid.standardized_name for lipid in component_lipids[:5]])  # Print first 5 lipids for inspection
            else:
                print(f"No lipids found in dataset for component {getattr(main_reactant, 'compound_name', 'None')}.")

            main_product = next((compound for compound in reaction.products if getattr(compound, "compound_type", None) == "lm_main"), None)
            product_lipids = dataset.get_lipids_for_component(main_product) if dataset and main_product else []
            if product_lipids:
                print(f"Dataset has {len(product_lipids)} lipids for product {getattr(main_product, 'compound_name', 'None')}.")
                print([lipid.standardized_name for lipid in product_lipids[:5]])  # Print first 5 lipids for inspection
            else:
                print(f"No lipids found in dataset for product {getattr(main_product, 'compound_name', 'None')}.")
        #STEP 3: For each reactant and product, determine possible species based on headgroup reaction rules and dataset species under the same generic_lm_id

        # Build lookup for dataset species by generic_lm_id
        dataset_lipids = getattr(dataset, "lipids", []) if dataset else []
        generic_id_to_species = {}
        for lipid in dataset_lipids:
            generic_id = getattr(lipid, "generic_lm_id", None)
            if generic_id:
                generic_id_to_species.setdefault(generic_id, []).append(lipid)

        def _get_val(compound, field: str):
            if compound is None:
                return None
            if isinstance(compound, dict):
                return compound.get(field)
            return getattr(compound, field, None)

        def _get_headgroup_from_compound(compound) -> Optional[str]:
            if not compound:
                return None
            generic_lm_id = _get_val(compound, "compound_generic_lm_id")
            if generic_lm_id:
                hg = lm_id_to_headgroup.get(generic_lm_id)
                if hg:
                    return hg
            lm_id = _get_val(compound, "compound_lm_id")
            if lm_id:
                hg = lm_id_to_headgroup.get(lm_id)
                if hg:
                    return hg
            hg = _get_val(compound, "compound_headgroup")
            if hg:
                return hg
            return None

        def _get_linkage_from_compound(compound) -> Optional[str]:
            if not compound:
                return None
            hg = _get_headgroup_from_compound(compound) or _get_val(compound, "compound_headgroup") or _get_val(compound, "compound_abbrev") or  _get_val(compound, "compound_name") or ""
            if hg:
                if "P-" in hg:
                    return "ether_vinyl"
                if "O-" in hg:
                    return "ether_alkyl"
                hg_lower = hg.lower()
                if hg_lower.startswith("cer") or hg_lower.startswith("sm") or hg_lower.startswith("dhcer"):
                    return "amide"
            return "ester"

        def _chain_count_from_compound(compound) -> Optional[int]:
            if not compound:
                return None
            # Prefer explicit chain count if present
            cab = _get_val(compound, "compound_abbrev_chains")
            if cab:
                try:
                    return int(cab)
                except Exception:
                    pass
            # Use chain parser for robust chain counting
            from .chain_parser import parse_lipid
            text = _get_val(compound, "compound_abbrev") or _get_val(compound, "compound_name") or _get_val(compound, "compound_full_struct") or ""
            if text:
                parsed = parse_lipid(text)
                if parsed and hasattr(parsed, "chains"):
                    return len(parsed.chains)
            return None

        # For each reactant, check which dataset species under its generic_lm_id can participate
        for react in reaction.reactants:
            gid = _get_val(react, "compound_generic_lm_id")
            if not gid or gid not in generic_id_to_species:
                continue
            hg = _get_headgroup_from_compound(react)
            rules = self.rules.get("headgroup_reactions", {}).get(hg, {})
            conv = rules.get("conversion_rules", {}) or {}
            # For each product, check conversion
            for prod in reaction.products:
                hp = _get_headgroup_from_compound(prod)
                rule = conv.get(hp)
                if not rule:
                    continue
                # For each species under this generic_lm_id, check if it matches rule requirements
                for lip in generic_id_to_species[gid]:
                    # Build a fake comp dict from lipid species for checking
                    comp = lip.__dict__ if hasattr(lip, "__dict__") else lip
                    # Check acyl chain count if required
                    req_ch = rule.get("required_acyl_chains")
                    if req_ch is not None:
                        rc = _chain_count_from_compound(comp)
                        if rc is None or rc != req_ch:
                            continue
                    # Check linkage if required
                    if rule.get("require_same_linkage", True):
                        lr = _get_linkage_from_compound(comp) or rules.get("linkage_type")
                        lp = _get_linkage_from_compound(prod) or self.rules.get("headgroup_reactions", {}).get(hp, {}).get("linkage_type")
                        if lr and lp and lr != lp:
                            continue
                    # If all checks pass, add species lm_id
                    possible_species["reactants"].append(getattr(lip, "lm_id", str(lip)))

        # Repeat for products
        for prod in reaction.products:
            gid = _get_val(prod, "compound_generic_lm_id")
            if not gid or gid not in generic_id_to_species:
                continue
            hp = _get_headgroup_from_compound(prod)
            rules = self.rules.get("headgroup_reactions", {}).get(hp, {})
            conv = rules.get("conversion_rules", {}) or {}
            for react in reaction.reactants:
                hr = _get_headgroup_from_compound(react)
                rule = conv.get(hr)
                if not rule:
                    continue
                for lip in generic_id_to_species[gid]:
                    comp = lip.__dict__ if hasattr(lip, "__dict__") else lip
                    req_ch = rule.get("required_acyl_chains")
                    if req_ch is not None:
                        pc = _chain_count_from_compound(comp)
                        if pc is None or pc != req_ch:
                            continue
                    if rule.get("require_same_linkage", True):
                        lp = _get_linkage_from_compound(comp) or rules.get("linkage_type")
                        lr = _get_linkage_from_compound(react) or self.rules.get("headgroup_reactions", {}).get(hr, {}).get("linkage_type")
                        if lr and lp and lr != lp:
                            continue
                    possible_species["products"].append(getattr(lip, "lm_id", str(lip)))

        evaluation: Dict[str, Any] = {"possible": bool(possible_species["reactants"] or possible_species["products"])}
        evaluation["possible_species"] = possible_species
        if not evaluation["possible"]:
            evaluation["explanation"] = "No matching species in dataset for this reaction."
        else:
            evaluation["explanation"] = f"Possible reactant species: {possible_species['reactants']}; Possible product species: {possible_species['products']}"
        return evaluation

    def all_generic_ids_present(self, reaction: ReactionData, dataset: Optional[Any]) -> bool:
        """
        Return True if all generic_lm_id's in reactants and products are present in the dataset.
        """
        if not dataset or not hasattr(dataset, "lipids"):
            return False
        dataset_gids = {getattr(l, "generic_lm_id", None) for l in dataset.lipids if getattr(l, "generic_lm_id", None)}
        reaction_gids = set()
        for comp in (getattr(reaction, "reactants", []) + getattr(reaction, "products", [])):
            gid = getattr(comp, "compound_generic_lm_id", None)
            if gid:
                reaction_gids.add(gid)
        return reaction_gids.issubset(dataset_gids)
