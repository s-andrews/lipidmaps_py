import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from ..models.reaction import ReactionData, CompoundComponent
from .chain_parser import ChainParser, LipidStructure
from .lipid_reaction_rules import lipid_reaction_rules

logger = logging.getLogger(__name__)


@dataclass
class _ReactionEvaluationResult:
    possible: bool
    reasons: List[str] = field(default_factory=list)
    pairs_info: List[str] = field(default_factory=list)
    possible_species: Dict[str, List[str]] = field(
        default_factory=lambda: {"reactants": [], "products": []}
    )
    explanation: Optional[str] = None

    def add_reason(self, reason: str) -> None:
        self.reasons.append(reason)

    def add_reasons(self, reasons: List[str]) -> None:
        self.reasons.extend(reasons)

    def add_pair(self, reactant: str, product: str, pair_info: str) -> None:
        self.pairs_info.append(pair_info)
        self.possible_species["reactants"].append(reactant)
        self.possible_species["products"].append(product)

    def finalize_possible_species(self) -> None:
        self.possible_species = {
            "reactants": sorted(set(self.possible_species["reactants"])),
            "products": sorted(set(self.possible_species["products"])),
        }

    def to_dict(self) -> Dict[str, Any]:
        return {
            "possible": self.possible,
            "reasons": self.reasons,
            "pairs_info": self.pairs_info,
            "possible_species": self.possible_species,
            "explanation": self.explanation
            if self.explanation is not None
            else ("rule_match" if self.possible else ",".join(self.reasons)),
        }


class ReactionEvaluator:
    """Evaluate ReactionData objects against rule metadata and optional dataset evidence.

    Minimal, self-contained evaluator that annotates ReactionData with
    `possible` (bool) and `possible_explanation` (str).
    """

    def __init__(self, rules: Optional[Dict[str, Any]] = None):
        self.rules = rules or lipid_reaction_rules

    @staticmethod
    def _base_headgroup(headgroup: Optional[str]) -> Optional[str]:
        if not headgroup:
            return None
        return headgroup.replace(" O-", "").replace(" P-", "").strip()

    @staticmethod
    def _component_name(component: CompoundComponent) -> Optional[str]:
        if component.compound_abbrev:
            return component.compound_abbrev
        if component.compound_full_struct and component.compound_headgroup:
            return f"{component.compound_headgroup} {component.compound_full_struct}"
        if component.compound_name:
            return component.compound_name
        return component.compound_headgroup

    def _component_structure(self, component: CompoundComponent) -> Optional[LipidStructure]:
        component_name = self._component_name(component)
        if not component_name:
            return None
        try:
            return ChainParser.parse(component_name)
        except Exception:
            return None

    @staticmethod
    def _has_chain_info(structure: Optional[LipidStructure]) -> bool:
        if structure is None:
            return False
        return bool(structure.backbone or structure.chains or structure.total_carbons or structure.total_double_bonds)

    @staticmethod
    def _acyl_chain_count(structure: Optional[LipidStructure], headgroup_rule: Dict[str, Any]) -> Optional[int]:
        if structure is None:
            return headgroup_rule.get("acyl_chains")
        if structure.backbone is not None:
            return len(structure.chains)
        if structure.chains:
            return len(structure.chains)
        return headgroup_rule.get("acyl_chains")

    @staticmethod
    def _is_sphingoid_present(structure: Optional[LipidStructure]) -> bool:
        return bool(structure and structure.backbone is not None)

    @staticmethod
    def _is_fa_like(component: CompoundComponent) -> bool:
        headgroup = (component.compound_headgroup or "").strip().lower()
        return headgroup in {"fa", "acyl coa"}

    def _resolve_conversion_rule(
        self,
        reactant_headgroup: Optional[str],
        product_headgroup: Optional[str],
    ) -> Tuple[Optional[Dict[str, Any]], Optional[str], Optional[str]]:
        headgroup_rules = self.rules.get("headgroup_reactions", {})
        if not reactant_headgroup or not product_headgroup:
            return None, None, None

        reactant_rule = headgroup_rules.get(reactant_headgroup)
        if reactant_rule:
            direct_rule = reactant_rule.get("conversion_rules", {}).get(product_headgroup)
            if direct_rule is not None:
                return direct_rule, reactant_headgroup, product_headgroup

        base_reactant = self._base_headgroup(reactant_headgroup)
        base_product = self._base_headgroup(product_headgroup)
        if not base_reactant or not base_product:
            return None, None, None

        reactant_rule = headgroup_rules.get(base_reactant)
        if reactant_rule:
            base_rule = reactant_rule.get("conversion_rules", {}).get(base_product)
            if base_rule is not None:
                return base_rule, base_reactant, base_product

        return None, None, None

    def _build_evaluation(
        self,
        possible: bool,
        reasons: List[str],
        pairs_info: Optional[List[str]] = None,
        possible_species: Optional[Dict[str, List[str]]] = None,
        explanation: Optional[str] = None,
    ) -> _ReactionEvaluationResult:
        result = _ReactionEvaluationResult(
            possible=possible,
            reasons=list(reasons),
            pairs_info=list(pairs_info or []),
            possible_species=possible_species or {"reactants": [], "products": []},
            explanation=explanation,
        )
        return result

    @staticmethod
    def _component_display_name(component: CompoundComponent, structure: Optional[LipidStructure] = None) -> str:
        if structure and structure.input_name:
            return structure.input_name
        return (
            component.compound_abbrev
            or component.compound_full_struct
            or component.compound_name
            or component.compound_headgroup
            or "Unknown"
        )

    @staticmethod
    def _total_signature(structure: Optional[LipidStructure]) -> Tuple[int, int]:
        if structure is None:
            return 0, 0
        return structure.total_carbons, structure.total_double_bonds

    @staticmethod
    def _has_comparable_total_info(
        reactant_structure: Optional[LipidStructure],
        product_structure: Optional[LipidStructure],
    ) -> bool:
        if reactant_structure is None or product_structure is None:
            return False
        if not reactant_structure.chains or not product_structure.chains:
            return False
        return (reactant_structure.backbone is None) == (product_structure.backbone is None)

    def _compound_delta_available(
        self,
        dataset: Optional[Any],
        headgroup: str,
        carbon_diff: int,
        bond_diff: int,
    ) -> bool:
        if dataset is None:
            return carbon_diff >= 0 and bond_diff >= 0
        return dataset.check_structure_exist(
            headgroup=headgroup,
            total_carbons=carbon_diff,
            total_double_bonds=bond_diff,
        )

    def _pair_rule_reasons(
        self,
        reactant_structure: Optional[LipidStructure],
        product_structure: Optional[LipidStructure],
        reactant_rule: Dict[str, Any],
        product_rule: Dict[str, Any],
        conversion_rule: Dict[str, Any],
        dataset: Optional[Any] = None,
    ) -> List[str]:
        reasons: List[str] = []
        reactant_linkage = getattr(reactant_structure, "linkage_type", None) or reactant_rule.get("linkage_type")
        product_linkage = getattr(product_structure, "linkage_type", None) or product_rule.get("linkage_type")
        if (
            conversion_rule.get("require_same_linkage")
            and reactant_linkage
            and product_linkage
            and reactant_linkage != product_linkage
        ):
            reasons.append("linkage_mismatch")

        if (reactant_rule.get("has_sphingoid") or product_rule.get("has_sphingoid")) and not (
            self._is_sphingoid_present(reactant_structure)
            or self._is_sphingoid_present(product_structure)
        ):
            reasons.append("missing_sphingoid")

        expected_chain_change = conversion_rule.get("acyl_chain_change")
        reactant_chain_count = self._acyl_chain_count(reactant_structure, reactant_rule)
        product_chain_count = self._acyl_chain_count(product_structure, product_rule)
        if (
            expected_chain_change is not None
            and reactant_chain_count is not None
            and product_chain_count is not None
            and (product_chain_count - reactant_chain_count) != expected_chain_change
        ):
            reasons.append("acyl_chain_mismatch")

        reactant_total_carbons, reactant_total_double_bonds = self._total_signature(reactant_structure)
        product_total_carbons, product_total_double_bonds = self._total_signature(product_structure)
        delta_carbons = reactant_total_carbons - product_total_carbons
        delta_double_bonds = reactant_total_double_bonds - product_total_double_bonds
        external_compounds = conversion_rule.get("reaction_requirements", {}).get("external_compounds", [])

        if expected_chain_change == 0 and self._has_comparable_total_info(reactant_structure, product_structure) and (
            reactant_total_carbons != product_total_carbons
            or reactant_total_double_bonds != product_total_double_bonds
        ):
            reasons.append("acyl_balance_mismatch")
        elif expected_chain_change == -1 and not external_compounds:
            if not self._compound_delta_available(dataset, "FA", delta_carbons, delta_double_bonds):
                reasons.append("missing_fa")
        elif expected_chain_change == 1 and external_compounds and external_compounds[0] == "acylcoa":
            if not self._compound_delta_available(dataset, "acyl CoA", -delta_carbons, -delta_double_bonds):
                reasons.append("missing_acylcoa")

        return reasons

    def _evaluate_without_dataset(self, reaction: ReactionData) -> Dict[str, Any]:
        headgroup_rules = self.rules.get("headgroup_reactions", {})
        reasons: List[str] = []
        main_reactants = [component for component in reaction.reactants if not self._is_fa_like(component)]
        main_products = [component for component in reaction.products if not self._is_fa_like(component)]

        if len(main_reactants) != 1 or len(main_products) != 1:
            return self._build_evaluation(False, ["unsupported_reaction_shape"], explanation="unsupported_reaction_shape").to_dict()

        main_reactant = main_reactants[0]
        main_product = main_products[0]
        reactant_structure = self._component_structure(main_reactant)
        product_structure = self._component_structure(main_product)
        reactant_headgroup = getattr(reactant_structure, "headgroup", None) or main_reactant.compound_headgroup
        product_headgroup = getattr(product_structure, "headgroup", None) or main_product.compound_headgroup

        conversion_rule, rule_reactant_headgroup, rule_product_headgroup = self._resolve_conversion_rule(
            reactant_headgroup,
            product_headgroup,
        )
        if conversion_rule is None:
            return self._build_evaluation(False, ["unsupported_conversion"], explanation="unsupported_conversion").to_dict()

        reactant_rule = headgroup_rules.get(rule_reactant_headgroup or reactant_headgroup, {})
        product_rule = headgroup_rules.get(rule_product_headgroup or product_headgroup, {})
        reasons.extend(
            self._pair_rule_reasons(
                reactant_structure,
                product_structure,
                reactant_rule,
                product_rule,
                conversion_rule,
            )
        )

        explicit_side_components = (len(reaction.reactants) + len(reaction.products)) > 2
        if explicit_side_components:
            all_structures = [self._component_structure(component) for component in reaction.reactants + reaction.products]
            if any(not self._has_chain_info(structure) for structure in all_structures):
                reasons.append("missing_chain_info")
            else:
                reactant_carbons = sum(structure.total_carbons for structure in all_structures[: len(reaction.reactants)] if structure)
                product_carbons = sum(structure.total_carbons for structure in all_structures[len(reaction.reactants):] if structure)
                if reactant_carbons != product_carbons:
                    reasons.append("carbon_mismatch")

        possible = not reasons
        return self._build_evaluation(possible, reasons).to_dict()

    def annotate_reactions(self, reactions: List[ReactionData], dataset: Optional[Any] = None) -> None:
        for r in reactions:
            try:
                res = self.evaluate_reaction(r, dataset=dataset)
                r.evaluation = res

            except Exception:
                logger.exception("Failed to evaluate reaction %s", getattr(r, "reaction_id", r))
                r.evaluation = self._build_evaluation(False, ["evaluation_error"]).to_dict()
    
    def evaluate_reaction(self, reaction: ReactionData, dataset: Optional[Any] = None) -> Dict[str, Any]:
        """
        For each reactant and product generic_lm_id, check which species in the dataset can participate in the reaction.
        Return a dict with possible species for each reactant/product.
        """
        if dataset is None:
            return self._evaluate_without_dataset(reaction)

        reasons: List[str] = []
        result = _ReactionEvaluationResult(possible=False)

        #STEP 1: Build lookup for dataset species by generic_lm_id
        reaction_gids = set(reaction.list_generic_lm_ids())
        reaction_lmids = set(reaction.list_lm_ids())
        #TODO: acyl CoA species currently not included in generic_lm_id list, so we need to also check lm_ids for those. We should update the generic_lm_id assignment logic to include acyl CoA species as well.
        reaction_gids_except_fa_coa = {gid for gid in reaction_gids if not (gid.startswith("LMFA0705"))}   
        generic_lm_ids_exist = dataset.generic_lm_ids_exist(list(reaction_gids_except_fa_coa)) if dataset else False
        if not generic_lm_ids_exist:
            return self._build_evaluation(
                False,
                ["missing_generic_lm_id"],
                explanation=f"Missing generic_lm_id for reactants/products in dataset.{reaction_gids}{reaction_lmids}",
            ).to_dict()

        #STEP 2: For reactions with 1 reactant lipid and 1 product lipid
        if len(reaction.list_nonfa_noncoa_reactant_lm_ids()) == 1 and len(reaction.list_nonfa_noncoa_product_lm_ids()) == 1:
            main_reactant = next((compound for compound in reaction.reactants if getattr(compound, "compound_type", None) == "lm_main" and compound.compound_lm_id != "LMFA01010000" and compound.compound_lm_id != "LMFA07050000"), None)
            reactant_lipids = dataset.get_lipids_for_component(main_reactant) if dataset and main_reactant else []

            main_product = next((compound for compound in reaction.products if getattr(compound, "compound_type", None) == "lm_main" and compound.compound_lm_id != "LMFA01010000" and compound.compound_lm_id != "LMFA07050000"), None)
            product_lipids = dataset.get_lipids_for_component(main_product) if dataset and main_product else []


            for reactant_lipid in reactant_lipids:
                reactant_structure = getattr(reactant_lipid, "structure", None)
                reactant_headgroup = getattr(reactant_structure, "headgroup", None)
                reactant_input_name = getattr(reactant_structure, "input_name", None)
                reactant_rule = self.rules.get("headgroup_reactions", {}).get(reactant_headgroup, {})

                for product_lipid in product_lipids:
                    product_structure = getattr(product_lipid, "structure", None)
                    product_headgroup = getattr(product_structure, "headgroup", None)
                    if product_headgroup in reactant_rule.get("can_convert_to", []):
                        product_input_name = getattr(product_structure, "input_name", None)
                        product_rule = reactant_rule.get("conversion_rules", {}).get(product_headgroup, None)
                        pair_reasons = self._pair_rule_reasons(
                            reactant_structure,
                            product_structure,
                            reactant_rule,
                            self.rules.get("headgroup_reactions", {}).get(product_headgroup, {}),
                            product_rule or {},
                            dataset=dataset,
                        )
                        if not pair_reasons:
                            reasons.append(f"{reactant_input_name} can convert to {product_input_name} based on headgroup reaction rules.")
                            result.possible = True
                            result.add_pair(
                                reactant_lipid.standardized_name or reactant_lipid.input_name,
                                product_lipid.standardized_name or product_lipid.input_name,
                                f"{reactant_lipid.standardized_name} -> {product_lipid.standardized_name}",
                            )
                        else:
                            reasons.extend(pair_reasons)
        else: # Here for if we have multiple lipids in reactants or products
            pass

        result.reasons = reasons
        if not result.possible:
            result.explanation = "No matching species in dataset for this reaction."
            return result.to_dict()

        result.finalize_possible_species()
        result.explanation = (
            f"Possible reactant species: {result.possible_species['reactants']}; "
            f"Possible product species: {result.possible_species['products']}"
        )
        return result.to_dict()


