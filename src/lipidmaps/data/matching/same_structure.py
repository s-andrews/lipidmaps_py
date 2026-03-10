"""
Same-structure reaction matcher.

Handles reactions where the acyl chains remain unchanged,
e.g., PC -> PA, PE -> PC, PA -> DG.
"""

from typing import List, Set

from ..models.species_reaction import (
    ClassReaction,
    ReactionMatchResult,
    ReactionType,
)
from ..utils.chain_parser import LipidSpecies
from .base import MatcherContext, ReactionMatcher


class SameStructureMatcher(ReactionMatcher):
    """Matches species pairs where acyl composition is conserved.
    
    This is the simplest reaction type - if PC(34:1) exists,
    we look for PA(34:1). No FA/FACoA compound is required.
    
    Corresponds to R's `get_reaction_same_struct()` function.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.SAME_STRUCTURE
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Find pairs where reactant and product have same chain composition."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        
        if not reactants or not products:
            return result
        
        # Index products by total composition for O(1) lookup
        products_by_composition: dict[tuple[int, int], List[LipidSpecies]] = {}
        for prod in products:
            key = (prod.total_carbons, prod.total_double_bonds)
            if key not in products_by_composition:
                products_by_composition[key] = []
            products_by_composition[key].append(prod)
        
        # Match reactants to products with same composition
        for reactant in reactants:
            key = (reactant.total_carbons, reactant.total_double_bonds)
            matching_products = products_by_composition.get(key, [])
            
            for product in matching_products:
                pair = self._create_pair(reactant, product, class_reaction)
                result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result


class DesaturationMatcher(ReactionMatcher):
    """Matches desaturation reactions (dhCer -> Cer, dhSM -> SM).
    
    These are similar to same-structure but involve adding a double bond
    to the sphingoid backbone.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.DESATURATION
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match desaturation pairs.
        
        For dhCer(d18:0/16:0) -> Cer(d18:1/16:0):
        - Same FA chain
        - Product backbone has +1 double bond
        """
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        
        if not reactants or not products:
            return result
        
        # For desaturation, we need to match:
        # - Same total carbons
        # - Product has +1 bond (the desaturation)
        for reactant in reactants:
            for product in products:
                # Same carbons, product has one more double bond
                if (reactant.total_carbons == product.total_carbons and
                    product.total_double_bonds == reactant.total_double_bonds + 1):
                    
                    # Additional check: if we have full structure info,
                    # the difference should be in the backbone, not FA
                    if reactant.backbone and product.backbone:
                        # Check backbone differs by 1 bond
                        if product.backbone.double_bonds != reactant.backbone.double_bonds + 1:
                            continue
                    
                    pair = self._create_pair(reactant, product, class_reaction)
                    result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result
