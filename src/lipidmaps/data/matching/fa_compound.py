"""
FA compound release reaction matcher.

Handles reactions where one acyl chain is released as a free fatty acid,
e.g., PC -> LPC, PE -> LPE, DG -> MG.
"""

from typing import Dict, List, Optional

from ..models.species_reaction import (
    ClassReaction,
    ReactionMatchResult,
    ReactionType,
)
from ..utils.chain_parser import AcylChain, LipidStructure
from .base import MatcherContext, ReactionMatcher


class FACompoundMatcher(ReactionMatcher):
    """Matches pairs where an FA chain is released.
    
    For PC 34:1 -> LPC 16:0:
    1. Product must have fewer carbons/bonds than reactant
    2. Difference (34-16=18, 1-0=1) must match an available FA
    3. FA 18:1 must exist in the dataset
    
    Corresponds to R's `is_valid_reaction()` with acyl_add=FALSE.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.FA_RELEASE
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Find pairs where reactant releases FA to form product."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        result.available_fa = set(context.fa_species.keys())
        
        if not reactants or not products:
            return result
        
        # For each reactant, find products that differ by an available FA
        for reactant in reactants:
            for product in products:
                # Product must be smaller (chain released)
                if product.total_carbons >= reactant.total_carbons:
                    continue
                
                carbon_diff, bond_diff = self._get_chain_difference(reactant, product)
                
                # Check if we have an FA matching this difference
                if carbon_diff <= 0 or bond_diff < 0:
                    continue
                
                fa = context.fa_for_difference(carbon_diff, bond_diff)
                if fa is not None:
                    pair = self._create_pair(
                        reactant, product, class_reaction,
                        required_compound=fa
                    )
                    result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result


class FACompoundMatcherFullStructure(ReactionMatcher):
    """FA release matcher with full structure awareness.
    
    When lipids are reported at full structure level (PC 16:0_18:1),
    we can directly check if the released FA matches a specific chain.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.FA_RELEASE
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match using full chain information when available."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        result.available_fa = set(context.fa_species.keys())
        
        if not reactants or not products:
            return result
        
        for reactant in reactants:
            matches = self._find_fa_release_matches(reactant, products, context)
            for product, fa in matches:
                pair = self._create_pair(
                    reactant, product, class_reaction,
                    required_compound=fa
                )
                result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result
    
    def _find_fa_release_matches(
        self,
        reactant: LipidStructure,
        products: List[LipidStructure],
        context: MatcherContext,
    ) -> List[tuple[LipidStructure, AcylChain]]:
        """Find products that can result from FA release from reactant."""
        matches = []
        
        # If reactant has full structure, we know exactly which chains can be released
        if reactant.chains and len(reactant.chains) > 1:
            for i, released_chain in enumerate(reactant.chains):
                # Check if this FA exists in dataset
                if not context.has_fa(released_chain):
                    continue
                
                # Calculate remaining composition
                remaining_carbons = sum(
                    c.carbons for j, c in enumerate(reactant.chains) if j != i
                )
                remaining_bonds = sum(
                    c.double_bonds for j, c in enumerate(reactant.chains) if j != i
                )
                
                # Find matching products
                for product in products:
                    if (product.total_carbons == remaining_carbons and
                        product.total_double_bonds == remaining_bonds):
                        matches.append((product, released_chain))
        else:
            # Fall back to sum composition matching
            for product in products:
                if product.total_carbons >= reactant.total_carbons:
                    continue
                
                carbon_diff = reactant.total_carbons - product.total_carbons
                bond_diff = reactant.total_double_bonds - product.total_double_bonds
                
                if carbon_diff > 0 and bond_diff >= 0:
                    fa = context.fa_for_difference(carbon_diff, bond_diff)
                    if fa is not None:
                        matches.append((product, fa))
        
        return matches
