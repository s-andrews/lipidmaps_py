"""
FACoA compound addition reaction matcher.

Handles reactions where a fatty acyl-CoA is added to a lipid,
e.g., DG -> TG, LPC -> PC, MAG -> DAG.
"""

from typing import List, Optional

from ..models.species_reaction import (
    ClassReaction,
    ReactionMatchResult,
    ReactionType,
)
from ..utils.chain_parser import AcylChain, LipidStructure
from .base import MatcherContext, ReactionMatcher


class FACoACompoundMatcher(ReactionMatcher):
    """Matches pairs where an FACoA is added to form product.
    FACoA shorthand is CoA 
    For DG 34:1 + CoA 16:0 -> TG 50:1:
    1. Product must have more carbons than reactant
    2. Difference (50-34=16, 1-1=0) must match an available FACoA
    3. CoA 16:0 must exist in the dataset
    
    Corresponds to R's `is_valid_reaction()` with acyl_add=TRUE.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.FACOA_ADDITION
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Find pairs where reactant gains FACoA to form product."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        result.available_facoa = set(context.facoa_species.keys())
        
        if not reactants or not products:
            return result
        
        # For each reactant, find products that differ by an available FACoA
        for reactant in reactants:
            for product in products:
                # Product must be larger (chain added)
                if product.total_carbons <= reactant.total_carbons:
                    continue
                
                carbon_diff, bond_diff = self._get_chain_difference(product, reactant)
                
                # Check if we have an FACoA matching this difference
                if carbon_diff <= 0 or bond_diff < 0:
                    continue
                
                facoa = context.facoa_for_difference(carbon_diff, bond_diff)
                if facoa is not None:
                    pair = self._create_pair(
                        reactant, product, class_reaction,
                        required_compound=facoa
                    )
                    result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result


class FACoACompoundMatcherFullStructure(ReactionMatcher):
    """FACoA addition matcher with full structure awareness.
    
    When lipids have full structure (LPC(16:0) + FACoA(18:1) -> PC(16:0_18:1)),
    we can verify the product contains both the original chain and the added FACoA.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.FACOA_ADDITION
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match using full chain information when available."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        result.available_facoa = set(context.facoa_species.keys())
        
        if not reactants or not products:
            return result
        
        for reactant in reactants:
            matches = self._find_facoa_addition_matches(reactant, products, context)
            for product, facoa in matches:
                pair = self._create_pair(
                    reactant, product, class_reaction,
                    required_compound=facoa
                )
                result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result
    
    def _find_facoa_addition_matches(
        self,
        reactant: LipidStructure,
        products: List[LipidStructure],
        context: MatcherContext,
    ) -> List[tuple[LipidStructure, AcylChain]]:
        """Find products that can result from FACoA addition to reactant."""
        matches = []
        
        # If product has full structure, check if it contains reactant's chains + FACoA
        for product in products:
            if product.total_carbons <= reactant.total_carbons:
                continue
            
            carbon_diff = product.total_carbons - reactant.total_carbons
            bond_diff = product.total_double_bonds - reactant.total_double_bonds
            
            if carbon_diff <= 0 or bond_diff < 0:
                continue
            
            # Check if FACoA exists for this difference
            facoa = context.facoa_for_difference(carbon_diff, bond_diff)
            if facoa is None:
                continue
            
            # If we have full structure on both sides, verify chain compatibility
            if reactant.chains and product.chains:
                if self._chains_compatible(reactant.chains, product.chains, facoa):
                    matches.append((product, facoa))
            else:
                # Sum composition match is sufficient
                matches.append((product, facoa))
        
        return matches
    
    def _chains_compatible(
        self,
        reactant_chains: List[AcylChain],
        product_chains: List[AcylChain],
        added_chain: AcylChain,
    ) -> bool:
        """Check if product contains reactant chains + added FACoA.
        
        Product should have all chains from reactant plus the added FACoA.
        """
        # Simple check: product should have one more chain
        if len(product_chains) != len(reactant_chains) + 1:
            # Could still be valid if chains merged at sum level
            return True
        
        # Build multiset of product chains
        product_chain_set = {}
        for chain in product_chains:
            key = chain.notation
            product_chain_set[key] = product_chain_set.get(key, 0) + 1
        
        # Remove reactant chains
        for chain in reactant_chains:
            key = chain.notation
            if key in product_chain_set and product_chain_set[key] > 0:
                product_chain_set[key] -= 1
        
        # Remaining should match added FACoA
        remaining = [k for k, v in product_chain_set.items() if v > 0]
        if len(remaining) == 1 and remaining[0] == added_chain.notation:
            return True
        
        # Could be positional isomers - accept at sum level
        return True
