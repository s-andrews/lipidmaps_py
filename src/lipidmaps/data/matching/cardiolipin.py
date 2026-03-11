"""
Cardiolipin reaction matcher.

Handles the special case of 2 PG -> 1 CL reactions, where two
phosphatidylglycerol molecules combine to form one cardiolipin.

This uses "all-to-all" matching because there's no way to determine
which specific PG pair formed which CL from sum composition alone.
"""

from itertools import combinations
from typing import Dict, List, Set, Tuple

from ..models.species_reaction import (
    ClassReaction,
    ReactionMatchResult,
    ReactionType,
    SpeciesReactionPair,
)
from ..utils.chain_parser import LipidStructure
from .base import MatcherContext, ReactionMatcher


class CardiolipinMatcher(ReactionMatcher):
    """Matches PG -> CL and CL -> PG reactions.
    
    Cardiolipin synthesis: 2 PG -> CL
    - CL has 4 acyl chains from 2 PG molecules
    - At sum level, we can't know which PGs combined
    - Matching every PG to every CL - all to all
    
    Cardiolipin degradation: CL -> PG
    - Similar issue: which PG comes from which CL?
    - Uses all-to-all matching
    
    This is necessary because:
    - CL(72:4) could come from PG(36:2) + PG(36:2)
    - Or from PG(34:1) + PG(38:3)
    - We have no way to distinguish without molecular species level data
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.CARDIOLIPIN
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match all PG-CL pairs (all-to-all).
        
        Unlike other reactions, we don't filter by composition
        because we can't determine the correct pairing at sum level.
        """
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        
        if not reactants or not products:
            return result
        
        # All-to-all matching
        for reactant in reactants:
            for product in products:
                # For CL reactions, create pair without composition validation
                # because we can't determine the correct pairing
                pair = self._create_pair(reactant, product, class_reaction)
                pair.validation_notes.append(
                    "All-to-all match: CL composition cannot be resolved to specific PG pairs"
                )
                result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result


class CardiolipinCompositionMatcher(ReactionMatcher):
    """Alternative CL matcher that attempts composition-aware matching.
    
    This provides a more constrained matching when we want to filter
    impossible combinations, while still accepting ambiguous ones.
    
    For 2 PG -> CL:
    - CL carbons should equal sum of 2 PG carbons
    - CL bonds should equal sum of 2 PG bonds
    
    This is still ambiguous but reduces impossible combinations.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.CARDIOLIPIN
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match PG-CL pairs with composition constraints."""
        
        # Determine direction
        if "CL" in class_reaction.product_class:
            return self._match_pg_to_cl(class_reaction, context)
        else:
            return self._match_cl_to_pg(class_reaction, context)
    
    def _match_pg_to_cl(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match 2 PG -> 1 CL with composition awareness."""
        
        pg_species = context.get_species(class_reaction.reactant_class)
        cl_species = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, pg_species, cl_species)
        
        if not pg_species or not cl_species:
            return result
        
        # Index CL by composition
        cl_by_composition: Dict[Tuple[int, int], List[LipidStructure]] = {}
        for cl in cl_species:
            key = (cl.total_carbons, cl.total_double_bonds)
            if key not in cl_by_composition:
                cl_by_composition[key] = []
            cl_by_composition[key].append(cl)
        
        # For each PG, find CL that could result from 2x this PG
        # or this PG + another PG
        seen_pairs: Set[str] = set()
        
        for pg in pg_species:
            # Case 1: 2x same PG
            double_key = (pg.total_carbons * 2, pg.total_double_bonds * 2)
            for cl in cl_by_composition.get(double_key, []):
                pair_id = f"{pg.full_name}->{cl.full_name}"
                if pair_id not in seen_pairs:
                    pair = self._create_pair(pg, cl, class_reaction)
                    pair.validation_notes.append("Could form from 2x this PG")
                    result.pairs.append(pair)
                    seen_pairs.add(pair_id)
            
            # Case 2: This PG + any other PG
            for other_pg in pg_species:
                combined_key = (
                    pg.total_carbons + other_pg.total_carbons,
                    pg.total_double_bonds + other_pg.total_double_bonds
                )
                for cl in cl_by_composition.get(combined_key, []):
                    pair_id = f"{pg.full_name}->{cl.full_name}"
                    if pair_id not in seen_pairs:
                        pair = self._create_pair(pg, cl, class_reaction)
                        pair.validation_notes.append(
                            f"Could form with PG partner: {other_pg.full_name}"
                        )
                        result.pairs.append(pair)
                        seen_pairs.add(pair_id)
        
        result.pairs_matched = len(result.pairs)
        return result
    
    def _match_cl_to_pg(
        self,
        class_reaction: ClassReaction, 
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match CL -> PG (degradation) with composition awareness."""
        
        cl_species = context.get_species(class_reaction.reactant_class)
        pg_species = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, cl_species, pg_species)
        
        if not cl_species or not pg_species:
            return result
        
        # Index PG by composition
        pg_by_composition: Dict[Tuple[int, int], List[LipidStructure]] = {}
        for pg in pg_species:
            key = (pg.total_carbons, pg.total_double_bonds)
            if key not in pg_by_composition:
                pg_by_composition[key] = []
            pg_by_composition[key].append(pg)
        
        seen_pairs: Set[str] = set()
        
        for cl in cl_species:
            # Find PGs that could result from CL degradation
            # CL splits into 2 PGs, so look for PGs where 2x composition = CL
            
            # Exact half
            if cl.total_carbons % 2 == 0 and cl.total_double_bonds % 2 == 0:
                half_key = (cl.total_carbons // 2, cl.total_double_bonds // 2)
                for pg in pg_by_composition.get(half_key, []):
                    pair_id = f"{cl.full_name}->{pg.full_name}"
                    if pair_id not in seen_pairs:
                        pair = self._create_pair(cl, pg, class_reaction)
                        result.pairs.append(pair)
                        seen_pairs.add(pair_id)
            
            # Asymmetric split: any PG that could be one of the two products
            for (c, b), pgs in pg_by_composition.items():
                other_c = cl.total_carbons - c
                other_b = cl.total_double_bonds - b
                
                # Check if other half exists as a valid PG composition
                if other_c > 0 and other_b >= 0:
                    other_key = (other_c, other_b)
                    if other_key in pg_by_composition:
                        for pg in pgs:
                            pair_id = f"{cl.full_name}->{pg.full_name}"
                            if pair_id not in seen_pairs:
                                pair = self._create_pair(cl, pg, class_reaction)
                                result.pairs.append(pair)
                                seen_pairs.add(pair_id)
        
        result.pairs_matched = len(result.pairs)
        return result
